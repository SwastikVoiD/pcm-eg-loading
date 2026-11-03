"""
pcm1d.py  -  Python reference implementation (line-for-line mirror of the MATLAB code).

1D radial (or planar) enthalpy-method finite-volume solver for a composite PCM
annulus heated at the inner radius by a heat-transfer fluid with a time-varying
temperature. Implicit Euler in time; Voller-Swaminathan (1991) source-based
linearisation of the liquid fraction; batched Thomas algorithm so many cases
(e.g. every filler loading) are solved simultaneously.
"""
import numpy as np
from scipy.special import erf, erfc
from scipy.optimize import brentq

# ----------------------------------------------------------------------------
# Material data
# ----------------------------------------------------------------------------
def pcm_properties(name):
    if name == 'erythritol':
        return dict(name='erythritol', Tm=117.0, L=344e3, cp=2070.0, rho=1390.0,
                    k=0.72, dTm=1.0)
    if name == 'paraffin':   # n-docosane
        return dict(name='paraffin', Tm=44.0, L=256e3, cp=1950.0, rho=800.0,
                    k=0.22, dTm=1.0)
    if name == 'li2020_salt':  # LiNO3-NaNO3-KCl, Li et al. (2020) Table 1
        return dict(name='li2020_salt', Tm=160.0, L=266e3, cp=1330.0, rho=2297.0,
                    k=0.88, dTm=1.0)
    raise ValueError(name)

EG = dict(rho=2200.0, cp=710.0)

# Liquid properties for the natural-convection model (Raithby-Hollands annulus
# correlation). Viscosity and expansion coefficient are typical literature
# values [VERIFY]; their uncertainty is bracketed with conv['ra_mult'].
LIQUID = {
    'erythritol': dict(k=0.326, rho=1300.0, cp=2760.0, mu=0.016, beta=5.0e-4),
    'paraffin':   dict(k=0.15,  rho=778.0,  cp=2380.0, mu=0.0036, beta=9.0e-4),
}

def convection_params(pcm_name, w, w_s=0.10, ra_mult=1.0):
    L = LIQUID[pcm_name]
    nu = L['mu'] / L['rho']; alpha = L['k'] / (L['rho'] * L['cp'])
    return dict(Pr=nu / alpha, Ra_coef=9.81 * L['beta'] / (nu * alpha) * ra_mult,
                w=np.atleast_1d(np.asarray(w, float)), w_s=w_s)

def convection_multiplier(T, f, grid, Tm, conv):
    """Effective-conductivity factor for the molten layer (B,) and liquid-cell count.
    Raithby & Hollands (1975) correlation for a horizontal concentric annulus,
    applied quasi-steadily to the melt layer between the tube and the melt front.
    EG suppresses convection: enhancement falls linearly to zero at w = w_s."""
    liq = f >= 0.5
    n_liq = np.cumprod(liq, axis=1).sum(axis=1)                # contiguous liquid cells from the wall
    ri = grid['rf'][0]; s = grid['rf'][n_liq]; delta = s - ri
    dT = np.maximum(T[:, 0] - Tm, 0.0)
    Di = 2 * ri; Do = 2 * s
    with np.errstate(divide='ignore', invalid='ignore'):
        RaL = conv['Ra_coef'] * dT * delta**3
        Rac = np.log(Do / Di)**4 / (delta**3 * (Di**-0.6 + Do**-0.6)**5) * RaL
        ratio = 0.386 * (conv['Pr'] / (0.861 + conv['Pr']))**0.25 * np.maximum(Rac, 0.0)**0.25
    ratio = np.where(n_liq > 0, np.maximum(np.nan_to_num(ratio, nan=1.0), 1.0), 1.0)
    phi = np.clip(1.0 - conv['w'] / conv['w_s'], 0.0, 1.0)
    return 1.0 + phi * (ratio - 1.0), n_liq

K_SCENARIOS = {'low': 6.38, 'mid': 28.28, 'high': 117.9}   # W/(m K) per unit mass fraction

def composite_properties(pcm, w, beta):
    w = np.atleast_1d(np.asarray(w, float))
    rho = 1.0 / ((1 - w) / pcm['rho'] + w / EG['rho'])
    cp = (1 - w) * pcm['cp'] + w * EG['cp']
    L = (1 - w) * pcm['L']
    k = pcm['k'] + beta * w
    return dict(rho=rho, cp=cp, L=L, k=k, Tm=pcm['Tm'], dTm=pcm['dTm'])

def k_maxwell(km, kf, phi):
    return km * (kf + 2 * km + 2 * phi * (kf - km)) / (kf + 2 * km - phi * (kf - km))

# ----------------------------------------------------------------------------
# Source profiles
# ----------------------------------------------------------------------------
def source_temperature(t, Tmean, A, P, shape_id, phase):
    """t scalar [s]; A, P, phase, shape_id arrays (B,).
    shape_id: 0 = constant, 1 = sine, 2 = square wave (50 % duty)."""
    s = np.sin(2 * np.pi * t / P + phase)
    sq = np.where(s >= 0, 1.0, -1.0)
    return Tmean + A * np.where(shape_id == 1, s, np.where(shape_id == 2, sq, 0.0))

# ----------------------------------------------------------------------------
# Grid
# ----------------------------------------------------------------------------
def make_grid(geom, r_in, r_out, N):
    rf = np.linspace(r_in, r_out, N + 1)
    rc = 0.5 * (rf[:-1] + rf[1:])
    dr = rf[1] - rf[0]
    if geom == 'radial':
        V = np.pi * (rf[1:] ** 2 - rf[:-1] ** 2)
        Af = 2 * np.pi * rf
    else:  # planar, unit area
        V = np.full(N, dr)
        Af = np.ones(N + 1)
    return dict(rf=rf, rc=rc, dr=dr, V=V, Af=Af, N=N, geom=geom)

# ----------------------------------------------------------------------------
# Batched Thomas algorithm: a (sub), b (diag), c (super), d (rhs), all (B,N)
# ----------------------------------------------------------------------------
def tdma_batch(a, b, c, d):
    B, N = b.shape
    cp_ = np.empty((B, N)); dp = np.empty((B, N))
    cp_[:, 0] = c[:, 0] / b[:, 0]
    dp[:, 0] = d[:, 0] / b[:, 0]
    for i in range(1, N):
        m = b[:, i] - a[:, i] * cp_[:, i - 1]
        cp_[:, i] = c[:, i] / m
        dp[:, i] = (d[:, i] - a[:, i] * dp[:, i - 1]) / m
    x = np.empty((B, N))
    x[:, -1] = dp[:, -1]
    for i in range(N - 2, -1, -1):
        x[:, i] = dp[:, i] - cp_[:, i] * x[:, i + 1]
    return x

# ----------------------------------------------------------------------------
# Core solver
# ----------------------------------------------------------------------------
def solve_pcm(props, grid, h, T0, src, t_end, dt, record_every=60.0,
              tol=1e-9, max_iter=100, bypass=False, conv=None):
    """
    props : dict of (B,) arrays rho, cp, L, k and scalars Tm, dTm
    src   : dict Tmean, A(B,), P(B,), shape, phase(B,)
    Returns dict with time vector and per-case histories.
    """
    rho = props['rho'][:, None]; cp = props['cp'][:, None]
    L = props['L'][:, None]; k = props['k'][:, None]
    B = rho.shape[0]; N = grid['N']
    Ts = props['Tm'] - props['dTm'] / 2; Tl = props['Tm'] + props['dTm'] / 2
    dTm = props['dTm']
    V = grid['V'][None, :]; dr = grid['dr']; Af = grid['Af']

    Gint = Af[None, 1:-1] * k / dr                  # (B,N-1) interior faces
    Gb = Af[0] / (1.0 / np.asarray(h, float) + (dr / 2) / k[:, 0])   # (B,) boundary conductance; h scalar or (B,)
    C = rho * cp * V / dt
    Lat = rho * L * V / dt

    T = np.full((B, N), float(T0))
    f = np.clip((T - Ts) / dTm, 0, 1)

    nsteps = int(round(t_end / dt))
    rec_stride = max(1, int(round(record_every / dt)))
    times, E_hist, fl_hist, Qin_hist, Thtf_hist = [0.0], [], [], [np.zeros(B)], []
    E0 = stored_energy(T, f, props, grid, T0)
    E_hist.append(E0); fl_hist.append(melt_fraction(f, grid))
    Qin = np.zeros(B)
    Thtf_hist.append(source_temperature(0.0, src['Tmean'], src['A'], src['P'], src['shape_id'], src['phase']))
    iters_total = 0; max_it_used = 0

    for n in range(1, nsteps + 1):
        t = n * dt
        Thtf = source_temperature(t, src['Tmean'], src['A'], src['P'], src['shape_id'], src['phase'])
        Told, fold = T.copy(), f.copy()
        Tk, fk = T.copy(), f.copy()
        if conv is not None:   # lagged natural-convection enhancement of the melt layer
            eps, n_liq = convection_multiplier(Told, fold, grid, props['Tm'], conv)
            j = np.arange(N - 1)[None, :]
            M = np.where(j < (n_liq[:, None] - 1), eps[:, None], 1.0)
            Gint_s = Gint * M
            m0 = np.where(n_liq >= 1, eps, 1.0)
            Gb_s = Af[0] / (1.0 / np.asarray(h, float) + (dr / 2) / (k[:, 0] * m0))
        else:
            Gint_s, Gb_s = Gint, Gb
        # bypass rule: HTF flow stops while the source is colder than the store wall
        Gb_eff = Gb_s * (Thtf > Told[:, 0]) if bypass else Gb_s
        for it in range(max_iter):
            S = np.where(((fk > 0) & (fk < 1)) | ((fk <= 0) & (Tk > Ts)) | ((fk >= 1) & (Tk < Tl)),
                         1.0 / dTm, 0.0)
            Tfk = Ts + fk * dTm
            a = np.zeros((B, N)); c = np.zeros((B, N))
            a[:, 1:] = -Gint_s; c[:, :-1] = -Gint_s
            diag = C + Lat * S
            diag[:, 1:] += Gint_s; diag[:, :-1] += Gint_s
            diag[:, 0] += Gb_eff
            rhs = C * Told + Lat * (fold - fk + S * Tfk)
            rhs[:, 0] += Gb_eff * Thtf
            Tn = tdma_batch(a, diag, c, rhs)
            fn = np.clip(fk + S * (Tn - Tfk), 0.0, 1.0)
            dmax = np.max(np.abs(fn - fk)); tmax = np.max(np.abs(Tn - Tk))
            Tk, fk = Tn, fn
            if dmax < tol and tmax < 1e-7:
                break
        iters_total += it + 1; max_it_used = max(max_it_used, it + 1)
        T, f = Tk, fk
        Qin += Gb_eff * (Thtf - T[:, 0]) * dt
        if n % rec_stride == 0:
            times.append(t)
            E_hist.append(stored_energy(T, f, props, grid, T0))
            fl_hist.append(melt_fraction(f, grid))
            Qin_hist.append(Qin.copy())
            Thtf_hist.append(Thtf)
    E = np.array(E_hist).T - E0[:, None]
    return dict(t=np.array(times), E=E, fl=np.array(fl_hist).T, Qin=np.array(Qin_hist).T,
                Thtf=np.array(Thtf_hist).T, T=T, f=f,
                mean_iter=iters_total / nsteps, max_iter=max_it_used)

def stored_energy(T, f, props, grid, T0):
    rho = props['rho'][:, None]; cp = props['cp'][:, None]; L = props['L'][:, None]
    return np.sum(rho * grid['V'][None, :] * (cp * (T - T0) + f * L), axis=1)

def melt_fraction(f, grid):
    return np.sum(f * grid['V'][None, :], axis=1) / np.sum(grid['V'])

# ----------------------------------------------------------------------------
# Neumann two-phase analytical solution (equal properties in both phases)
# ----------------------------------------------------------------------------
def neumann_lambda(St_l, St_s):
    g = lambda lam: (St_l * np.exp(-lam**2) / erf(lam) - St_s * np.exp(-lam**2) / erfc(lam)
                     - lam * np.sqrt(np.pi))
    return brentq(g, 1e-6, 5.0)

def neumann_solution(x, t, Tw, Ti, Tm, alpha, lam):
    eta = x / (2 * np.sqrt(alpha * t))
    s = 2 * lam * np.sqrt(alpha * t)
    T = np.where(x <= s, Tw - (Tw - Tm) * erf(eta) / erf(lam),
                 Ti + (Tm - Ti) * erfc(eta) / erfc(lam))
    return T, s
