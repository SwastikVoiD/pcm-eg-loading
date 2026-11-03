"""Mirror of summarize_sweep.m / optimum_loading.m"""
import numpy as np, pandas as pd
W_FS = 0.10

def optimum_loading(w, E, w_min=-np.inf):
    """Maximum of a shape-preserving (pchip) interpolant on a 1e-4 grid."""
    from scipy.interpolate import PchipInterpolator
    w = np.asarray(w, float); E = np.asarray(E, float)
    wf = np.arange(max(w[0], w_min), w[-1] + 1e-12, 1e-4)
    Ef = PchipInterpolator(w, E)(wf)
    i = int(np.argmax(Ef))
    return wf[i], Ef[i], 0 < i < len(wf) - 1

def summarize(d, t_c):
    j = int(np.argmin(np.abs(d['t'] - t_c)))
    E = d['E'][:, j]; fl = d['fl'][:, j]
    cases = d['cases']; w_all = d['w']; cidx_all = d['case']
    w = np.unique(w_all)
    keys = []
    for row in cases[:, :3]:
        if not any(np.allclose(row, k) for k in keys): keys.append(row)
    rows = []; Ebars = []
    for key in keys:
        cidx = [i for i in range(len(cases)) if np.allclose(cases[i, :3], key)]
        Ew = []; Fw = []; wq_list = []
        for q in cidx:
            sel = cidx_all == q; order = np.argsort(w_all[sel])
            e = E[sel][order]; f = fl[sel][order]
            Ew.append(e); Fw.append(f); wq_list.append(optimum_loading(w, e)[0])
        Ebar = np.mean(Ew, 0); Fbar = np.mean(Fw, 0); Ebars.append(Ebar)
        wo, Eo, inter = optimum_loading(w, Ebar)
        wfs = optimum_loading(w, Ebar, W_FS)[0]
        rows.append(dict(shape=['constant', 'sine', 'square'][int(key[0])], A_K=key[1], P_min=key[2] / 60,
                         w_opt=wo, w_opt_fs=wfs, w_opt_lo=min(wq_list), w_opt_hi=max(wq_list),
                         E_opt_kJ_per_m=Eo / 1e3, fl_opt=np.interp(wo, w, Fbar), interior=inter))
    from scipy.interpolate import PchipInterpolator
    w0 = rows[0]['w_opt']; E0 = rows[0]['E_opt_kJ_per_m']
    for r, Ebar in zip(rows, Ebars):
        r['dw'] = r['w_opt'] - w0
        r['penalty_pct'] = 100 * (1 - PchipInterpolator(w, Ebar)(w0) / 1e3 / r['E_opt_kJ_per_m'])
        r['dE_pct'] = 100 * (r['E_opt_kJ_per_m'] / E0 - 1)
    return pd.DataFrame(rows), w, Ebars
