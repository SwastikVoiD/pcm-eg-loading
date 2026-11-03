function out = solve_pcm_1d(props, grid, h, T0, src, t_end, dt, record_every, opts)
%SOLVE_PCM_1D  Implicit enthalpy-method solver for B composite-PCM cases at once.
%
%   out = solve_pcm_1d(props, grid, h, T0, src, t_end, dt, record_every)
%   out = solve_pcm_1d(..., record_every, opts)   opts.bypass = true stops the
%         HTF flow (h -> 0) during any step where T_htf is below the wall-cell
%         temperature, so a cold source cannot pull heat back out of the store.
%         opts.conv = convection_params(...) adds the natural-convection
%         enhancement of the melt layer (convection_multiplier.m), lagged one step.
%
%   Governing equation (per unit tube length, radial or planar):
%       rho*cp*dT/dt + rho*L*df/dt = (1/r) d/dr( r k dT/dr )
%   Inner wall : -k dT/dr = h (T_htf(t) - T)   (HTF-side resistance in series
%                with half a cell of conduction)
%   Outer wall : adiabatic
%   Liquid fraction: f = (T - Ts)/dTm clipped to [0,1] over a mushy band dTm.
%
%   Time: implicit Euler. The latent term is linearised with the source-based
%   scheme of Voller & Swaminathan (1991):  f = f_k + S (T - F^-1(f_k)),
%   S = 1/dTm for cells that are (or are crossing into) the mushy band, else 0.
%   The resulting tridiagonal systems are solved with tdma_batch and iterated
%   until f and T stop changing.
%
%   Inputs
%     props : struct with B x 1 vectors rho, cp, L, k and scalars Tm, dTm
%     grid  : from make_grid
%     h     : HTF-side heat transfer coefficient [W/m2K], scalar or B x 1
%     T0    : initial (and energy reference) temperature [degC]
%     src   : struct Tmean (scalar) and B x 1 A, P, shape_id, phase
%     t_end, dt, record_every : [s]
%   Output fields
%     t (1 x nrec), E, fl, Qin, Thtf (B x nrec), T, f (final fields),
%     mean_iter, max_iter

if nargin < 8, record_every = 60; end
if nargin < 9, opts = struct(); end
bypass = isfield(opts, 'bypass') && opts.bypass;
useconv = isfield(opts, 'conv') && ~isempty(opts.conv);
tol_f = 1e-9; tol_T = 1e-7; max_iter = 100;

rho = props.rho(:); cp = props.cp(:); L = props.L(:); k = props.k(:);
B = numel(rho); N = grid.N;
dTm = props.dTm; Ts = props.Tm - dTm/2; Tl = props.Tm + dTm/2;
V = grid.V(:)'; Af = grid.Af(:)'; dr = grid.dr;

% conductances and capacities
Gint = k .* (Af(2:end-1) / dr);             % B x (N-1), interior faces
Gb   = Af(1) ./ (1 ./ h(:) + (dr/2) ./ k);   % B x 1, inner boundary (h scalar or B x 1)
C    = (rho .* cp) .* V / dt;                % B x N
Lat  = (rho .* L)  .* V / dt;                % B x N

% constant part of the tridiagonal matrix
a = zeros(B, N); c = zeros(B, N);
a(:, 2:end)   = -Gint;
c(:, 1:end-1) = -Gint;
diag0 = C;
diag0(:, 2:end)   = diag0(:, 2:end)   + Gint;
diag0(:, 1:end-1) = diag0(:, 1:end-1) + Gint;   % boundary term added per step

% initial state
T = T0 * ones(B, N);
f = min(max((T - Ts) / dTm, 0), 1);

nsteps = round(t_end / dt);
stride = max(1, round(record_every / dt));
nrec   = floor(nsteps / stride) + 1;
out.t    = zeros(1, nrec);
out.E    = zeros(B, nrec);
out.fl   = zeros(B, nrec);
out.Qin  = zeros(B, nrec);
out.Thtf = zeros(B, nrec);

E0 = stored_energy(T, f, props, grid, T0);
out.fl(:, 1)   = melt_fraction(f, grid);
out.Thtf(:, 1) = source_temperature(0, src);
Qin = zeros(B, 1);
irec = 1; it_total = 0; it_max = 0;

for n = 1:nsteps
    t = n * dt;
    Thtf = source_temperature(t, src);
    Told = T; fold = f; Tk = T; fk = f;
    if useconv
        [epsc, nliq] = convection_multiplier(Told, fold, grid, props.Tm, opts.conv);
        mask = (1:N-1) < nliq;                      % faces with liquid on both sides
        M = ones(B, N-1); Erep = repmat(epsc, 1, N-1); M(mask) = Erep(mask);
        Gint_s = Gint .* M;
        m0 = ones(B, 1); m0(nliq >= 1) = epsc(nliq >= 1);
        Gb_s = Af(1) ./ (1 ./ h(:) + (dr/2) ./ (k .* m0));
        a = zeros(B, N); c = zeros(B, N);
        a(:, 2:end) = -Gint_s; c(:, 1:end-1) = -Gint_s;
        diag0 = C;
        diag0(:, 2:end)   = diag0(:, 2:end)   + Gint_s;
        diag0(:, 1:end-1) = diag0(:, 1:end-1) + Gint_s;
    else
        Gb_s = Gb;
    end
    if bypass
        Gb_eff = Gb_s .* (Thtf > Told(:, 1));
    else
        Gb_eff = Gb_s;
    end
    for it = 1:max_iter
        S = double(((fk > 0) & (fk < 1)) | ((fk <= 0) & (Tk > Ts)) | ((fk >= 1) & (Tk < Tl))) / dTm;
        Tfk = Ts + fk * dTm;
        dg  = diag0 + Lat .* S;
        dg(:, 1) = dg(:, 1) + Gb_eff;
        rhs = C .* Told + Lat .* (fold - fk + S .* Tfk);
        rhs(:, 1) = rhs(:, 1) + Gb_eff .* Thtf;
        Tn = tdma_batch(a, dg, c, rhs);
        fn = min(max(fk + S .* (Tn - Tfk), 0), 1);
        dmax = max(abs(fn(:) - fk(:)));
        tmax = max(abs(Tn(:) - Tk(:)));
        Tk = Tn; fk = fn;
        if dmax < tol_f && tmax < tol_T
            break
        end
    end
    it_total = it_total + it;
    it_max = max(it_max, it);
    T = Tk; f = fk;
    Qin = Qin + Gb_eff .* (Thtf - T(:, 1)) * dt;

    if mod(n, stride) == 0
        irec = irec + 1;
        out.t(irec)       = t;
        out.E(:, irec)    = stored_energy(T, f, props, grid, T0) - E0;
        out.fl(:, irec)   = melt_fraction(f, grid);
        out.Qin(:, irec)  = Qin;
        out.Thtf(:, irec) = Thtf;
    end
end
out.T = T; out.f = f;
out.mean_iter = it_total / nsteps;
out.max_iter  = it_max;
end
