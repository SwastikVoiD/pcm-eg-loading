function R = run_li2020_comparison(cfg)
%RUN_LI2020_COMPARISON  Calibrated consistency check against Li et al. (2020).
%   Li et al. modelled the same annulus and salt PCM (Table 1) with 2D CFD
%   including natural convection, heated by a gas HTF whose heat transfer
%   coefficient is not reported. Step 1 calibrates h so that this model matches
%   their constant-source melting time (49.7 min). Step 2 predicts the melting
%   times for their fluctuating sources with no further tuning and compares.
%   Li et al. values are read from their text (Sections 4.1-4.2).

pcm  = pcm_properties('li2020_salt');
grid = make_grid('radial', 0.0125, 0.025, cfg.N);
T0 = pcm.Tm - 10; Tmean = 300;                % 573 K

% ---- step 1: calibrate h on the constant source
hs = logspace(log10(20), log10(400), 40)'; B = numel(hs);
props = composite_properties(pcm, zeros(B,1), 0);
src = struct('Tmean', Tmean, 'A', zeros(B,1), 'P', 3600*ones(B,1), 'shape_id', zeros(B,1), 'phase', zeros(B,1));
res = solve_pcm_1d(props, grid, hs, T0, src, 3*3600, cfg.dt, 10);
tm = nan(B,1);
for b = 1:B
    j = find(res.fl(b,:) >= 0.999, 1, 'first');
    if ~isempty(j), tm(b) = res.t(j)/60; end
end
ok = ~isnan(tm);
h_cal = interp1(tm(ok), hs(ok), 49.7);
fprintf('Calibrated h = %.1f W/m2K\n', h_cal);

% ---- step 2: predict the fluctuating cases
names = ["constant"; "P=1 min"; "P=6 min"; "P=20 min"; "P=60 min"; "A=50 K"; "A=150 K"];
A = [0 100 100 100 100 50 150]'; P = [3600 60 360 1200 3600 3600 3600]';
li_t = [49.7 49.8 49.3 NaN 37.3 44.2 32.3]';
B = numel(A);
props = composite_properties(pcm, zeros(B,1), 0);
src = struct('Tmean', Tmean, 'A', A, 'P', P, 'shape_id', [0; ones(B-1,1)], 'phase', zeros(B,1));
res = solve_pcm_1d(props, grid, h_cal, T0, src, 3*3600, 1.0, 5);
t_melt = nan(B,1); E_melt = nan(B,1);
for b = 1:B
    j = find(res.fl(b,:) >= 0.999, 1, 'first');
    t_melt(b) = res.t(j)/60; E_melt(b) = res.E(b,j)/1e3;
end
R = table(names, A, P/60, li_t, t_melt, 100*(t_melt - li_t)./li_t, E_melt, ...
    'VariableNames', {'case', 'A_K', 'P_min', 'Li2020_t_melt_min', 'model_t_melt_min', 'diff_pct', 'model_E_melt_kJ_per_m'});
R.h_cal = repmat(h_cal, B, 1);
disp(R);
writetable(R, fullfile(cfg.results_dir, 'li2020_comparison.csv'));
end
