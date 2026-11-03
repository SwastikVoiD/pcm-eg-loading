function T = run_independence(cfg)
%RUN_INDEPENDENCE  Grid and time-step independence on the production geometry.
%   Erythritol, mid-k scenario, sine source (A = 20 K, P = 30 min), h = 1000,
%   loadings 0, 5, 10, 20 wt%, stored energy after 60 min.

pcm = pcm_properties('erythritol');
w = [0 0.05 0.10 0.20];
props = composite_properties(pcm, w, k_scenario_beta('mid'));
B = numel(w);
src = struct('Tmean', pcm.Tm + cfg.dT_super, 'A', 20*ones(B,1), 'P', 1800*ones(B,1), ...
             'shape_id', ones(B,1), 'phase', zeros(B,1));
runs = [40 2; 80 2; 160 2; 80 1; 80 5];
E = zeros(size(runs,1), B);
for r = 1:size(runs, 1)
    grid = make_grid(cfg.geom, cfg.r_in, cfg.r_out, runs(r,1));
    res = solve_pcm_1d(props, grid, cfg.h_list(1), pcm.Tm - cfg.dT_sub, src, 3600, runs(r,2), 60);
    E(r, :) = res.E(:, end)';
end
ref = E(3, :);                                % N = 160, dt = 2 s
dev_pct = 100 * max(abs(E - ref) ./ ref, [], 2);
T = table(runs(:,1), runs(:,2), E/1e3, dev_pct, ...
    'VariableNames', {'N', 'dt_s', 'E60_kJ_per_m', 'max_dev_pct_vs_N160'});
disp(T);
writetable(T, fullfile(cfg.results_dir, 'independence.csv'));
end
