function V = run_validation_neumann(cfg)
%RUN_VALIDATION_NEUMANN  Verify the solver against the two-phase Neumann solution.
%   Pure erythritol, planar semi-infinite slab (0.3 m domain, 600 cells), wall
%   held at Tm + 20 K via a very large h, initial Tm - 10 K. Compared at 30 and
%   120 min (Table III of the paper). V is a 1 x 2 struct array.

pcm = pcm_properties('erythritol');
props = composite_properties(pcm, 0, 0);
Tw = pcm.Tm + 20; Ti = pcm.Tm - 10;
alpha = pcm.k / (pcm.rho * pcm.cp);
lam = neumann_lambda(pcm.cp*(Tw - pcm.Tm)/pcm.L, pcm.cp*(pcm.Tm - Ti)/pcm.L);
grid = make_grid('planar', 0, 0.3, 600);
src = struct('Tmean', Tw, 'A', 0, 'P', 1, 'shape_id', 0, 'phase', 0);
near = grid.rc < 0.08;
times = [1800 7200];
for i = 1:numel(times)
    res = solve_pcm_1d(props, grid, 1e9, Ti, src, times(i), cfg.dt, 600);
    [Ta, s] = neumann_solution(grid.rc, times(i), Tw, Ti, pcm.Tm, alpha, lam);
    xf = melt_front_position(grid.rc, res.f(1, :));
    err = res.T(1, :) - Ta;
    V(i) = struct('t_min', times(i)/60, 'lambda', lam, 'front_num_mm', 1e3*xf, ...
        'front_exact_mm', 1e3*s, 'front_err_pct', 100*(xf - s)/s, ...
        'maxTerr_K', max(abs(err(near))), 'rmsTerr_K', sqrt(mean(err(near).^2)), ...
        'energy_balance_rel', (res.E(1,end) - res.Qin(1,end)) / res.Qin(1,end), ...
        'x', grid.rc, 'T_num', res.T(1, :), 'T_exact', Ta); %#ok<AGROW>
    fprintf('t = %3d min | front %.2f mm (exact %.2f mm, err %.2f%%) | max |dT| %.3f K | energy balance %.1e\n', ...
        times(i)/60, V(i).front_num_mm, V(i).front_exact_mm, V(i).front_err_pct, V(i).maxTerr_K, V(i).energy_balance_rel);
end
fprintf('lambda = %.4f\n', lam);
save(fullfile(cfg.results_dir, 'validation_neumann.mat'), 'V');
end
