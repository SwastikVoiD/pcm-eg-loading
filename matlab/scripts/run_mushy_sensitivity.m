function T = run_mushy_sensitivity(cfg)
%RUN_MUSHY_SENSITIVITY  Effect of the mushy-band width (0.5, 1, 2 K).
%   Both PCMs, low route, base h; constant source plus sine and square waves
%   (A = 20 K, P = 30 min, four start phases); windows of 30 and 60 min.
%   Output: results/mushy_sensitivity.csv

ph = cfg.phases;
cases = [0 0 3600 0; repmat([1 20 1800], numel(ph), 1), ph(:); repmat([2 20 1800], numel(ph), 1), ph(:)];
groups = {'sine', 1 + (1:numel(ph)); 'square', 1 + numel(ph) + (1:numel(ph))};
grid = make_grid(cfg.geom, cfg.r_in, cfg.r_out, cfg.N);
w = cfg.w(:); nc = size(cases, 1);
[ww, cc] = ndgrid(w, 1:nc); ww = ww(:); cc = cc(:);
rows = {};
for ip = 1:numel(cfg.pcms)
    for dTm = [0.5 1.0 2.0]
        pcm = pcm_properties(cfg.pcms{ip});
        pcm.dTm = dTm;
        props = composite_properties(pcm, ww, k_scenario_beta('low'));
        src = struct('Tmean', pcm.Tm + cfg.dT_super, 'shape_id', cases(cc,1), 'A', cases(cc,2), ...
                     'P', cases(cc,3), 'phase', cases(cc,4));
        res = solve_pcm_1d(props, grid, cfg.h_list(1), pcm.Tm - cfg.dT_sub, src, 3600, cfg.dt, 60);
        for t_c = [1800 3600]
            E = energy_at_window(res, t_c);
            Ec = phase_mean(E, ww, cc, 1);
            wc = optimum_loading(w, Ec);
            for g = 1:size(groups, 1)
                Ef = phase_mean(E, ww, cc, groups{g, 2});
                [wf, Efs] = optimum_loading(w, Ef);
                pen = 100 * (1 - interp1(w, Ef, wc, 'pchip') / Efs);
                rows(end+1, :) = {string(pcm.name), dTm, t_c/60, string(groups{g,1}), wc, wf, pen}; %#ok<AGROW>
            end
        end
    end
end
T = cell2table(rows, 'VariableNames', {'pcm', 'dTm_K', 't_c_min', 'shape', 'w_const', 'w_fluct', 'penalty_pct'});
disp(T);
writetable(T, fullfile(cfg.results_dir, 'mushy_sensitivity.csv'));
end

function Em = phase_mean(E, ww, cc, idx)
Em = 0;
for q = idx(:)'
    sel = cc == q;
    [~, o] = sort(ww(sel)); e = E(sel);
    Em = Em + e(o);
end
Em = Em / numel(idx);
end
