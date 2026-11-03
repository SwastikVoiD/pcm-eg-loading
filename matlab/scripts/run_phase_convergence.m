function D = run_phase_convergence(cfg)
%RUN_PHASE_CONVERGENCE  Start-phase convergence (P = 60 min; t_c = 30 and 60 min).
%   Four start phases are converged when the charging window spans at least
%   one source period, but not when P = 60 min and t_c = 30 min (square waves
%   sampled at 0 and pi sit on their switching points). This reruns that regime
%   for every PCM and route with 16 phases and reports 4/8/16-phase results,
%   at t_c = 30 min (P > t_c) and t_c = 60 min (P = t_c, the P <= t_c regime).
%   Output: results/phase_convergence_Pgt_tc.csv

ph16 = (0:15) * pi/8;
cases = [0 0 3600 0];
for sid = [1 2]
    for A = cfg.amps
        cases = [cases; repmat([sid A 3600], 16, 1), ph16(:)]; %#ok<AGROW>
    end
end
grid = make_grid(cfg.geom, cfg.r_in, cfg.r_out, cfg.N);
w = cfg.w(:); nw = numel(w); nc = size(cases, 1);
rows = {};
for ip = 1:numel(cfg.pcms)
    pcm = pcm_properties(cfg.pcms{ip});
    for ik = 1:numel(cfg.k_scenarios)
        ks = cfg.k_scenarios{ik};
        [ww, cc] = ndgrid(w, 1:nc); ww = ww(:); cc = cc(:);
        props = composite_properties(pcm, ww, k_scenario_beta(ks));
        src = struct('Tmean', pcm.Tm + cfg.dT_super, 'shape_id', cases(cc,1), 'A', cases(cc,2), ...
                     'P', cases(cc,3), 'phase', cases(cc,4));
        res = solve_pcm_1d(props, grid, cfg.h_list(1), pcm.Tm - cfg.dT_sub, src, 3600, cfg.dt, 60);
        for t_c = [1800 3600]
            E = energy_at_window(res, t_c);
            curve = @(idx) mean(cell2mat(arrayfun(@(q) sortby(E(cc==q), ww(cc==q)), idx(:)', 'UniformOutput', false)), 2);
            Ec = curve(1);
            [wc, Ecs] = optimum_loading(w, Ec);
            g = 0;
            for sid = [1 2]
                for A = cfg.amps
                    i16 = 1 + 16*g + (1:16); g = g + 1;
                    for nph = [4 8 16]
                        idx = i16(1:16/nph:end);
                        Ef = curve(idx);
                        [wf, Efs] = optimum_loading(w, Ef);
                        pen = 100 * (1 - interp1(w, Ef, wc, 'pchip') / Efs);
                        rows(end+1, :) = {string(pcm.name), string(ks), t_c/60, sid, A, nph, wc, wf, wf - wc, pen, 100*(Efs/Ecs - 1)}; %#ok<AGROW>
                    end
                end
            end
        end
        fprintf('%s %s done\n', pcm.name, ks);
    end
end
D = cell2table(rows, 'VariableNames', {'pcm', 'k', 't_c_min', 'shape_id', 'A_K', 'n_phases', 'w_const', 'w_fluct', 'dw', 'penalty_pct', 'dE_pct'});
writetable(D, fullfile(cfg.results_dir, 'phase_convergence_Pgt_tc.csv'));
for t_c = [30 60]
    for nph = [4 8 16]
        s = D(D.n_phases == nph & D.t_c_min == t_c, :);
        fprintf('t_c = %d min, %2d phases: max penalty %.2f %%, max |dw| %.3f, min dE %.1f %%\n', ...
            t_c, nph, max(s.penalty_pct), max(abs(s.dw)), min(s.dE_pct));
    end
end
end

function e = sortby(e, w)
[~, o] = sort(w); e = e(o);
end
