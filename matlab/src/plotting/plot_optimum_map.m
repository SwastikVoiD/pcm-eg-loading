function fig = plot_optimum_map(cfg, ALL)
%PLOT_OPTIMUM_MAP  Optimum loading vs charging window, constant source, with
%   the spread across all fluctuating sources shown as error bars.
fig = figure('Color', 'w', 'Position', [100 100 900 380]);
for ip = 1:numel(cfg.pcms)
    subplot(1, numel(cfg.pcms), ip); hold on; box on;
    for ik = 1:numel(cfg.k_scenarios)
        for h = cfg.h_list(1)
            rows = ALL.pcm == cfg.pcms{ip} & ALL.k_scenario == cfg.k_scenarios{ik} & ALL.h == h;
            T = ALL(rows, :);
            tc = unique(T.t_c_min);
            [wc, lo, hi] = deal(zeros(size(tc)));
            for i = 1:numel(tc)
                r = T.t_c_min == tc(i);
                wc(i) = T.w_opt(r & T.shape == "constant");
                lo(i) = min(T.w_opt(r)); hi(i) = max(T.w_opt(r));
            end
            errorbar(tc, 100*wc, 100*(wc - lo), 100*(hi - wc), '-o', 'LineWidth', 1.3, ...
                'DisplayName', cfg.k_scenarios{ik});
        end
    end
    yline(100*cfg.w_formstable, 'k:', 'HandleVisibility', 'off');
    xlabel('Charging window (min)'); ylabel('Optimum EG fraction w^* (%)');
    title(cfg.pcms{ip}, 'FontWeight', 'normal'); legend('Location', 'northeast');
end
save_fig(fig, cfg, 'fig5_optimum_map');
end
