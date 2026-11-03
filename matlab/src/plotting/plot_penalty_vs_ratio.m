function fig = plot_penalty_vs_ratio(cfg, ALL)
%PLOT_PENALTY_VS_RATIO  Design penalty of the steady-source optimum and change
%   in achievable energy vs the source-period / charging-window ratio, for every
%   fluctuating case at the base h without bypass.
T = ALL(ALL.shape ~= "constant" & ~ALL.bypass & ALL.h == cfg.h_list(1), :);
ratio = T.P_min ./ T.t_c_min;
fig = figure('Color', 'w', 'Position', [100 100 900 360]);
pcms = unique(T.pcm); ks = {'low', 'mid', 'high'}; mk = {'o', 's'}; cols = lines(3);
for panel = 1:2
    subplot(1, 2, panel); hold on; box on;
    for ip = 1:numel(pcms)
        for ik = 1:3
            r = T.pcm == pcms(ip) & T.k_scenario == ks{ik};
            x = ratio(r) .* exp(0.16*(rand(nnz(r),1) - 0.5));     % jitter for visibility
            if panel == 1, y = T.penalty_pct(r); else, y = T.dE_pct(r); end
            scatter(x, y, 14, cols(ik,:), mk{ip}, 'filled', 'MarkerFaceAlpha', 0.6, ...
                'DisplayName', sprintf('%s, %s', pcms(ip), ks{ik}));
        end
    end
    set(gca, 'XScale', 'log'); xline(1, 'k:', 'HandleVisibility', 'off');
    xticks([1/12 1/4 1/2 1 2]); xticklabels({'1/12', '1/4', '1/2', '1', '2'});
    xlabel('Source period / charging window, P/t_c');
    if panel == 1
        ylabel('Design penalty of steady-source w^* (%)'); legend('Location', 'northwest', 'FontSize', 7);
    else
        ylabel('Change in achievable energy (%)');
    end
end
save_fig(fig, cfg, 'fig6_penalty_vs_ratio');
end
