function fig = plot_bypass(cfg, ALL)
%PLOT_BYPASS  Bypass rule vs HTF always flowing, 60 min window, every fluctuating
%   case (both PCMs) at the base h. Left: change in achievable energy against the
%   constant source. Right: design penalty of the steady-source optimum.
T = ALL(ALL.h == cfg.h_list(1) & ALL.t_c_min == 60 & ALL.shape ~= "constant", :);
ks = {'low', 'mid', 'high'};
fig = figure('Color', 'w', 'Position', [100 100 900 380]);
cols = lines(2); fld = {'dE_pct', 'penalty_pct'};
ylab = {'Change in achievable energy vs constant source (%)', 'Design penalty of steady-source w^* (%)'};
for panel = 1:2
    subplot(1, 2, panel); hold on; box on;
    for ik = 1:3
        a = T(T.k_scenario == ks{ik} & ~T.bypass, :);
        b = T(T.k_scenario == ks{ik} & T.bypass, :);
        if isempty(b), continue; end
        plot(ik - 0.15 + 0.1*(rand(height(a),1) - 0.5), a.(fld{panel}), 'o', 'Color', cols(1,:), ...
            'DisplayName', 'HTF always flowing', 'HandleVisibility', vis(ik));
        plot(ik + 0.15 + 0.1*(rand(height(b),1) - 0.5), b.(fld{panel}), 's', 'Color', cols(2,:), ...
            'DisplayName', 'Bypass when cold', 'HandleVisibility', vis(ik));
    end
    yline(0, 'k:', 'HandleVisibility', 'off'); xlim([0.5 3.5]); xticks(1:3); xticklabels(ks);
    xlabel('Conductivity route'); ylabel(ylab{panel});
    if panel == 1, legend('Location', 'northwest', 'FontSize', 8); end
end
save_fig(fig, cfg, 'fig7_bypass');
end

function v = vis(ik)
if ik == 1, v = 'on'; else, v = 'off'; end
end
