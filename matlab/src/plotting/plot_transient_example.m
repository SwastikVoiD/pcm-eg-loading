function fig = plot_transient_example(cfg)
%PLOT_TRANSIENT_EXAMPLE  Source profiles, melt-fraction and stored-energy response.
%   Erythritol, 5 wt% EG, low-k route, h = base, P = 60 min, A = 20 K.
%   Dashed vertical lines mark the 30 and 60 min charging windows.
pcm = pcm_properties('erythritol');
props = composite_properties(pcm, [0.05 0.05 0.05], k_scenario_beta('low'));
src = struct('Tmean', pcm.Tm + cfg.dT_super, 'A', [0; 20; 20], 'P', [3600; 3600; 3600], ...
             'shape_id', [0; 1; 2], 'phase', [0; 0; 0]);
grid = make_grid(cfg.geom, cfg.r_in, cfg.r_out, cfg.N);
res = solve_pcm_1d(props, grid, cfg.h_list(1), pcm.Tm - cfg.dT_sub, src, cfg.t_end, cfg.dt, 30);
fig = figure('Color', 'w', 'Position', [100 100 560 700]);
names = {'constant', 'sine', 'square'};
ax1 = subplot(3,1,1); hold on; box on;
for b = 1:3, plot(res.t/60, res.Thtf(b,:), 'LineWidth', 1.3, 'DisplayName', names{b}); end
yline(pcm.Tm, 'k:', 'T_m', 'HandleVisibility', 'off'); ylabel('HTF temperature (^\circC)');
legend('Location', 'southwest');
ax2 = subplot(3,1,2); hold on; box on;
for b = 1:3, plot(res.t/60, res.fl(b,:), 'LineWidth', 1.3, 'HandleVisibility', 'off'); end
ylabel('Melt fraction (-)');
ax3 = subplot(3,1,3); hold on; box on;
for b = 1:3, plot(res.t/60, res.E(b,:)/1e3, 'LineWidth', 1.3, 'HandleVisibility', 'off'); end
xlabel('Time (min)'); ylabel('Stored energy (kJ m^{-1})');
for ax = [ax1 ax2 ax3]
    for tw = [30 60], xline(ax, tw, 'k--', 'HandleVisibility', 'off'); end
    xlim(ax, [0 res.t(end)/60]);
end
save_fig(fig, cfg, 'fig3_transient_example');
end
