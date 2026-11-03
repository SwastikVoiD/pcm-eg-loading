function fig = plot_validation(cfg)
%PLOT_VALIDATION  Numerical vs Neumann temperature profiles at 30 and 120 min.
L = load(fullfile(cfg.results_dir, 'validation_neumann.mat')); V = L.V;
fig = figure('Color', 'w', 'Position', [100 100 560 380]); hold on; box on;
cols = lines(numel(V));
for i = 1:numel(V)
    near = V(i).x < 0.06;
    plot(1e3*V(i).x(near), V(i).T_exact(near), '-', 'Color', cols(i,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Neumann, t = %d min', V(i).t_min));
    idx = find(near); idx = idx(1:6:end);
    plot(1e3*V(i).x(idx), V(i).T_num(idx), 'o', 'Color', cols(i,:), 'MarkerSize', 4, ...
        'DisplayName', sprintf('Model, t = %d min', V(i).t_min));
end
yline(pcm_properties('erythritol').Tm, 'k:', 'HandleVisibility', 'off');
xlabel('Distance from heated wall (mm)'); ylabel('Temperature (^\circC)');
legend('Location', 'northeast', 'FontSize', 8);
save_fig(fig, cfg, 'fig2_validation');
end
