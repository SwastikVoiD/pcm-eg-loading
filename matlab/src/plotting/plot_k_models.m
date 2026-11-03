function fig = plot_k_models(cfg)
%PLOT_K_MODELS  Conductivity gain (k - k_pcm) vs EG loading: the three
%   preparation-route slopes, the Maxwell lower bound, and published data.
[beta, data] = fit_k_slope();
eg = eg_properties(); pcm = pcm_properties('erythritol');
w = linspace(0.002, 0.25, 200);
fig = figure('Color', 'w', 'Position', [100 100 560 420]); hold on; box on;
cols = lines(3); names = {'low', 'mid', 'high'};
for i = 1:3
    b = beta.(names{i});
    plot(100*w, b*w, '-', 'Color', cols(i,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%s: \\beta = %.1f W/mK', names{i}, b));
end
phi = mass_to_volume_fraction(w, pcm.rho, eg.rho);
plot(100*w, k_maxwell(pcm.k, eg.k_particle, phi) - pcm.k, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Maxwell, dispersed spheres');
plot(100*data.sari.w(2:end), data.sari.k(2:end) - data.sari.k(1), 'o', 'Color', cols(1,:), 'MarkerFaceColor', cols(1,:), 'DisplayName', 'Sari & Karaipekli 2007 (paraffin)');
plot(100*data.li2024.w(2), diff(data.li2024.k), 's', 'Color', cols(2,:), 'MarkerFaceColor', cols(2,:), 'DisplayName', 'Li et al. 2024 (paraffin)');
plot(100*data.yuan2019.w(2), diff(data.yuan2019.k), '^', 'Color', cols(3,:), 'MarkerFaceColor', cols(3,:), 'DisplayName', 'Yuan et al. 2019 (erythritol)');
set(gca, 'YScale', 'log'); ylim([0.01 50]); xlim([0 25]);
xlabel('EG mass fraction w (%)'); ylabel('Conductivity gain k - k_{pcm} (W m^{-1} K^{-1})');
legend('Location', 'southeast', 'FontSize', 8);
save_fig(fig, cfg, 'fig1_k_models');
end
