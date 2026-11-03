function fig = plot_energy_vs_loading(cfg, pcm_name, kscen, h)
%PLOT_ENERGY_VS_LOADING  Phase-averaged stored energy vs EG loading.
%   Constant vs sine vs square (A = 20 K) at 30 and 60 min windows, with the
%   source period equal to the window (P = t_c), where four phases are converged.
L = load(fullfile(cfg.results_dir, sprintf('sweep_%s_%s_h%d_%s.mat', pcm_name, kscen, h, cfg.mode)));
fig = figure('Color', 'w', 'Position', [100 100 900 380]);
wins = [1800 3600];
names = {'constant', 'sine', 'square'};
for iw = 1:2
    keys = [0 0 3600; 1 20 wins(iw); 2 20 wins(iw)];
    subplot(1,2,iw); hold on; box on;
    E = energy_at_window(L.res, wins(iw));
    for g = 1:3
        cidx = find(ismember(L.meta.cases(:,1:3), keys(g,:), 'rows'));
        Ew = zeros(numel(cfg.w), numel(cidx));
        for q = 1:numel(cidx)
            sel = L.meta.case == cidx(q);
            [~, o] = sort(L.meta.w(sel)); e = E(sel); Ew(:,q) = e(o);
        end
        Ebar = mean(Ew, 2);
        [wo, Eo] = optimum_loading(cfg.w, Ebar);
        p = plot(100*cfg.w, Ebar/1e3, 'LineWidth', 1.4, 'DisplayName', names{g});
        plot(100*wo, Eo/1e3, 'o', 'Color', p.Color, 'MarkerFaceColor', p.Color, 'HandleVisibility', 'off');
    end
    xline(100*cfg.w_formstable, 'k:', 'form-stable', 'HandleVisibility', 'off');
    xlabel('EG mass fraction w (%)'); ylabel('Stored energy (kJ m^{-1})');
    title(sprintf('Charging window %d min', wins(iw)/60), 'FontWeight', 'normal');
    if iw == 1, legend('Location', 'southwest'); end
end
save_fig(fig, cfg, sprintf('fig4_E_vs_w_%s_%s', pcm_name, kscen));
end
