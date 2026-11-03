function fig = plot_convection_bound(cfg, ALL) %#ok<INUSD>
%PLOT_CONVECTION_BOUND  Stored energy vs loading with and without the
%   natural-convection model (erythritol, low route, base h): constant source and
%   square wave (A = 20 K, P = 30 min), 30 and 60 min windows.
f0 = fullfile(cfg.results_dir, sprintf('sweep_erythritol_low_h%d_%s.mat', cfg.h_list(1), cfg.mode));
f1 = fullfile(cfg.results_dir, sprintf('sweep_erythritol_low_h%d_conv_%s.mat', cfg.h_list(1), cfg.mode));
if ~exist(f1, 'file'), warning('plot_convection_bound: %s not found, skipped.', f1); fig = []; return; end
D = {load(f0), load(f1)}; styles = {'-', '--'}; labs = {'conduction only', 'with convection'};
keys = [0 0 3600; 2 20 1800]; knames = {'constant', 'square, A=20 K, P=30 min'};
cols = lines(2); wins = [1800 3600];
fig = figure('Color', 'w', 'Position', [100 100 900 380]);
for iw = 1:2
    subplot(1, 2, iw); hold on; box on;
    for d = 1:2
        L = D{d}; E = energy_at_window(L.res, wins(iw));
        for g = 1:2
            cidx = find(ismember(L.meta.cases(:,1:3), keys(g,:), 'rows'));
            Ew = zeros(numel(cfg.w), numel(cidx));
            for q = 1:numel(cidx)
                sel = L.meta.case == cidx(q);
                [~, o] = sort(L.meta.w(sel)); e = E(sel); Ew(:, q) = e(o);
            end
            Ebar = mean(Ew, 2);
            [wo, Eo] = optimum_loading(cfg.w, Ebar);
            plot(100*cfg.w, Ebar/1e3, styles{d}, 'Color', cols(g,:), 'LineWidth', 1.3, ...
                'DisplayName', sprintf('%s, %s', knames{g}, labs{d}));
            plot(100*wo, Eo/1e3, 'o', 'Color', cols(g,:), 'HandleVisibility', 'off');
        end
    end
    xline(100*cfg.w_formstable, 'k:', 'HandleVisibility', 'off');
    xlabel('EG mass fraction w (%)'); ylabel('Stored energy (kJ m^{-1})');
    title(sprintf('Charging window %d min', wins(iw)/60), 'FontWeight', 'normal');
    if iw == 1, legend('Location', 'southeast', 'FontSize', 7); end
end
save_fig(fig, cfg, 'fig8_convection_bound');
end
