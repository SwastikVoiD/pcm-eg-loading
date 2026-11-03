function save_fig(fig, cfg, name)
%SAVE_FIG  Save a figure as PNG (300 dpi) and MATLAB .fig.
exportgraphics(fig, fullfile(cfg.figures_dir, [name '.png']), 'Resolution', 300);
savefig(fig, fullfile(cfg.figures_dir, [name '.fig']));
end
