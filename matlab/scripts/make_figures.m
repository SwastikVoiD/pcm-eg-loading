function make_figures(cfg)
%MAKE_FIGURES  Produce every figure used in the paper.
ALL = load_summary(cfg);
plot_k_models(cfg);
plot_validation(cfg);
plot_transient_example(cfg);
plot_energy_vs_loading(cfg, 'erythritol', 'low', cfg.h_list(1));
plot_energy_vs_loading(cfg, 'paraffin', 'mid', cfg.h_list(1));
base = ALL(~ALL.bypass & ~ALL.convection, :);
plot_optimum_map(cfg, base);
plot_penalty_vs_ratio(cfg, ALL(~ALL.convection, :));
plot_bypass(cfg, ALL(~ALL.convection, :));
plot_convection_bound(cfg, ALL(~ALL.bypass, :));
end
