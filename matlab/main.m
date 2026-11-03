%MAIN  Run the whole study: validation, independence, sweep, analysis, figures.
%
%   Composite PCM filler loading under intermittent heat supply (InDecarb 2026).
%   Set MODE = 'quick' for a fast smoke test, 'full' for the paper results.
%   Full mode is ~18 solver batches of ~2,500 cases each; runtime depends on the
%   machine (expect tens of minutes). Results are cached in results/*.mat, so
%   re-running skips finished batches.

clear; clc;
setup_paths();
MODE = 'full';                       % 'quick' or 'full'
cfg = default_config(MODE);

fprintf('=== 1. Verification against the Neumann solution ===\n');
run_validation_neumann(cfg);

fprintf('\n=== 2. Grid and time-step independence ===\n');
run_independence(cfg);

fprintf('\n=== 3. Parametric sweep ===\n');
run_sweep(cfg);

fprintf('\n=== 4. Optimum-loading analysis ===\n');
run_analysis(cfg);

fprintf('\n=== 5. Consistency check against Li et al. (2020) ===\n');
run_li2020_comparison(cfg);

fprintf('\n=== 6. Start-phase convergence (P > t_c) ===\n');
run_phase_convergence(cfg);

fprintf('\n=== 7. Mushy-band sensitivity ===\n');
run_mushy_sensitivity(cfg);

fprintf('\n=== 8. Figures ===\n');
make_figures(cfg);
fprintf('\nDone. Tables in results/, figures in figures/.\n');
