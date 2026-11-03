# MATLAB code

Produces every result and figure in the paper. Requires MATLAB R2020a or newer (uses
`exportgraphics`, string arrays and implicit expansion); no toolboxes.

## Run

```matlab
cd matlab
main          % set MODE at the top: 'quick' (smoke test) or 'full' (paper results)
```

`main.m` runs eight steps in order:

1. `run_validation_neumann`: Neumann analytical solution at 30 and 120 min
2. `run_independence`: grid and time-step independence
3. `run_sweep`: parametric sweep (core, bypass and convection-bound batches)
4. `run_analysis`: builds `results/summary_all.csv`
5. `run_li2020_comparison`: calibrated consistency check against Li et al. (2020)
6. `run_phase_convergence`: 4, 8 and 16 start phases at P = 60 min, t_c = 30 and 60 min
7. `run_mushy_sensitivity`: mushy-band width 0.5, 1 and 2 K
8. `make_figures`: all figures into `figures/`

Finished sweep batches are cached in `results/*.mat`, so an interrupted run resumes and later steps
can be rerun on their own:

```matlab
setup_paths(); cfg = default_config('full');
make_figures(cfg);
```

## Files

| Folder | File | Purpose |
|---|---|---|
| `.` | `main.m`, `setup_paths.m` | Run everything; add folders to the path |
| `config` | `default_config.m` | Every study parameter (geometry, sweep ranges, sources, windows, batches) |
| `src/materials` | `pcm_properties.m` | Pure PCM data with sources (erythritol, n-docosane, Li 2020 salt) |
| | `composite_properties.m` | Effective ρ, c_p, L and k of PCM/EG composites |
| | `k_scenario_beta.m`, `fit_k_slope.m` | Conductivity slopes for the three routes, and their fit to published data |
| | `liquid_properties.m`, `convection_params.m` | Melt properties and inputs for the convection bound |
| | `k_maxwell.m`, `mass_to_volume_fraction.m`, `eg_properties.m` | Helpers |
| `src/numerics` | `solve_pcm_1d.m` | Batched implicit enthalpy solver (options: bypass, convection) |
| | `tdma_batch.m` | Thomas algorithm for many tridiagonal systems at once |
| | `convection_multiplier.m` | Raithby-Hollands effective-conductivity factor for the melt layer |
| | `make_grid.m`, `stored_energy.m`, `melt_fraction.m` | Grid and post-processing |
| `src/sources` | `source_temperature.m`, `build_cases.m`, `build_batch.m` | HTF profiles, case table, batch assembly |
| `src/analysis` | `optimum_loading.m` | Optimum from a pchip interpolant, optional form-stable bound |
| | `summarize_sweep.m`, `energy_at_window.m` | Phase-averaged optimum, shift, penalty, capacity change |
| | `load_summary.m` | Reads `summary_all.csv` with the logical columns restored |
| `src/validation` | `neumann_lambda.m`, `neumann_solution.m`, `melt_front_position.m` | Analytical benchmark |
| `src/plotting` | `plot_*.m`, `save_fig.m` | One function per paper figure |

## Outputs

`results/` holds the cached `.mat` batches and these CSVs: `summary_all.csv`, `independence.csv`,
`li2020_comparison.csv`, `phase_convergence_Pgt_tc.csv`, `mushy_sensitivity.csv`. The CSVs from the
reported run are in [`../results/matlab/`](../results/matlab/). `figures/` receives PNG (300 dpi)
and `.fig` files.

## Notes

- A full run is 20 batches of 2,499 simulations each; runtime depends on the machine.
- Property values marked `[VERIFY]` in `pcm_properties.m` and `liquid_properties.m` are typical
  literature values still to be confirmed against primary sources.
