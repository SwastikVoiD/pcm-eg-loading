# Python mirror

A line-for-line Python implementation of the MATLAB solver and analysis. It was used to develop and
cross-check the study; the reported results are the MATLAB ones in
[`../results/matlab/`](../results/matlab/), and the two agree to round-off.

```bash
pip install -r requirements.txt
```

## Cross-check (from the repository root)

```bash
python python/compare_matlab_python.py
```

Rewrites [`results/comparison/matlab_vs_python.md`](../results/comparison/matlab_vs_python.md) from
the CSVs in `results/`.

## Running the solver (from this folder)

Scripts write to `results/` and `figures/` inside this folder.

| Script | Purpose |
| --- | --- |
| `pcm1d.py` | Solver, materials, sources, convection model, Neumann solution |
| `analysis.py` | Optimum loading (pchip) and per-source summary |
| `val_neumann.py` | Neumann benchmark at three grid sizes |
| `indep.py` | Grid and time-step independence |
| `sweep.py` | `run(pcm, route, h, bypass=False, conv=False)` for one batch; `python sweep.py 1` runs the first batch |
| `sweep_conv.py` | The six convection-bound batches |
| `build_summary.py` | Builds `results/summary_all.csv` from the batch files |
| `li2020_cal.py` | Calibrated consistency check against Li et al. (2020) |
| `phase16.py` | 4/8/16 start phases for P = 60 min, t_c = 30 min |
| `review_checks.py` | Mushy-band and phase-count checks |
| `figures.py` | Figure functions (`fig1_k_models()` ... `fig8_convection_bound()`) |
| `compare_matlab_python.py` | MATLAB vs Python cross-check |

Each batch (2,499 simulations) takes about 3-7 minutes on one core. Example:

```python
from sweep import run
run('erythritol', 'low', 1000.0)          # writes results/sweep_erythritol_low_h1000.npz
```
