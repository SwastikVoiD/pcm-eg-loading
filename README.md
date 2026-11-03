# Optimal expanded-graphite loading of composite PCM storage under intermittent waste-heat supply

Code, results and findings for the InDecarb 2026 paper (1st International Industry-Academia
Conference on Industrial Decarbonization, MIT Manipal with Forbes Marshall; Track 6, Thermal
Energy Utilization and Storage).

**Authors:** Swastik Debesh Patnaik, Navyaa Nayak, Rushi Mahajan, Arushi and N Govindha Rasu
(corresponding author, ngrasu@vit.ac.in), School of Mechanical Engineering, Vellore Institute of
Technology, Vellore, India.

## The question

Industrial waste heat is often intermittent, but composite phase change material (PCM) stores are
usually sized for a steady source. Expanded graphite (EG) raises the conductivity of the PCM but
displaces latent heat, so there is a filler loading that maximises the energy stored within a
charging window. This study asks whether an intermittent source moves that optimum.

## Main findings

Erythritol/EG and n-docosane/EG composites, three EG preparation routes, constant, sinusoidal and
square-wave heat-transfer-fluid (HTF) temperatures, 30-120 min charging windows. Full tables are in
[docs/FINDINGS.md](docs/FINDINGS.md).

1. **With continuous HTF flow and a charging window at least one source period long, fluctuation
   barely moves the optimum.** A loading chosen for the steady mean source loses at most 0.45 % of
   the achievable stored energy with conduction only, and 1.09 % with the natural-convection upper
   bound.
2. **Source periods longer than the window matter.** The design penalty rises to 6.3 % and the
   achievable energy can fall by up to 18 % (16 start phases).
3. **The charging window, preparation route and melt convection set the optimum.** For erythritol
   on the low-conductivity route with conduction only, the optimum falls from 18 wt% at 30 min to
   zero at 120 min. With convection, unconstrained optima lie at 0-2.5 wt%.
4. **For form-stable composites (at least 10 wt% EG) the optimum does not depend on convection.**
5. **An idealised bypass rule** (stop the HTF when it is colder than the store) lets fluctuating
   sources store up to 12.4 % more energy, but the steady-source design then loses up to 5.2 %.

## Repository layout

| Path | Contents |
| --- | --- |
| [`matlab/`](matlab/) | MATLAB project that produces every result and figure; run `main.m` |
| [`python/`](python/) | Line-for-line Python mirror of the solver and analysis, plus the cross-check script |
| [`results/matlab/`](results/matlab/) | MATLAB output CSVs (the reported results) |
| [`results/python_reference/`](results/python_reference/) | Python mirror output for the same cases |
| [`results/comparison/`](results/comparison/) | MATLAB vs Python cross-check report |
| [`figures/`](figures/) | Paper figures (Figs 1-8) |
| [`docs/`](docs/) | [Findings](docs/FINDINGS.md), [model and assumptions](docs/MODEL.md), [literature](docs/LITERATURE.md) |

## Quick start

**MATLAB** (R2020a or newer, no toolboxes):

```matlab
cd matlab
main          % set MODE at the top: 'quick' for a smoke test, 'full' for the paper results
```

Results go to `matlab/results/`, figures to `matlab/figures/`. Finished sweep batches are cached, so
an interrupted run resumes where it stopped. See [matlab/README.md](matlab/README.md).

**Python** (3.10 or newer):

```bash
pip install -r python/requirements.txt
python python/compare_matlab_python.py      # re-check MATLAB against the Python reference
```

See [python/README.md](python/README.md) for running the solver itself.

## Verification

- Neumann two-phase analytical solution: melt-front error 0.68 % at 30 min and 0.27 % at 120 min.
- Grid and time step: 80 cells at 2 s are within 0.008 % of 160 cells.
- Energy conserved to about 1e-11 (relative) in every run.
- Calibrated consistency check against the CFD of Li et al. (2020): melting times within 8 %.
- Mushy-band width 0.5-2 K: optimum changes by at most 0.5 wt%.
- MATLAB and Python agree to round-off on all 702 shared result rows
  ([report](results/comparison/matlab_vs_python.md)).

## Results version

The results and figures here are from the MATLAB run of 3 October 2026, with the material
properties listed in [docs/MODEL.md](docs/MODEL.md). The submitted extended abstract quotes slightly
different values from a later run (for example 1.04 % rather than 1.09 % for the convection-bound
penalty, and 13.4 % rather than 12.4 % for the bypass gain). Those values are not reproduced by the
code and data in this repository.

## Citation

See [CITATION.cff](CITATION.cff). Please cite the InDecarb 2026 paper once published.

## License

Code is released under the MIT License ([LICENSE](LICENSE)).
