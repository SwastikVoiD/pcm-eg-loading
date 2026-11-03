# Findings

All values come from [`results/matlab/`](../results/matlab/) unless stated, at the base HTF-side
coefficient h = 1000 W m⁻² K⁻¹, with the HTF always flowing and conduction only, unless a column
or row says otherwise.

**Definitions.** w* is the EG mass fraction that maximises energy stored within the charging
window t_c. The *design penalty* is the energy lost by using the steady-source optimum under a
fluctuating source: 1 − E_fluct(w*_const) / E_fluct(w*_fluct). The *change in achievable energy* is
E_fluct(w*_fluct) / E_const(w*_const) − 1. P is the source period. Fluctuating cases are averaged
over four start phases unless stated.

## 1. Fluctuation and the optimum

| Regime | Largest optimum shift | Largest design penalty | Largest loss of achievable energy |
| --- | --- | --- | --- |
| P ≤ t_c, conduction only | 1.5 wt% | 0.45 % | 2.3 % |
| P ≤ t_c, convection upper bound | 2.0 wt% | 1.09 % | 0.8 % |
| P > t_c (P = 60 min, t_c = 30 min), 16 phases | 10.0 wt% | 6.34 % | 17.9 % |

When the window spans at least one source period, a store sized for the steady mean source is
within about 1 % of one optimised for the actual fluctuating source. Only when a single source cycle
is longer than the window does timing matter. Fig. 6 shows every case against P/t_c.

## 2. What sets the optimum

Constant source, w* by PCM, conductivity route and charging window:

| PCM | Route | w* 30 min | w* 60 min | w* 120 min | w* 30 min, convection | w* 60 min, convection | Form-stable w* 30 min |
|---|---|---|---|---|---|---|---|
| Erythritol | Low | 18.0 % | 4.0 % | 0.0 % | 2.5 % | 0.0 % | 18.0 % |
| Erythritol | Mid | 6.0 % | 1.5 % | 0.0 % | 1.0 % | 0.0 % | 10.0 % |
| Erythritol | High | 2.0 % | 0.5 % | 0.0 % | 0.5 % | 0.0 % | 10.0 % |
| Paraffin | Low | 11.0 % | 4.0 % | 0.5 % | 0.5 % | 0.0 % | 11.0 % |
| Paraffin | Mid | 3.0 % | 1.0 % | 0.5 % | 0.5 % | 0.0 % | 10.0 % |
| Paraffin | High | 1.0 % | 0.5 % | 0.5 % | 0.5 % | 0.0 % | 10.0 % |

- Shorter windows need faster charging, so more EG.
- Better-connected networks (higher β) reach the needed conductivity with less graphite.
- With the convection upper bound, the melt does the work EG would otherwise do, and the
  unconstrained optimum falls to 0-2.5 wt%.
- Most unconstrained optima lie below the 10 wt% form-stability limit. With that limit enforced,
  the optimum is the same under both convection bounds (last column).

## 3. Convection bounds by route

Largest design penalty for P ≤ t_c:

| PCM | Route | Conduction only | Convection upper bound |
|---|---|---|---|
| Erythritol | Low | 0.10 % | 0.52 % |
| Erythritol | Mid | 0.24 % | 0.09 % |
| Erythritol | High | 0.14 % | 0.00 % |
| Paraffin | Low | 0.23 % | 0.09 % |
| Paraffin | Mid | 0.33 % | 0.48 % |
| Paraffin | High | 0.45 % | 1.09 % |

## 4. Bypass control

Idealised rule: HTF flow stops during any step in which the HTF is colder than the store wall.
Erythritol (all routes) and paraffin (mid route), P ≤ t_c:

| Operating mode | Change in achievable energy | Largest design penalty |
| --- | --- | --- |
| HTF always flowing | −2.3 % to 0.0 % | 0.45 % |
| Bypass when cold | −0.4 % to +12.4 % | 5.22 % |

Bypass turns the capacity loss into a gain because the store is charged only by the hot
half-cycles. The steady mean source then becomes a poorer design basis (optimum shift up to
5.0 wt%). The rule assumes instantaneous switching with no delay or minimum flow, so real gains
would be smaller.

## 5. Sensitivity to the HTF-side coefficient (low route)

| PCM | h (W m⁻² K⁻¹) | w* 30 min | w* 60 min | Largest design penalty, P ≤ t_c |
|---|---|---|---|---|
| Erythritol | 200 | 25.0 % (edge of range) | 17.0 % | 0.08 % |
| Erythritol | 1000 | 18.0 % | 4.0 % | 0.10 % |
| Erythritol | 5000 | 13.0 % | 2.0 % | 0.32 % |
| Paraffin | 200 | 22.5 % | 6.5 % | 0.38 % |
| Paraffin | 1000 | 11.0 % | 4.0 % | 0.23 % |
| Paraffin | 5000 | 9.5 % | 3.5 % | 0.23 % |

A weaker HTF-side coefficient raises the optimum: more EG pulls heat away from the wall, keeping the
wall cooler and the driving temperature difference larger. The main finding holds at every h.

## 6. Numerical checks

**Neumann analytical solution** (pure erythritol, wall at T_m + 20 K, initial T_m − 10 K, λ = 0.2209):

| Time | Front, model | Front, exact | Error | Max temperature error |
| --- | --- | --- | --- | --- |
| 30 min | 9.44 mm | 9.37 mm | 0.68 % | 0.21 K |
| 120 min | 18.80 mm | 18.75 mm | 0.27 % | 0.24 K |

These values come from the Python mirror; `matlab/scripts/run_validation_neumann.m` runs the same
check and prints its result.

**Grid and time step** ([`independence.csv`](../results/matlab/independence.csv)), stored energy
after 60 min against 160 cells at 2 s: 40 cells 0.033 %, 80 cells 0.008 %, 80 cells at 1 s
0.007 %, 80 cells at 5 s 0.022 %. Production runs use 80 cells and 2 s.

**Calibrated consistency check against Li et al. (2020)**
([`li2020_comparison.csv`](../results/matlab/li2020_comparison.csv)). The HTF coefficient is
calibrated once (h = 37.1 W m⁻² K⁻¹) to their constant-source melting time; the fluctuating cases
are then predicted without further tuning:

| Sinusoidal source (mean 300 °C) | Li et al. (min) | Model (min) | Difference |
| --- | --- | --- | --- |
| A = 100 K, P = 1 min | 49.8 | 49.6 | −0.4 % |
| A = 100 K, P = 6 min | 49.3 | 49.2 | −0.3 % |
| A = 100 K, P = 60 min | 37.3 | 36.8 | −1.3 % |
| A = 50 K, P = 60 min | 44.2 | 45.5 | +2.9 % |
| A = 150 K, P = 60 min | 32.3 | 29.9 | −7.4 % |

Because h was fitted to one of their cases, this is a consistency check, not independent
validation. The stored-energy trend does not match their reported value (−2.0 % against −9.5 % at
P = 60 min).

**Start phases** ([`phase_convergence_Pgt_tc.csv`](../results/matlab/phase_convergence_Pgt_tc.csv)),
P = 60 min, all PCMs, routes and fluctuating sources:

| Charging window | Start phases | Max penalty | Max optimum shift | Largest capacity loss |
|---|---|---|---|---|
| 30 min (P > t_c) | 4 | 6.47 % | 11.5 wt% | −21.0 % |
| 30 min (P > t_c) | 8 | 6.41 % | 8.5 wt% | −18.5 % |
| 30 min (P > t_c) | 16 | 6.34 % | 10.0 wt% | −17.9 % |
| 60 min (P = t_c) | 4 | 0.29 % | 1.0 wt% | −1.9 % |
| 60 min (P = t_c) | 8 | 0.24 % | 1.0 wt% | −1.6 % |
| 60 min (P = t_c) | 16 | 0.22 % | 1.0 wt% | −1.5 % |

Four phases are converged when the window spans a full period. When P > t_c only the worst-case
envelope converges; individual cases still move by up to 5.5 wt% between 8 and 16 phases.

**Mushy-band width** ([`mushy_sensitivity.csv`](../results/matlab/mushy_sensitivity.csv)), 0.5 and
2 K against the 1 K default: optimum changes by at most 0.5 wt%, design penalty by at most 0.14
percentage points.

**MATLAB vs Python**: all 702 shared rows agree to round-off
([report](../results/comparison/matlab_vs_python.md)).

## 7. Limitations

- One-dimensional cross-section with no axial HTF temperature drop; results apply to the inlet
  section.
- Melt convection is bracketed by a quasi-steady correlation, not resolved. Real units lie between
  the bounds.
- Each conductivity route's β is transferred between the two PCMs; the routes are sensitivity
  bounds, not material-specific predictions.
- Latent heat is taken as (1 − w) L_pcm. Measured erythritol/EG composites can fall a few percent
  below this.
- The form-stability limit depends on the preparation route (about 10 to 20 wt% in the
  literature); 10 wt% is used here.
- Source profiles are synthetic sinusoids and square waves, not measured plant data.
- Only charging is modelled; erythritol's supercooling on discharge is not.
- The bypass controller is idealised.
- Several property values are typical literature values still to be checked against primary
  sources (see [MODEL.md](MODEL.md)).
