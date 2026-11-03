# Model, parameters and assumptions

## Geometry

One radial cross-section of a shell-and-tube unit, with the geometry of Li et al. (2020): an HTF
tube of radius r_i = 12.5 mm inside a PCM annulus of outer radius r_o = 25 mm. The outer wall is
adiabatic. All energies are per metre of tube.

## Governing equations

Enthalpy form of the energy equation:

```
ρ c_p ∂T/∂t + ρ L ∂f/∂t = (1/r) ∂/∂r ( r k ∂T/∂r ),        r_i < r < r_o
−k ∂T/∂r = h (T_HTF(t) − T)   at r = r_i
∂T/∂r = 0                     at r = r_o
```

The liquid fraction f rises linearly from 0 to 1 over a 1 K mushy band centred on the melting
point. The PCM starts 10 K below its melting point.

## Numerics

- 80 uniform finite volumes, implicit Euler, Δt = 2 s.
- Latent term linearised with the source-based method of Voller & Swaminathan (1991); each step
  iterates tridiagonal solves (Thomas algorithm) until f changes by less than 1e-9.
- Every filler loading and source case is solved together as one batch (2,499 cases per batch).
- The optimum is the maximum of a shape-preserving (pchip) interpolant through the 51 loading
  values at 1e-4 resolution; the same interpolant gives the design penalty.

## Composite properties

At EG mass fraction w:

```
L   = (1 − w) L_pcm
c_p = (1 − w) c_p,pcm + w c_p,EG
1/ρ = (1 − w)/ρ_pcm + w/ρ_EG
k   = k_pcm + β w
```

## Conductivity routes

Maxwell's dispersed-sphere model underpredicts every EG measurement found (Fig. 1), because EG forms
a connected network. The linear slope β is set by how the composite is made:

| Route | β (W m⁻¹ K⁻¹) | Anchor data |
| --- | --- | --- |
| Low: loose mixing and absorption | 6.38 | Paraffin/EG, 0-10 wt%, least-squares fit (Sari & Karaipekli 2007) |
| Mid: vacuum impregnation | 28.28 | Paraffin/EG, 6.0 W m⁻¹ K⁻¹ at 20 wt% (Li Y. et al. 2024) |
| High: impregnation, compression, sintering | 117.9 | Erythritol/EG, 12.5 W m⁻¹ K⁻¹ at 10 wt% (Yuan et al. 2019) |

`matlab/src/materials/fit_k_slope.m` reproduces these slopes from the published data. Each route is
applied to both PCMs as a modelling assumption, so the routes are sensitivity bounds spanning an
18-fold range of β. Independent data fall inside the bracket: Hirschey et al. (2021) β ≈ 8-12,
Yang et al. (2026) erythritol/EG β ≈ 15, Yu & Tao (2022) β ≈ 25, Fang et al. (2024) β ≈ 69 and
Mallow et al. (2017) β ≈ 70.

## Material properties

Constant, phase-averaged properties, as in Li et al. (2020). Values marked † are typical literature
values still to be confirmed against primary sources.

| Property | Erythritol | n-Docosane | EG |
| --- | --- | --- | --- |
| Melting point (°C) | 117 | 44 | - |
| Latent heat (kJ kg⁻¹) | 344 | 256 | - |
| Specific heat (J kg⁻¹ K⁻¹) | 2070 † | 1950 † | 710 |
| Density (kg m⁻³) | 1390 † | 800 † | 2200 |
| Conductivity (W m⁻¹ K⁻¹) | 0.72 | 0.22 | - |

Melt properties used only by the convection bound:

| Property | Erythritol | n-Docosane |
| --- | --- | --- |
| Conductivity (W m⁻¹ K⁻¹) | 0.326 | 0.15 |
| Density (kg m⁻³) | 1300 | 778 |
| Specific heat (J kg⁻¹ K⁻¹) | 2760 | 2380 |
| Viscosity (Pa s) | 0.016 † | 0.0036 † |
| Thermal expansion (K⁻¹) | 5.0e-4 † | 9.0e-4 † |

Yang et al. (2026), citing Wang et al. (2016), give 3.95e-4 K⁻¹ for erythritol's thermal expansion.
Using it would lower the convection enhancement by about 6 %, since the correlation scales with the
Rayleigh number to the power 1/4.

## Heat source

```
T_HTF(t) = T_m + 20 K + A · s(2πt/P + φ),     s = sin or sign(sin)
```

- Amplitude A = 10 and 20 K; period P = 10, 30 and 60 min; plus a constant source.
- The 20 K mean superheat matches the erythritol shell-and-tube tests of Agyenim et al. (2011)
  (HTF at 140 °C) and keeps erythritol below its 160 °C decomposition point.
- Four start phases φ = 0, π/2, π, 3π/2, results averaged; 16 phases for the P > t_c regime.
- h = 1000 W m⁻² K⁻¹ base; 200 and 5000 W m⁻² K⁻¹ as sensitivities.

## Design objective

```
E(w, t_c) = ∫ ρ [ c_p (T − T_0) + f L ] dV,      t_c = 30, 60, 120 min
```

Without a time limit pure PCM always stores the most, so the objective is the energy stored within
a fixed charging window. EG loading is swept from 0 to 25 wt% in 0.5 wt% steps (51 values).

## Natural-convection bound

The conduction-only model is the lower bound. For the upper bound, conductivity in the molten layer
between the tube wall and the melt front is multiplied by the Raithby-Hollands correlation for a
horizontal concentric annulus, evaluated each step from the melt-layer thickness L and wall
superheat:

```
k_eff / k = 0.386 (Pr / (0.861 + Pr))^(1/4) Ra_c^(1/4),     k_eff / k ≥ 1
Ra_c = [ln(D_o/D_i)]^4 / ( L^3 (D_i^(−3/5) + D_o^(−3/5))^5 ) · Ra_L
```

Applied to a full annular melt layer this overestimates real convection, since in a horizontal tube
the solid sinks and the lower region melts mainly by conduction. EG is assumed to suppress
convection, with the enhancement falling linearly to zero at the 10 wt% form-stable loading. Mallow
et al. (2017) found convection negligible at about 6 wt%, so this assumption is conservative.

## Bypass control

In the bypass runs the HTF-side coefficient is set to zero during any step in which T_HTF is below
the wall-cell temperature. Switching is instantaneous, with no delay or minimum flow.

## Simulation count

Each batch is 51 loadings × 49 source cases (1 constant + 2 shapes × 2 amplitudes × 3 periods ×
4 phases) = 2,499 transient simulations. The full MATLAB configuration runs 20 batches (6 core,
4 h-sensitivity, 4 bypass, 6 convection), about 50,000 simulations, each read at the three charging
windows.
