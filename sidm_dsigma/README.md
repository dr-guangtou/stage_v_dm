# SIDM vs CDM Stage-V Weak-Lensing Forecast

`sidm_dsigma` is a lightweight, transparent forecasting pipeline for comparing CDM and SIDM lensing signals using stacked halo ensembles.

## Current Scope (Implemented)

The codebase supports three modeling tiers:

- Tier-1: SIDM-modified inner halo only.
- Tier-2: SIDM inner halo + DK14-like outer profile stitching.
- Tier-3: Tier-2 hybrid profile + empirical SIDM outskirts correction.

This is a controlled science-case forecast, not a full survey realism pipeline.
Tier-3 empirical outskirts corrections are disabled by default in the forecast runner and are intended for explicit sensitivity studies.

## Core Capabilities

- Ensemble generation with two modes:
  - `HMF` (halo-mass-selected, cluster-like)
  - `SHMR` (stellar-mass-selected, dwarf-like)
- SIDM profile generation with `parametricSIDM` backend support:
  - effective cross-section mode
  - velocity-dependent mode (`sigma0/m`, `w`)
- CDM baseline options:
  - NFW baseline
  - parametric SIDM tiny-cross-section reference mode
- Numerical projection to `Sigma(R)` and `DeltaSigma(R)` using local `numpy/scipy` integrals.
- Weighted stacking on a common radial grid.
- Toy distinguishability metrics (`Delta chi^2`, required fractional precision).
- Tier-3 diagnostics and redshift-overlay summary plotting.

## Main Entry Points

- Ensemble forecast runner:
  - `scripts/run_ensemble_forecast.py`
  - pass `--enable-tier3-empirical` only when intentionally running Tier-3 empirical correction scenarios
- Tier-3 redshift-overlay summary builder:
  - `scripts/build_tier3_redshift_overlay_summaries.py`
- Optional validation tools:
  - `scripts/run_reference_crosscheck.py`
  - `scripts/run_precision_sweep.py`
  - `scripts/run_cdm_engine_convergence_check.py`

## Configuration

YAML-driven runs are parsed by `src/sidm_stagev_forecast/ensemble_yaml.py`.

Key config blocks:

- `mode`: `HMF` or `SHMR`
- `sidm`:
  - `parameterization`: `effective` or `velocity_dependent`
  - `sigma_over_m_grid` or `sigma0_over_m_grid`
  - `w_km_s`, `time_model`, `mass_definition`
  - optional `cdm_reference` block
- `tier2`: DK14-like outer profile + stitch controls
- `tier3`: empirical outskirts correction controls

Split ensemble configs currently included:

- Cluster: `docs/cluster1_z0p1_ensemble_config.yaml`, `docs/cluster1_z0p5_ensemble_config.yaml`, `docs/cluster2_z0p1_ensemble_config.yaml`, `docs/cluster2_z0p5_ensemble_config.yaml`
- Dwarf: `docs/dwarf1_z0p05_ensemble_config.yaml`, `docs/dwarf1_z0p2_ensemble_config.yaml`, `docs/dwarf2_z0p05_ensemble_config.yaml`, `docs/dwarf2_z0p2_ensemble_config.yaml`

## Current Figure Set Policy

`outputs/figures/` is intentionally pruned in this phase to keep only default forecast summary figures:

- `cluster_ensemble_summary.png`
- `dwarf_ensemble_summary.png`

Captions for retained figures are documented in `outputs/figures/CAPTION.md`.

## Units and Conventions

- Core halo mass convention: `M200c`.
- Tier-2 matching may use `R200m` if configured.
- Core profile tables use physical `Msun`, `kpc` conventions.
- Tier-3 overlay summary lensing panel is displayed in `h Msun / pc^2`.

## Non-Goals (Current)

Not implemented as first-class physics in this stage:

- 2-halo term
- splashback calibration from first-principles SIDM cosmological simulations
- baryons/HOD/subhalo modeling
- full survey covariance realism
