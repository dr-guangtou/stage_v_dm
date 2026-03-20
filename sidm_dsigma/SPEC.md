# SPEC.md

## Project Scope

`sidm_dsigma` provides a modular, first-order SIDM lensing forecast pipeline centered on stacked halo ensembles.

Implemented modeling tiers:

- Tier-1: SIDM inner-halo modifications only.
- Tier-2: SIDM inner halo stitched to DK14-like outskirts.
- Tier-3: Tier-2 profile plus empirical SIDM-dependent outer correction.

The pipeline is designed for interpretability and fast iteration, not full survey realism.

## Scientific Framing

### Tier Interpretation

- Tier-1: inner halo response to SIDM.
- Tier-2: inner SIDM + realistic outer steepening baseline.
- Tier-3: empirical, simulation-informed nuisance extension for SIDM outskirts sensitivity studies.

Tier-3 is explicitly an approximate correction framework, not a first-principles cosmological SIDM splashback solver.

### Regimes and Selection Modes

- Cluster-like runs: `mode=HMF` (halo-mass-selected).
- Dwarf-like runs: `mode=SHMR` (stellar-mass-selected, central-only assumption in Tier-1 design).

## Architecture

Core package: `src/sidm_stagev_forecast/`

- `config.py`: benchmark/config dataclasses and defaults.
- `cosmology.py`: cosmology helpers and `M200c`/`Rdelta` support.
- `ensemble.py`: HMF and SHMR samplers, redshift/selection/scatter logic.
- `ensemble_yaml.py`: YAML normalization into runtime config containers.
- `profiles.py`: NFW baseline, parametric SIDM wrapper integration, Tier-2/Tier-3 builders.
- `outer_profiles.py`: DK14-like outer profile components.
- `stitch.py`: smooth log-density stitching utilities.
- `outer_corrections.py`: empirical Tier-3 outer corrections.
- `calibration.py`: correction preset/calibration helpers.
- `projection.py`: local numerical projection (`Sigma`, `DeltaSigma`).
- `stacking.py`: interpolation and weighted stacking.
- `forecast.py`: toy distinguishability metrics.
- `plotting.py`: figure generation, including Tier-3 redshift-overlay summary with inset distributions.
- `io.py`: output writing and caption/inventory appenders.
- `velocity_dependence.py`: SIDM velocity-dependent helper functions.

Primary scripts:

- `scripts/run_ensemble_forecast.py`
- `scripts/build_tier3_redshift_overlay_summaries.py`
- `scripts/run_cdm_engine_convergence_check.py`
- `scripts/run_reference_crosscheck.py`
- `scripts/run_precision_sweep.py`

## Configuration Contract

YAML is the operational interface for ensemble runs.

Required high-level blocks:

- `mode`: `HMF` or `SHMR`
- `ensemble`, `redshift`, `projection`, `stacking`, `concentration`
- `sidm`
- optional `tier2`
- optional `tier3`

### SIDM Block

Supported parameterizations:

- `effective`: uses `sigma_over_m_grid`
- `velocity_dependent`: uses `sigma0_over_m_grid` with `w_km_s`

Supported fields include:

- `parameterization`
- `w_km_s`
- `time_model`
- `mass_definition`
- optional `cdm_reference`:
  - `profile_source`
  - `sigma0_over_m`
  - `w_km_s`
  - `time_model`

### Tier-2 Block

Controls DK14-like outskirts + stitch model:

- `enabled`
- `outer_profile_model`
- `stitch_method`
- `r_match_mode`
- `r_match_value`
- `smooth_width_dex`
- `continuity`
- optional regime overrides

### Tier-3 Block

Controls empirical outer correction layer:

- `enabled`
- `correction_model`
- `sigma_pivot`
- `apply_to_regimes`
- `calibration_mode`
- `preset`
- model-specific parameter sub-blocks (`rt_shift`, `gamma_shift`, `beta_shift`, `outer_window`)

Runtime policy:

- Tier-3 empirical correction is opt-in at execution time.
- Default forecast runs do not apply Tier-3 empirical correction unless explicitly enabled via CLI.

## Numerical and Physical Conventions

- Core mass convention: `M200c`.
- `R200m` is used only when requested by Tier-2 matching mode.
- Projection is local numerical integration with `numpy/scipy`.
- External packages (`colossus`, `pyccl`) are optional references/validation, not mandatory runtime requirements.

## Output Contract

Standard outputs are stored in:

- `outputs/intermediate/`
- `outputs/tables/`
- `outputs/figures/`

Operational documentation requirements:

- Every generated figure must have an entry in `outputs/figures/CAPTION.md`.
- Generated intermediate/table artifacts are tracked in `outputs/INVENTORY.md`.

## Current Figure Retention Policy (Phase-1 Organization)

To reduce clutter and preserve canonical outputs, only default ensemble summary figures are retained:

- `cluster_ensemble_summary.png`
- `dwarf_ensemble_summary.png`

## Explicit Non-Goals

Still out of scope:

- 2-halo term
- baryonic/HOD/subhalo realism
- full covariance realism
- first-principles cosmological SIDM outskirts/splashback simulation in this repository
