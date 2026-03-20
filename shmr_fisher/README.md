# SHMR Forecast

Fisher matrix forecast tool for predicting how well spectroscopic surveys (combined with Stage-IV imaging lensing) can constrain the stellar-halo mass relation (SHMR).

## What it does

Given a spectroscopic survey configuration (area, redshift range, galaxy count, mass completeness), this tool computes:

- **Galaxy-galaxy lensing signal** DeltaSigma(R) via the halo model (HOD + HMF + NFW 1-halo + analytic 2-halo from colossus)
- **Clustering summary statistics** (effective halo bias, galaxy number density) per stellar mass bin
- **Fisher information matrix** for the Moster+2013 SHMR with optional mass-dependent scatter (Cao & Tinker 2020)
- **Marginalized parameter constraints** with optional systematic errors (shear calibration, photo-z bias, fractional floor)

The tool is designed for **relative survey comparisons** (Stage-III vs IV vs V), not absolute error bars.

## Quick start

```bash
# Install dependencies
uv sync

# Run forecast from YAML config (preferred)
uv run python scripts/run_forecast.py --config configs/default.yaml

# Generate all science figures
uv run python scripts/generate_science_figures.py --config configs/default.yaml

# 2-survey comparison
uv run python scripts/run_forecast.py --config configs/stage4_vs_stage5.yaml

# Parameter sweep
uv run python scripts/run_sweep.py --config configs/default.yaml --survey-key stage4_low_z --param area_deg2 --values 1000,5000,10000,14000

# Legacy mode (still works)
uv run python scripts/run_forecast.py --surveys stage4_low_z stage5_wide --systematics
```

## Package structure

```
shmr_fisher/
    __init__.py          # Sets Planck18 cosmology at import
    config.py            # SHMRParams, SurveyConfig, LensingConfig, ForecastConfig, NuisanceConfig
    config_io.py         # YAML config loader: RunConfig, load_run_config()
    shmr_model.py        # Moster+2013 SHMR: forward, inverse, scatter
    halo_model.py        # HOD, DeltaSigma (1h+2h), b_eff, n_gal via HMF integrals
    covariance.py        # Shape noise, Sigma_crit, survey volume, clustering noise
    fisher.py            # Fisher matrix computation, derivatives, marginalization
    systematics.py       # Systematic error helpers (floor, nuisance derivatives)
    survey_configs.py    # 5 predefined surveys + make_custom_survey() factory (legacy)
    validate.py          # Phase 1 validation checks
    plot_results.py      # All visualization functions (dynamic survey colors)
configs/
    default.yaml         # 4-survey comparison (Stage-III through Stage-V)
    stage4_vs_stage5.yaml  # 2-survey head-to-head
    dwarf_regime.yaml    # Single low-z survey for dwarf science
scripts/
    run_forecast.py      # Main driver: --config YAML or legacy CLI
    run_sweep.py         # Parameter sweep driver: --config or --base
    generate_science_figures.py  # Full figure generation from YAML config
outputs/
    {run_name}/          # Per-run directory (e.g., outputs/default/)
        config.yaml      # Copy of input config for reproducibility
        forecast_results.json, .npz
        *.png            # Figures (300 dpi)
        CAPTION.md       # Figure captions for this run
    phase1/...phase4/    # Legacy phase-based outputs (from development)
    CAPTION.md           # Master figure captions
```

## YAML configuration

Survey configurations are defined in YAML files under `configs/`. Each config can define 1-4 surveys. The config filename becomes the run name, and all outputs go to `outputs/{run_name}/`.

Example (`configs/stage4_vs_stage5.yaml`):
```yaml
shmr_params:
  use_mass_dependent_scatter: true

forecast:
  systematic_floor_fraction: 0.05
  include_nuisance_params: true
  vary_z_evolution: false

nuisance:
  sigma_m: 0.02
  sigma_dz_source: 0.03

surveys:
  stage4_low_z:
    name: "Stage-IV Low-z"
    area_deg2: 14000
    z_min: 0.05
    z_max: 0.4
    n_gal_total: 10000000
    log_Mstar_min: 9.0

  stage5_wide:
    name: "Stage-V Wide"
    area_deg2: 10000
    z_min: 0.05
    z_max: 1.0
    n_gal_total: 50000000
    log_Mstar_min: 8.5
    log_Mstar_min_z_dep:
      base: 8.5
      slope: 1.0
```

## SHMR parameterization

Uses the Moster+2013 double power-law with z/(1+z) evolution:

    M*/Mh = 2N(z) [(Mh/M1(z))^(-beta(z)) + (Mh/M1(z))^(gamma(z))]^(-1)

Core parameters: log_M1_0, N_0, beta_0, gamma_0 (z=0 shape), nu_M1, nu_N, nu_beta, nu_gamma (evolution).

Scatter: either constant (sigma_logMs = 0.15 dex) or mass-dependent following Cao & Tinker (2020), parameterized as sigma_high + sigma_rise * [1 - tanh((log_Mh - log_Mh_break) / delta)].

## Systematic error options

Both are toggleable via YAML config or `ForecastConfig`:

- **Fractional floor** (`systematic_floor_fraction`): adds (f_sys * DS_fid)^2 to the lensing variance per radial bin. Captures baryonic effects, miscentering, etc.
- **Nuisance marginalization** (`include_nuisance_params`): adds shear calibration bias (m) and source photo-z bias (dz) to the Fisher matrix with Gaussian priors from `NuisanceConfig`.

## Predefined surveys

| Name | Area | z-range | N_gal | log M*_min |
|------|------|---------|-------|-----------|
| Stage-III Shallow Wide | 7500 deg^2 | 0.02-0.2 | 700k | 9.5 |
| Stage-IV Low-z | 14000 deg^2 | 0.05-0.4 | 10M | 9.0 |
| Stage-IV High-z | 14000 deg^2 | 0.4-1.0 | 5M | 10.8 |
| Stage-V Wide | 10000 deg^2 | 0.05-1.0 | 50M | 8.5 (z-dep) |
| Stage-V Deep | 3000 deg^2 | 0.1-1.5 | 20M | 8.0 (z-dep) |

## Key dependencies

- **colossus** (Diemer 2018): cosmology, HMF, halo bias, NFW profiles, c-M relation
- **numpy**, **scipy**: numerics, root-finding, integration
- **matplotlib**: visualization
- **astropy**: physical constants for Sigma_crit
- **pyyaml**: YAML configuration parsing

## References

- Moster, Naab & White (2013), MNRAS, 428, 3121
- Cao & Tinker (2020), MNRAS, 498, 5080
- Leauthaud et al. (2012), ApJ, 744, 159
- Tinker et al. (2008, 2010), ApJ
- Diemer (2018), ApJS, 239, 35
- Oguri & Takada (2011), PRD, 83, 023008
