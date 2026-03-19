# SHMR Forecast

Fisher matrix forecast tool for predicting how well spectroscopic surveys (combined with Stage-IV imaging lensing) can constrain the stellar-halo mass relation (SHMR).

## What it does

Given a spectroscopic survey configuration (area, redshift range, galaxy count, mass completeness), this tool computes:

- **Galaxy-galaxy lensing signal** DeltaSigma(R) via the halo model (HOD + HMF + NFW profiles from colossus)
- **Clustering summary statistics** (effective halo bias, galaxy number density) per stellar mass bin
- **Fisher information matrix** for the 9-parameter Moster+2013 SHMR (4 z=0 shape + 4 evolution + 1 scatter)
- **Marginalized parameter constraints** with optional systematic errors (shear calibration, photo-z bias, fractional floor)

The tool is designed for **relative survey comparisons** (Stage-III vs IV vs V), not absolute error bars.

## Quick start

```bash
# Install dependencies
uv add colossus numpy scipy matplotlib astropy

# Run all predefined surveys (stat-only)
uv run python scripts/run_forecast.py

# Run with realistic systematic errors
uv run python scripts/run_forecast.py --systematics

# Run a specific subset
uv run python scripts/run_forecast.py --surveys stage4_low_z stage5_wide

# Parameter sweep
uv run python scripts/run_sweep.py --base stage5_wide --param n_gal_total --values 1e6,5e6,1e7,5e7,1e8
```

## Package structure

```
shmr_fisher/
    __init__.py          # Sets Planck18 cosmology at import
    config.py            # SHMRParams, SurveyConfig, LensingConfig, ForecastConfig, NuisanceConfig
    shmr_model.py        # Moster+2013 SHMR: forward, inverse, scatter
    halo_model.py        # HOD, DeltaSigma, b_eff, n_gal via HMF integrals
    covariance.py        # Shape noise, Sigma_crit, survey volume, clustering noise
    fisher.py            # Fisher matrix computation, derivatives, marginalization
    survey_configs.py    # 5 predefined surveys + make_custom_survey() factory
    validate.py          # Phase 1 validation checks
    plot_results.py      # All visualization functions
scripts/
    run_forecast.py      # Main driver (CLI)
    run_sweep.py         # Parameter sweep driver (CLI)
outputs/
    CAPTION.md           # Detailed captions for all figures
    *.png                # Generated figures (300 dpi PNG)
    *.json, *.npz        # Saved Fisher matrices and constraints
```

## SHMR parameterization

Uses the Moster+2013 double power-law with z/(1+z) evolution:

    M*/Mh = 2N(z) [(Mh/M1(z))^(-beta(z)) + (Mh/M1(z))^(gamma(z))]^(-1)

9 free parameters: log_M1_0, N_0, beta_0, gamma_0 (z=0 shape), nu_M1, nu_N, nu_beta, nu_gamma (evolution), sigma_logMs (scatter).

## Systematic error options

Both are toggleable via `ForecastConfig`:

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

## References

- Moster, Naab & White (2013), MNRAS, 428, 3121
- Leauthaud et al. (2012), ApJ, 744, 159
- Tinker et al. (2008, 2010), ApJ
- Diemer (2018), ApJS, 239, 35
- Oguri & Takada (2011), PRD, 83, 023008
