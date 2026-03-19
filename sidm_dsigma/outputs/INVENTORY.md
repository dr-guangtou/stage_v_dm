# Outputs Inventory

Last updated: 2026-03-19 (MST)

This file records generated artifacts in `outputs/tables/` and `outputs/intermediate/`.

## Naming legend
- `dwarf` -> benchmark halo with `M200 = 1e10 Msun`.
- `cluster` -> benchmark halo with `M200 = 1e14 Msun`.
- `sigma_0p200` -> `sigma/m = 0.200 cm^2/g` (where `p` encodes decimal point).
- `profile` files -> 3D density profile comparison outputs (`rho(r)`).
- `delta_sigma` files -> projected lensing comparison outputs (`DeltaSigma(R)`).

## Tables (`outputs/tables`)

- `benchmark_summary.csv`: Main benchmark summary table for CDM vs SIDM runs (per benchmark and `sigma/m`) with `Delta chi^2`, SNR, required `3sigma` uniform precision, and max ratio shift.
- `nfw_reference_backend_status.csv`: Availability/execution status of optional reference backends (`colossus`, `pyccl`) used during NFW projection cross-checking.
- `nfw_reference_crosscheck.csv`: Radius-by-radius comparison of local numeric NFW `DeltaSigma` against analytic NFW reference (and optional external reference if available), including fractional differences.
- `precision_sweep.csv`: Required uniform per-bin fractional precision for target separations (`2sigma`, `3sigma`, `5sigma`) across radial windows (`full`, `inner`, `outer`) and SIDM grid points.
- `ensemble_delta_chi2_summary.csv`: Tier-1 halo-ensemble stacked forecast summary over `sigma/m = [0.2, 0.5, 1.0, 2.0]` with `Delta chi^2`, sigma separation (`5%` and `10%` toy precision scenarios), required uniform `3sigma` precision, max stacked ratio shift, and sampled-ensemble mass statistics.
- `cluster_ensemble_delta_chi2_summary.csv`: YAML-configured cluster (HMF mode) stacked forecast summary with `Delta chi^2`, sigma separation, required `3sigma` precision, and mass-statistics diagnostics.
- `cluster_ensemble_validation_checks.csv`: Validation checks for cluster YAML run (median mass target, radial-bin consistency, reproducibility with fixed seed, and `N_halos=1` single-halo limit).
- `dwarf_ensemble_delta_chi2_summary.csv`: YAML-configured dwarf (SHMR mode) stacked forecast summary with `Delta chi^2`, sigma separation, required `3sigma` precision, and mass-statistics diagnostics.
- `dwarf_ensemble_validation_checks.csv`: Validation checks for dwarf YAML run (median mass target, radial-bin consistency, reproducibility with fixed seed, and `N_halos=1` single-halo limit).
- `cluster_dwarf_ensemble_delta_chi2_summary.csv`: Combined table concatenating cluster + dwarf YAML forecast summaries for direct cross-regime comparison.

## Intermediate (`outputs/intermediate`)

- `cluster_delta_sigma_sigma_0p000.csv`: Cluster `DeltaSigma(R)` comparison table for CDM baseline (`sigma/m = 0.0`); includes CDM and model columns.
- `cluster_delta_sigma_sigma_0p200.csv`: Cluster `DeltaSigma(R)` comparison table for `sigma/m = 0.2`.
- `cluster_delta_sigma_sigma_0p500.csv`: Cluster `DeltaSigma(R)` comparison table for `sigma/m = 0.5`.
- `cluster_delta_sigma_sigma_1p000.csv`: Cluster `DeltaSigma(R)` comparison table for `sigma/m = 1.0`.
- `cluster_delta_sigma_sigma_2p000.csv`: Cluster `DeltaSigma(R)` comparison table for `sigma/m = 2.0`.

- `cluster_profile_sigma_0p000.csv`: Cluster 3D density profile comparison table (`rho(r)`) for CDM baseline (`sigma/m = 0.0`).
- `cluster_profile_sigma_0p200.csv`: Cluster 3D density profile comparison table (`rho(r)`) for `sigma/m = 0.2`.
- `cluster_profile_sigma_0p500.csv`: Cluster 3D density profile comparison table (`rho(r)`) for `sigma/m = 0.5`.
- `cluster_profile_sigma_1p000.csv`: Cluster 3D density profile comparison table (`rho(r)`) for `sigma/m = 1.0`.
- `cluster_profile_sigma_2p000.csv`: Cluster 3D density profile comparison table (`rho(r)`) for `sigma/m = 2.0`.

- `dwarf_delta_sigma_sigma_0p000.csv`: Dwarf `DeltaSigma(R)` comparison table for CDM baseline (`sigma/m = 0.0`); includes CDM and model columns.
- `dwarf_delta_sigma_sigma_0p200.csv`: Dwarf `DeltaSigma(R)` comparison table for `sigma/m = 0.2`.
- `dwarf_delta_sigma_sigma_0p500.csv`: Dwarf `DeltaSigma(R)` comparison table for `sigma/m = 0.5`.
- `dwarf_delta_sigma_sigma_1p000.csv`: Dwarf `DeltaSigma(R)` comparison table for `sigma/m = 1.0`.
- `dwarf_delta_sigma_sigma_2p000.csv`: Dwarf `DeltaSigma(R)` comparison table for `sigma/m = 2.0`.

- `dwarf_profile_sigma_0p000.csv`: Dwarf 3D density profile comparison table (`rho(r)`) for CDM baseline (`sigma/m = 0.0`).
- `dwarf_profile_sigma_0p200.csv`: Dwarf 3D density profile comparison table (`rho(r)`) for `sigma/m = 0.2`.
- `dwarf_profile_sigma_0p500.csv`: Dwarf 3D density profile comparison table (`rho(r)`) for `sigma/m = 0.5`.
- `dwarf_profile_sigma_1p000.csv`: Dwarf 3D density profile comparison table (`rho(r)`) for `sigma/m = 1.0`.
- `dwarf_profile_sigma_2p000.csv`: Dwarf 3D density profile comparison table (`rho(r)`) for `sigma/m = 2.0`.

- `ensemble_halo_catalog.csv`: Sampled Tier-1 halo catalog used for stacking, with one row per halo containing `m200_msun`, `z`, `c200`, and stacking `weight`.
- `ensemble_stacked_delta_sigma_profiles.csv`: Common-grid stacked `DeltaSigma(R)` table for CDM and SIDM ensemble forecasts; includes one SIDM column per benchmark `sigma/m`.
- `ensemble_stacked_rho_profiles.csv`: Common-grid stacked 3D density `rho(r)` table for CDM and SIDM ensemble forecasts; includes one SIDM column per benchmark `sigma/m`.
- `cluster_ensemble_halo_catalog.csv`: Halo catalog sampled from `docs/cluster_ensemble_config.yaml` (HMF mode), including `m200_msun`, `z`, `c200`, and equal-weight stack weights.
- `cluster_ensemble_stacked_delta_sigma_profiles.csv`: Cluster YAML stacked `DeltaSigma(R)` profiles (CDM + SIDM grid) on configured radial bins.
- `cluster_ensemble_stacked_rho_profiles.csv`: Cluster YAML stacked `rho(r)` profiles (CDM + SIDM grid) on common 3D radius grid.
- `dwarf_ensemble_halo_catalog.csv`: Halo catalog sampled from `docs/dwarf_ensemble_config.yaml` (SHMR mode), including inferred halo masses, redshift, concentration, and weights.
- `dwarf_ensemble_stacked_delta_sigma_profiles.csv`: Dwarf YAML stacked `DeltaSigma(R)` profiles (CDM + SIDM grid) on configured radial bins.
- `dwarf_ensemble_stacked_rho_profiles.csv`: Dwarf YAML stacked `rho(r)` profiles (CDM + SIDM grid) on common 3D radius grid.

## Maintenance rule
When new files are added under `outputs/intermediate/` or `outputs/tables/`, append a one-line description here immediately to preserve reproducibility context.
