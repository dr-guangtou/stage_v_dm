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

## Maintenance rule
When new files are added under `outputs/intermediate/` or `outputs/tables/`, append a one-line description here immediately to preserve reproducibility context.
