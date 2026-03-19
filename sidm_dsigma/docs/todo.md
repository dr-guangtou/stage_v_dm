# todo

## session plan
- [x] Read README.md, SPEC.md, AGENTS.md, PLAN.md.
- [x] Bootstrap package structure and project tooling.
- [x] Implement config and cosmology helpers with explicit M200c conventions.
- [x] Implement NFW baseline and SIDM wrapper interface.
- [x] Implement projection module for Sigma, Sigma_bar, DeltaSigma.
- [x] Implement toy forecast metrics and plotting utilities.
- [x] Implement end-to-end benchmark script and output writers.
- [x] Add projection validation test against analytic NFW behavior.
- [x] Run lint, tests, and benchmark script.
- [x] Add completion review section.

## completion review
- Implemented full v1 package stack under `src/sidm_stagev_forecast/` with explicit physical-unit conventions and `M200c` mass-definition consistency.
- Added a stable `parametricSIDM` wrapper contract with informative runtime failure for missing/upstream-changed dependency.
- Added local projection integrals and validated `DeltaSigma` against analytic NFW behavior in `tests/test_projection_nfw.py`.
- Produced end-to-end benchmark outputs (figures, intermediate CSVs, summary CSV) using both backends, including a successful default `parametric` run via the local parametricSIDM compatibility shim.
- Verified code quality with Ruff and pytest (`1 passed`).

## phase 2 progress
- [x] Step 1: Added mass-conservation and profile sanity tests (`tests/test_profiles_sanity.py`).
- [x] Step 2: Added optional NFW reference cross-check tooling (`src/sidm_stagev_forecast/reference_validation.py`, `scripts/run_reference_crosscheck.py`) with generated outputs:
  - `outputs/tables/nfw_reference_crosscheck.csv`
  - `outputs/tables/nfw_reference_backend_status.csv`
  - `outputs/figures/nfw_reference_crosscheck.png`
- [x] Step 3: Precision sweep (`2sigma/3sigma/5sigma`) over selected radial windows with generated outputs:
  - `outputs/tables/precision_sweep.csv`
  - `outputs/figures/dwarf_precision_sweep.png`
  - `outputs/figures/cluster_precision_sweep.png`

## tier 1 ensemble upgrade
- [x] Phase 6: Implemented `ensemble.py` with log-normal halo sampling, concentration relation + scatter, seed reproducibility, and ensemble summary statistics.
- [x] Phase 6 validation: Added `tests/test_ensemble.py` for reproducibility, weight normalization, and scatter recovery.
- [x] Phase 7: Implemented `stacking.py` with log-log interpolation to a common grid and weighted stacking utilities.
- [x] Phase 7 validation: Added `tests/test_stacking.py` including power-law interpolation exactness and single-halo-limit recovery.
- [x] Phase 8: Updated forecast metrics (`forecast.py`) with ensemble-oriented `Delta chi^2` helpers and scenario evaluation.
- [x] Phase 8: Added end-to-end pipeline script `scripts/run_ensemble_forecast.py` and generated Tier-1 outputs:
  - `outputs/figures/ensemble_stacked_delta_sigma.png`
  - `outputs/figures/ensemble_tier1_summary.png`
  - `outputs/tables/ensemble_delta_chi2_summary.csv`
  - `outputs/intermediate/ensemble_halo_catalog.csv`
  - `outputs/intermediate/ensemble_stacked_delta_sigma_profiles.csv`
  - `outputs/intermediate/ensemble_stacked_rho_profiles.csv`
- [x] Validation pass: Full test suite passes (`17 passed`), including projection validation and new ensemble/stacking tests.

## tier 1 ensemble extension (hmf + shmr)
- [x] Added unified ensemble API `generate_ensemble(mode, config_dict)` in `src/sidm_stagev_forecast/ensemble.py`.
- [x] Implemented `mode="HMF"` with:
  - `dn/dM`-based sampling (`hmf_model` or custom callable),
  - optional selection functions (`none`, `threshold`, `logistic`, or callable),
  - fixed or distributed redshift (`uniform` / `discrete`),
  - concentration model + scatter and reproducible seeds.
- [x] Implemented `mode="SHMR"` with:
  - stellar-mass sampling (Tier-1 lognormal),
  - SHMR mapping (`power_law`) plus halo-mass scatter,
  - central-only assumption,
  - concentration model + scatter and reproducible seeds.
- [x] Added config examples in `src/sidm_stagev_forecast/config.py`:
  - `CLUSTER_HMF_ENSEMBLE_CONFIG_EXAMPLE`
  - `DWARF_SHMR_ENSEMBLE_CONFIG_EXAMPLE`
- [x] Added validation tests:
  - HMF distribution + convergence + single-halo limit (`tests/test_ensemble_hmf.py`)
  - SHMR distribution + convergence + single-halo limit (`tests/test_ensemble_shmr.py`)
  - Example-config execution check (`tests/test_ensemble.py`)
- [x] Independent framework validation completed before integration:
  - `pytest tests/test_ensemble.py tests/test_ensemble_hmf.py tests/test_ensemble_shmr.py -q` -> pass
- [x] Integration with projection/stacking validated:
  - HMF smoke run through `scripts/run_ensemble_forecast.py` (surrogate backend, temp output)
  - SHMR smoke run through `scripts/run_ensemble_forecast.py` (surrogate backend, temp output)
- [x] Final regression:
  - full suite `24 passed`.

## yaml-configured ensembles (cluster + dwarf)
- [x] Added YAML parser/normalizer module: `src/sidm_stagev_forecast/ensemble_yaml.py`.
- [x] Added support for config-driven routing in `scripts/run_ensemble_forecast.py` via `--config-path`:
  - `mode: HMF` -> HMF sampler
  - `mode: SHMR` -> SHMR sampler
- [x] Added YAML-required sampling features:
  - Gaussian redshift distribution with truncation
  - Logistic selection mapping
  - Concentration scatter handling
  - Configurable SIDM grids and projection bins
- [x] Added optional HMF/concentration model hooks:
  - `Tinker08` requested in YAML uses optional `colossus` if available; otherwise deterministic fallback proxy with explicit runtime warning.
  - `DiemerJoyce19` concentration uses optional `colossus` if available; otherwise fallback to local relation with warning.
- [x] Added YAML validation tests (`tests/test_ensemble_yaml.py`) covering:
  - mode routing
  - cluster/dwarf median-mass targets
  - projection bin consistency
  - reproducibility with fixed seed
  - single-halo limit (`N_halos=1`)
- [x] Executed pipelines separately before combining:
  - `docs/cluster_ensemble_config.yaml` run complete
  - `docs/dwarf_ensemble_config.yaml` run complete
  - combined comparison table/figure produced afterward.
