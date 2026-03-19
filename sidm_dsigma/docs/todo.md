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

## tier 2 hybrid outskirts extension
- [x] Added optional Tier-2 modules:
  - `src/sidm_stagev_forecast/outer_profiles.py` (DK14-like outer structure)
  - `src/sidm_stagev_forecast/stitch.py` (smooth log-density stitching)
- [x] Extended cosmology helpers with `200m` support for stitch-scale controls while preserving `M200c` as core pipeline convention.
- [x] Added hybrid builder in `profiles.py`:
  - CDM Tier-2: CDM inner + DK14-like outskirts
  - SIDM Tier-2: SIDM inner + DK14-like outskirts
- [x] Extended YAML/config handling with optional `tier2` block and regime defaults.
- [x] Integrated Tier-2 into `scripts/run_ensemble_forecast.py` without changing Tier-1 default behavior (`tier2.enabled=false`).
- [x] Added Tier-2 figures:
  - single-halo diagnostics (cluster + dwarf)
  - stacked Tier-1 vs Tier-2 comparison (cluster + dwarf)
  - `Delta chi^2` Tier-1 vs Tier-2 comparison (cluster + dwarf)
- [x] Added Tier-2 tests:
  - `tests/test_outer_profiles.py`
  - `tests/test_stitch.py`
  - YAML Tier-2 parse/default checks in `tests/test_ensemble_yaml.py`
- [x] Validation status:
  - full test suite pass: `32 passed`
  - cluster median mass check pass (`~2.49e14 Msun`)
  - dwarf median mass check pass (`~1.03e10 Msun`)
  - projection bin, reproducibility, and single-halo checks pass in both modes

## sigma-grid regime update (paper-informed)
- [x] Updated regime-specific SIDM grids in YAML configs:
  - cluster: `[0.05, 0.1, 0.2, 0.3]` cm^2/g
  - dwarf: `[10.0, 20.0, 50.0, 100.0]` cm^2/g
- [x] Updated code defaults in `config.py` to keep benchmark runs consistent with regime intent.
- [x] Updated `README.md` and `SPEC.md` benchmark-grid documentation to match runtime configuration.
- [x] Re-ran cluster and dwarf ensemble forecast scripts and verified output artifacts/captions/inventory updates.
