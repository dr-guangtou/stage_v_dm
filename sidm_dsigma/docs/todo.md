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
