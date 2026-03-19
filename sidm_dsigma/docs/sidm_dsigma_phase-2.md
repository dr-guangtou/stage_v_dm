# SIDM DeltaSigma Forecast Journal — Phase 2

Date: 2026-03-19 (MST)  
Worktree: `/Users/mac/Desktop/stage_v_dm/sidm_dsigma/worktrees/ensemble_tier1`  
Branch: `codex/tier1-ensemble-stacking`

## 1. Scope and Goal

Phase 2 upgraded the completed single-fiducial-halo pipeline to a configurable, selection-aware halo-ensemble framework suitable for Tier-1 stacked weak-lensing forecasting.

The implementation target was:

1. Support both ensemble paradigms from `docs/ensemble-spec.md`:
   - HMF-based (halo-mass-selected)
   - SHMR-based (galaxy-selected)
2. Keep a unified API for ensemble generation.
3. Preserve modularity:
   - ensemble sampling
   - profile generation
   - projection to `DeltaSigma(R)`
   - stacking
   - toy `Delta chi^2` forecasting
4. Add YAML-configured runs for:
   - `docs/cluster_ensemble_config.yaml`
   - `docs/dwarf_ensemble_config.yaml`
5. Validate both pipelines separately before combined comparison.

Tier-1 non-goals were kept unchanged:

- no splashback
- no 2-halo term
- no subhalo modeling
- no baryons

## 2. Files Added and Updated

### 2.1 New/extended core modules

- `src/sidm_stagev_forecast/ensemble.py`
  - Introduced unified interface:
    - `generate_ensemble(mode, config_dict)`
  - Added HMF mode features:
    - mass grid sampling from `dn/dM`
    - optional selection (`none`, `threshold`, `logistic`, callable)
    - fixed/distributed redshift sampling (uniform/discrete/gaussian)
    - concentration model + log scatter
    - deterministic seed behavior
  - Added SHMR mode features:
    - stellar-mass sampling (Tier-1 lognormal)
    - SHMR mapping (`power_law`, lightweight `behroozi13` proxy)
    - halo-mass scatter
    - central-only Tier-1 behavior
  - Kept backward compatibility:
    - `sample_halo_ensemble(...)` still available
  - Added optional model routing:
    - `Tinker08` and `DiemerJoyce19` via optional `colossus`
    - deterministic fallback with explicit runtime warning if unavailable

- `src/sidm_stagev_forecast/stacking.py`
  - Added interpolation and weighted stack helpers used by ensemble runs.

- `src/sidm_stagev_forecast/forecast.py`
  - Added stacked-profile metric helpers:
    - `delta_chi2_with_fractional_error`
    - `sigma_separation_from_delta_chi2`
    - `evaluate_stacked_distinguishability`

- `src/sidm_stagev_forecast/plotting.py`
  - Added stacked profile and Tier-1 summary plotting routines.

- `src/sidm_stagev_forecast/projection.py`
  - Updated numerical integration call to `scipy.integrate.trapezoid` for compatibility.

- `src/sidm_stagev_forecast/config.py`
  - Added example configs:
    - `CLUSTER_HMF_ENSEMBLE_CONFIG_EXAMPLE`
    - `DWARF_SHMR_ENSEMBLE_CONFIG_EXAMPLE`

- `src/sidm_stagev_forecast/ensemble_yaml.py` (new)
  - Added YAML parsing and normalization:
    - maps YAML blocks to `generate_ensemble` runtime config
    - normalizes projection and SIDM grids

### 2.2 New/updated execution script

- `scripts/run_ensemble_forecast.py`
  - Added support for:
    - `--ensemble-mode HMF|SHMR`
    - `--config-path` (YAML-driven run)
  - Added runtime config resolver:
    - default config mode vs YAML-config mode
  - Added validation table generation per run:
    - median mass target check
    - radial-bin check
    - reproducibility check (fixed seed)
    - single-halo limit check (`N_halos=1`)
  - Output filenames are mode/config-label scoped:
    - `cluster_*`, `dwarf_*`, etc.

### 2.3 Test additions

- `tests/test_ensemble.py` (extended)
- `tests/test_ensemble_hmf.py` (new)
- `tests/test_ensemble_shmr.py` (new)
- `tests/test_ensemble_yaml.py` (new)
- Existing stacking/forecast tests retained and passing.

### 2.4 Dependency

- Added `pyyaml` to `pyproject.toml` for YAML parsing.

## 3. Validation Workflow

Validation was done incrementally, not only at final integration.

### 3.1 Framework-level validation (independent)

HMF and SHMR samplers were validated before projection/stacking integration:

1. distribution-shape sanity
2. convergence with increasing `N_halos`
3. single-halo limit behavior
4. seed reproducibility

Test run:

- `pytest tests/test_ensemble.py tests/test_ensemble_hmf.py tests/test_ensemble_shmr.py -q`
- pass

### 3.2 YAML parser and routing validation

YAML tests verify:

1. mode routing (`HMF` vs `SHMR`)
2. parsed projection bins match config
3. median masses land in expected Tier-1 ranges
4. reproducibility with fixed seed
5. single-halo limit

Test run:

- `pytest tests/test_ensemble_yaml.py -q`
- pass

### 3.3 Projection + stacking integration smoke runs

Separate smoke runs were executed through the full pipeline for both modes (first with surrogate backend for speed/robustness), then YAML-configured full runs into tracked outputs.

### 3.4 Full regression

Final full suite:

- `pytest -q`
- Result: `28 passed` (with warnings when optional `colossus` unavailable)

## 4. YAML-Driven Pipeline Runs

### 4.1 Cluster config run

Input:

- `docs/cluster_ensemble_config.yaml`

Key characteristics:

- mode: `HMF`
- N_halos: `400`
- redshift model: gaussian (truncated)
- HMF model requested: `Tinker08`
- selection: logistic
- concentration model requested: `DiemerJoyce19`
- projection bins: `R=50–2000 kpc`, `N_R=30`

Validation output:

- `outputs/tables/cluster_ensemble_validation_checks.csv`

Observed median halo mass:

- `2.4887592659980975e14 Msun` (passes target ~`3e14 Msun`)

### 4.2 Dwarf config run

Input:

- `docs/dwarf_ensemble_config.yaml`

Key characteristics:

- mode: `SHMR`
- N_halos: `800`
- redshift model: gaussian (truncated)
- stellar-mass lognormal sample
- SHMR model requested: `Behroozi13` (Tier-1 proxy)
- concentration model requested: `DiemerJoyce19`
- projection bins: `R=1–300 kpc`, `N_R=30`

Validation output:

- `outputs/tables/dwarf_ensemble_validation_checks.csv`

Observed median halo mass:

- `1.025845196086968e10 Msun` (passes target ~`1e10 Msun`)

### 4.3 Combined comparison

After separate runs, a combined summary was generated:

- table:
  - `outputs/tables/cluster_dwarf_ensemble_delta_chi2_summary.csv`
- figure:
  - `outputs/figures/cluster_dwarf_delta_chi2_comparison.png`

## 5. Output Artifacts (Phase 2)

### 5.1 Cluster YAML outputs

- `outputs/intermediate/cluster_ensemble_halo_catalog.csv`
- `outputs/intermediate/cluster_ensemble_stacked_delta_sigma_profiles.csv`
- `outputs/intermediate/cluster_ensemble_stacked_rho_profiles.csv`
- `outputs/tables/cluster_ensemble_delta_chi2_summary.csv`
- `outputs/tables/cluster_ensemble_validation_checks.csv`
- `outputs/figures/cluster_ensemble_stacked_delta_sigma.png`
- `outputs/figures/cluster_ensemble_tier1_summary.png`

### 5.2 Dwarf YAML outputs

- `outputs/intermediate/dwarf_ensemble_halo_catalog.csv`
- `outputs/intermediate/dwarf_ensemble_stacked_delta_sigma_profiles.csv`
- `outputs/intermediate/dwarf_ensemble_stacked_rho_profiles.csv`
- `outputs/tables/dwarf_ensemble_delta_chi2_summary.csv`
- `outputs/tables/dwarf_ensemble_validation_checks.csv`
- `outputs/figures/dwarf_ensemble_stacked_delta_sigma.png`
- `outputs/figures/dwarf_ensemble_tier1_summary.png`

### 5.3 Combined comparison outputs

- `outputs/tables/cluster_dwarf_ensemble_delta_chi2_summary.csv`
- `outputs/figures/cluster_dwarf_delta_chi2_comparison.png`

## 6. Documentation and Reproducibility Updates

Updated during this phase:

- `outputs/figures/CAPTION.md`
  - added clear captions + timestamps for new figures
- `outputs/INVENTORY.md`
  - added entries for new intermediate/table outputs
- `docs/todo.md`
  - appended detailed phase progress and validation checkpoints

## 7. Notes on Optional External HMF/Concentration Backends

The implementation checks optional external support:

- `Tinker08` HMF
- `DiemerJoyce19` concentration

Current environment status:

- `colossus` not installed
- code emits explicit warnings
- deterministic local fallback paths are used

This preserves pipeline functionality and reproducibility while signaling that requested reference models were approximated under Tier-1 constraints.

## 8. Known Caveats / Future Improvements

1. If exact `Tinker08` and `DiemerJoyce19` behavior is required, install and validate against `colossus` (or `pyccl` equivalent) in this environment.
2. SHMR `Behroozi13` branch is a lightweight Tier-1 proxy, not a full publication-grade fit implementation.
3. `sigma_crit^{-2}` weighting and survey-realistic redshift kernels are still optional/future work.
4. Current forecast covariance remains toy diagonal.

## 9. Final Status

Phase 2 objective (flexible Tier-1 ensemble support with YAML-configured cluster/dwarf pipelines) is completed in this worktree.

Both requested configurations now run end-to-end, with required validation checks passing and outputs saved in machine-readable and figure formats.
