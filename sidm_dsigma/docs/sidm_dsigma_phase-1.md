# sidm_dsigma Phase 1 Development Log

## 1. Phase objective and scope

Phase 1 implements a minimal, first-order forecasting pipeline for distinguishing CDM from SIDM using projected weak-lensing observables.

Primary target:
- Build an end-to-end, transparent workflow from halo density profiles to a toy distinguishability metric.

Scientific scope for this phase:
- Halo benchmarks only:
  - dwarf: `M200 = 1e10 Msun`
  - cluster: `M200 = 1e14 Msun`
- SIDM parameterization:
  - `sigma_over_m_eff` in `cm^2/g`
- Core signal chain:
  - `rho(r) -> Sigma(R) -> DeltaSigma(R)`
- Toy forecast layer:
  - fractional per-bin errors
  - `Delta chi^2` and simple required precision summary

Explicitly out of scope in Phase 1:
- HOD and galaxy-population realism
- 2-halo term
- miscentering
- survey mask/selection systematics
- realistic covariance modeling
- photo-z and IA propagation
- publication-grade inference machinery

## 2. Inputs that governed implementation

The implementation followed these project constraints and planning documents:
- `README.md`
- `SPEC.md`
- `AGENTS.md`
- `PLAN.md`

Operational mandates applied during development:
- Work executed in a new worktree and non-`main` branch.
- Code/data/docs written in English.
- Explicit unit conventions and mass-definition consistency.
- Local numerical projection implementation (`numpy`/`scipy`), no `cluster_toolkit` runtime dependency.
- Small, testable functions with clear interfaces.

## 3. Environment and branch/worktree setup

Execution workspace:
- `/Users/mac/Desktop/stage_v_dm/sidm_dsigma/worktrees/minimal_forecast/sidm_dsigma`

Branch:
- `codex/minimal-forecast-v1`

Reason for this setup:
- Satisfy non-`main` requirement.
- Isolate Phase 1 development artifacts from the primary working tree.

## 4. Project structure created in Phase 1

### 4.1 Source package

- `src/sidm_stagev_forecast/__init__.py`
- `src/sidm_stagev_forecast/config.py`
- `src/sidm_stagev_forecast/cosmology.py`
- `src/sidm_stagev_forecast/profiles.py`
- `src/sidm_stagev_forecast/projection.py`
- `src/sidm_stagev_forecast/forecast.py`
- `src/sidm_stagev_forecast/plotting.py`
- `src/sidm_stagev_forecast/io.py`

### 4.2 Runner and tests

- `scripts/run_benchmarks.py`
- `tests/test_projection_nfw.py`

### 4.3 Tooling and process docs

- `pyproject.toml`
- `.pre-commit-config.yaml`
- `docs/todo.md`
- `docs/lessons.md`

### 4.4 Runtime outputs

- `outputs/figures/`
- `outputs/tables/`
- `outputs/intermediate/`

### 4.5 External SIDM source checkout

- `third_party/parametricSIDM/` (Daneng Yang repository clone)

## 5. Module-by-module development details

## 5.1 `config.py`

Purpose:
- Centralize benchmark and forecast settings.

Implemented objects:
- `CosmologyConfig`
- `BenchmarkConfig`
- `ForecastConfig`

Default benchmark definitions:
- Dwarf:
  - `m200_msun=1e10`
  - `c200=15`
  - `z=0.3`
  - `sigma_over_m=[0.0, 0.2, 0.5, 1.0, 2.0]`
  - 3D radius range: `0.1-300 kpc`
  - lensing radius range: `3-300 kpc`
- Cluster:
  - `m200_msun=1e14`
  - `c200=4`
  - `z=0.3`
  - `sigma_over_m=[0.0, 0.2, 0.5, 1.0, 2.0]`
  - 3D radius range: `5-5000 kpc`
  - lensing radius range: `30-3000 kpc`

## 5.2 `cosmology.py`

Purpose:
- Keep low-level cosmology helpers explicit and unit-safe.

Conventions enforced:
- Radii in `kpc`
- Mass in `Msun`
- Density in `Msun/kpc^3`
- Default mass definition: `M200c`

Implemented functions:
- `e_z(z, cosmo)`
- `hubble_km_s_kpc(z, cosmo)`
- `rho_crit_z(z, cosmo)`
- `rdelta(mass_msun, z, cosmo, definition="200c")`

`rdelta` behavior:
- Supports `definition="200c"` only.
- Fails fast for any other mass-definition label.

## 5.3 `profiles.py`

Purpose:
- Provide CDM baseline profiles and a stable SIDM wrapper.

### CDM/NFW layer

Implemented:
- `NfwParameters` dataclass
- `nfw_parameters_from_m_c`
- `nfw_profile_from_m_c`
- `nfw_m_enclosed`
- `circular_velocity_km_s`
- `nfw_profile_bundle`

`nfw_profile_bundle` returns:
- `r_kpc`
- `rho_msun_kpc3`
- `m_enclosed_msun`
- `vcirc_km_s`
- metadata with `m200`, `c200`, `z`, and mass-definition tags.

### SIDM wrapper layer

Implemented:
- `_infer_parametric_sidm_function`
- `sidm_profile_from_parametric_model`

Design choices:
- Keep a stable local interface regardless of upstream package changes.
- Permit custom callable injection with `model_options['sidm_callable']`.
- If no callable is supplied, attempt importing `parametricSIDM` and resolving a known function.
- Provide informative runtime failure message when unavailable.

Output contract:
- Same profile bundle keys as CDM (`r`, `rho`, `M(<r)`, `Vc`) plus SIDM metadata.

## 5.4 `projection.py`

Purpose:
- Compute lensing observables from 3D density.

Implemented functions:
- `_loglog_interpolator`
- `sigma_of_R`
- `sigma_bar_of_R`
- `delta_sigma_of_R`
- `analytic_nfw_delta_sigma`

Numerical design details:
- Log-log interpolation in radius-density space.
- Line-of-sight integral uses substitution `r = sqrt(R^2 + z^2)`.
- `Sigma(R)` computed via numerical integration over `z`.
- `Sigma_bar(<R)` computed by cumulative area weighting.
- `DeltaSigma(R) = Sigma_bar - Sigma`.

Stability enhancement added:
- `delta_sigma_of_R` internally extends the projected-radius grid below requested `R_min` before computing `Sigma_bar`, then interpolates back.
- This prevents a biased first-bin `DeltaSigma` estimate from hard-truncating the inner integral.

Validation helper:
- `analytic_nfw_delta_sigma` provides a reference form for NFW.

## 5.5 `forecast.py`

Purpose:
- Compute toy distinguishability metrics.

Implemented:
- `fractional_error_model(R, regime, scenario)`
- `delta_chi2(model, reference, sigma)`
- `required_uniform_fractional_precision(model, reference, target_sigma_significance)`

Assumption set:
- Toy per-bin fractional error curves differ by regime (`dwarf` vs `cluster`) and scenario (`baseline` vs `conservative`).
- `Delta chi^2` uses diagonal per-bin uncertainties.
- Required precision estimate assumes a uniform fractional uncertainty per bin relative to the reference signal.

## 5.6 `plotting.py`

Purpose:
- Generate slide-ready comparison figures.

Implemented plot products:
- `plot_rho_profiles`:
  - left panel: `rho(r)` for CDM and SIDM grid
  - right panel: `rho_sidm / rho_cdm`
- `plot_delta_sigma_profiles`:
  - left panel: `DeltaSigma(R)` with toy `1 sigma` shaded CDM band
  - right panel: `DeltaSigma_sidm / DeltaSigma_cdm`

Implementation note:
- Mathtext labels were adjusted to valid raw-string usage after an escape/parsing issue surfaced during save.

## 5.7 `io.py`

Purpose:
- Keep output handling explicit and deterministic.

Implemented:
- `ensure_output_directories`
- `save_table`

This module standardizes where figures, tables, and debug intermediates are written.

## 5.8 `scripts/run_benchmarks.py`

Purpose:
- End-to-end benchmark execution for dwarf and cluster cases.

Pipeline actions per benchmark:
1. Build 3D and projected radius grids.
2. Generate CDM NFW profile.
3. Project CDM profile to `DeltaSigma`.
4. Loop over SIDM benchmark cross sections.
5. Generate SIDM profile (parametric backend or explicit fallback backend).
6. Project SIDM to `DeltaSigma`.
7. Compute metrics:
   - `Delta chi^2`
   - `SNR = sqrt(Delta chi^2)`
   - required uniform fractional precision for `3 sigma`
   - max absolute ratio shift (%)
8. Save intermediate profile and `DeltaSigma` tables.
9. Save benchmark-level figures.
10. Save summary table.

CLI options:
- `--output-root` (default: `outputs`)
- `--sidm-backend` (`parametric` or `surrogate`)

## 6. parametricSIDM integration details

## 6.1 Initial integration issue

The original wrapper expected a direct importable Python package named `parametricSIDM` with a specific function. The cloned upstream repository does not expose that package API directly.

## 6.2 Compatibility shim added

File:
- `src/parametricSIDM.py`

What it does:
- Exposes `density_profile_from_m200_c200(...)` as a stable callable name expected by local wrapper logic.
- Loads `third_party/parametricSIDM/parametricC4.py` dynamically.
- Converts `(M200c, c200, z, sigma_over_m)` into C4 profile parameters.

Behavior:
- `sigma_over_m <= 0`: return NFW-equivalent profile.
- `sigma_over_m > 0`:
  - compute collapse time scale via upstream `tc`
  - compute `tr = t/tc` with clipping to `[0, 1.1]`
  - evaluate `rhost`, `rst`, `rct`, then `frho`

Configurable override:
- `PARAMETRIC_SIDM_REPO` env var can redirect to a different local checkout path.

## 6.3 Why this approach

- Keeps local interface stable and explicit.
- Avoids editing upstream source files.
- Satisfies requirement to use Daneng Yang's code while preserving project modularity.

## 7. Test and validation work

## 7.1 Projection validation test

File:
- `tests/test_projection_nfw.py`

Validation strategy:
- Generate NFW profile numerically.
- Compute numeric `DeltaSigma` via local projection.
- Compare to analytic NFW `DeltaSigma` behavior.
- Assert maximum fractional deviation below tolerance.

Outcome:
- Test passed after fixing small-radius `Sigma_bar` handling and integration compatibility.

## 7.2 Lint and format

Tooling:
- Ruff for check/format.
- Pre-commit hook config based on Ruff.

Adjustment made:
- Added Ruff exclude for `third_party/parametricSIDM` to avoid linting vendored upstream scripts.

## 7.3 Runtime/compatibility issues encountered and fixed

1. `uv sync` panic in environment
- Symptom: runtime panic in `uv` executable.
- Resolution: proceeded with existing interpreter/tooling for validation commands.

2. NumPy integration method compatibility
- Symptom: `np.trapezoid` unavailable in active NumPy.
- Resolution: switched to `np.trapz`.

3. First-bin `DeltaSigma` mismatch in validation
- Symptom: large discrepancy at smallest `R` due to missing inner contribution in `Sigma_bar` integration.
- Resolution: internal projected-radius extension below requested `R_min` before interpolation back.

4. Matplotlib mathtext parse failure
- Symptom: escaped label strings invalid.
- Resolution: corrected math label strings and raw-string usage.

5. Font cache write warnings
- Symptom: non-writable default cache paths in sandbox.
- Resolution: run scripts with local `XDG_CACHE_HOME` and `MPLCONFIGDIR`.

## 8. Outputs produced in Phase 1

Figures:
- `outputs/figures/dwarf_rho_profiles.png`
- `outputs/figures/dwarf_delta_sigma_profiles.png`
- `outputs/figures/cluster_rho_profiles.png`
- `outputs/figures/cluster_delta_sigma_profiles.png`

Summary table:
- `outputs/tables/benchmark_summary.csv`

Intermediate debug tables:
- profile-level and `DeltaSigma`-level CSVs for each benchmark and each `sigma_over_m`.

## 9. Current benchmark summary (parametric backend run)

The latest generated summary table records the following columns:
- `benchmark`
- `m200_msun`
- `c200`
- `z`
- `sigma_over_m_cm2_g`
- `delta_chi2_baseline`
- `snr_baseline`
- `required_uniform_fraction_for_3sigma`
- `max_abs_delta_sigma_ratio_shift_percent`

Observed pattern in this Phase 1 run:
- CDM rows (`sigma_over_m = 0`) correctly return null distinguishability metrics.
- SIDM rows show increasing distinguishability with cross section in most benchmark cases.
- Cluster benchmark in this run gives stronger toy `Delta chi^2` separation than dwarf for higher cross sections under the adopted error model.

These are workshop-level indicators only and should not be treated as precision survey claims.

## 10. Reproducibility steps

From project root (`sidm_dsigma`):

```bash
PYTHONPATH=src pytest -q
ruff check .
XDG_CACHE_HOME=.cache MPLCONFIGDIR=.mplconfig PYTHONPATH=src python scripts/run_benchmarks.py --sidm-backend parametric
```

Optional fallback backend run:

```bash
XDG_CACHE_HOME=.cache MPLCONFIGDIR=.mplconfig PYTHONPATH=src python scripts/run_benchmarks.py --sidm-backend surrogate
```

## 11. Phase 1 acceptance checklist

Requested deliverable status:
1. Stable wrapper around `parametricSIDM`: implemented.
2. CDM baseline NFW profile generator: implemented.
3. Numerical projector for `Sigma`, `Sigma_bar`, `DeltaSigma`: implemented.
4. Toy forecast module with fractional errors and `Delta chi^2`: implemented.
5. One end-to-end script: implemented (`scripts/run_benchmarks.py`).
6. Slide-ready figures and machine-readable summary tables: implemented and generated.

Validation expectations status:
- Mass-definition consistency: enforced as `M200c` by default and metadata-tagged.
- Projection validation on NFW: implemented in test and passing.
- Intermediate output saving: implemented (`outputs/intermediate`).

## 12. Known limitations and next-phase priorities

Known limitations in this phase:
- Toy diagonal error model only.
- No covariance realism or observational selection effects.
- No 2-halo or baryonic profile components in the main path.
- `parametricSIDM` integration uses local compatibility shim rather than a pip-installable upstream package API.

Recommended Phase 2 priorities:
1. Add explicit unit tests for mass-conservation consistency in profile bundles.
2. Add optional cross-check mode against `colossus` or `pyccl` for standard NFW projection behavior.
3. Add scenario sweeps for required precision targets (`2 sigma`, `3 sigma`, `5 sigma`) and radial subranges.
4. Add a reproducible environment lock path once `uv` behavior in this environment is stable.

## 13. Phase conclusion

Phase 1 successfully established a clean, modular, and runnable minimal forecasting pipeline aligned with project constraints.

The implementation now supports:
- benchmark CDM/SIDM profile generation,
- local lensing projection,
- toy distinguishability quantification,
- and persistent figure/table outputs for workshop-level interpretation.

This phase is complete as a first-order forecasting baseline and is ready for controlled scientific extensions.
