# SIDM Stage-V Forecast Journal — Phase 3 (Tier-2 Hybrid Outskirts)

Date: 2026-03-19  
Worktree: `sidm_dsigma/worktrees/tier2_hybrid/sidm_dsigma`  
Branch: `codex/tier2-hybrid-dk14`

## Objective

Extend the Tier-1 ensemble weak-lensing forecast to an optional Tier-2 mode:
- preserve Tier-1 behavior by default,
- attach DK14-like outskirts to inner SIDM profiles,
- compare Tier-1 vs Tier-2 at single-halo and stacked-ensemble levels.

## Implementation Summary

### 1. New Tier-2 modules

- Added `src/sidm_stagev_forecast/outer_profiles.py`:
  - DK14-inspired profile structure:
    - inner baseline term,
    - transition/truncation term,
    - outer infall-like term.
  - returns diagnostic components:
    - `rho_total_msun_kpc3`,
    - `rho_inner_baseline_msun_kpc3`,
    - `rho_outer_term_msun_kpc3`,
    - `f_trans`,
    - metadata.

- Added `src/sidm_stagev_forecast/stitch.py`:
  - stitch-radius resolver with modes:
    - `fraction_r200m`,
    - `fraction_r200c`,
    - `fixed_kpc`.
  - smooth logistic blend in log-radius and log-density:
    - method: `logistic_logrho_blend`.

### 2. Core pipeline integration

- Updated `src/sidm_stagev_forecast/cosmology.py`:
  - added `omega_m_z`, `rho_mean_matter_z`, and `rdelta(..., definition='200m')`.
  - retained `M200c` as baseline mass-definition convention.

- Updated `src/sidm_stagev_forecast/profiles.py`:
  - added `build_hybrid_sidm_profile(...)`:
    - builds CDM inner, SIDM inner, DK14-like reference, and stitched hybrid outputs.

- Updated `src/sidm_stagev_forecast/config.py` and `ensemble_yaml.py`:
  - added optional Tier-2 config keys and defaults.

- Updated `scripts/run_ensemble_forecast.py`:
  - Tier-1 path preserved.
  - Tier-2 path activated only if `tier2.enabled=true`.
  - outputs Tier-1 and Tier-2 stacked `DeltaSigma` and `rho` for direct comparison.
  - writes Tier-aware summary rows (`tier1`, `tier2`) in `*_delta_chi2_summary.csv`.

### 3. Robustness updates

- Updated `src/sidm_stagev_forecast/ensemble.py`:
  - runtime fallback when `colossus` calculations fail (e.g., cache write restrictions in sandbox),
  - fallback to local deterministic proxy relation with explicit warnings.

## Validation and Testing

### New tests

- `tests/test_outer_profiles.py`
  - positivity/finiteness checks,
  - outer steepening diagnostic in log-slope.

- `tests/test_stitch.py`
  - inner/outer asymptotic behavior for stitched profiles,
  - finite, positive hybrid output.

- `tests/test_ensemble_yaml.py`
  - Tier-2 parsing assertions,
  - default-disabled behavior when `tier2` block is absent.

### Full suite

- `python -m pytest -q` -> `32 passed`.

### Tier-1 validation requirements rechecked

From generated validation tables:
- cluster median mass: `2.4888e14 Msun` (pass target near `3e14`),
- dwarf median mass: `1.0258e10 Msun` (pass target near `1e10`),
- radial bins, reproducibility, and `N_halos=1` checks: pass for both modes.

## Output Artifacts Generated

New Tier-2 figures (cluster + dwarf):
- `*_single_halo_tier2_diagnostic.png`
- `*_ensemble_tier1_vs_tier2_stacked.png`
- `*_ensemble_tier1_vs_tier2_delta_chi2.png`

Updated machine-readable outputs:
- `outputs/intermediate/*_ensemble_stacked_delta_sigma_profiles.csv`
- `outputs/intermediate/*_ensemble_stacked_rho_profiles.csv`
- `outputs/tables/*_ensemble_delta_chi2_summary.csv`

## Scientific Framing Notes

- Tier-2 is explicitly documented as a hybrid approximation:
  - inner SIDM modification + DK14-like outer attachment.
- It is not a self-consistent SIDM splashback simulation.
- No 2-halo term, baryons, subhalo dynamics, or splashback recalibration added in this stage.

## Detailed Development Timeline

### Step A — New worktree and branch isolation

- Created dedicated worktree and feature branch:
  - worktree: `sidm_dsigma/worktrees/tier2_hybrid/sidm_dsigma`
  - branch: `codex/tier2-hybrid-dk14`
- Kept pre-existing root-worktree modifications isolated from Tier-2 implementation.

### Step B — Tier-2 profile infrastructure

- Implemented DK14-like outskirts model with explicit components:
  - baseline inner (NFW proxy),
  - transition truncation (`f_trans`),
  - outer infall-like power-law term.
- Implemented stitch utility in log-space to reduce projection kinks:
  - logistic window in log radius,
  - logarithmic density blending.

### Step C — Mass-definition and matching-scale controls

- Preserved `M200c` as core mass convention to avoid Tier-1 regression risk.
- Added derived `R200m` support for matching-scale control (`r_match_mode=fraction_r200m`).
- Added `fraction_r200c` and `fixed_kpc` matching options for reproducibility and ablation tests.

### Step D — Pipeline and config integration

- Extended YAML parsing with optional `tier2` block and regime overrides.
- Added Tier-2 runtime path in ensemble script:
  - compute Tier-1 and Tier-2 stacks in one run when enabled,
  - write tier-labeled summary metrics,
  - generate Tier-1-vs-Tier-2 comparison figures.
- Kept legacy behavior unchanged when Tier-2 is disabled.

### Step E — Runtime robustness hardening

- Found runtime `colossus` failures due unwritable cache directories in sandbox.
- Added runtime exception fallbacks:
  - Tinker08 HMF -> power-law proxy fallback,
  - DiemerJoyce19 concentration -> maccio fallback.
- Preserved deterministic execution while still allowing optional `colossus` usage where available.

### Step F — Validation execution

- Ran targeted new tests for Tier-2 modules and parser logic.
- Ran full project suite:
  - `python -m pytest -q` -> `32 passed`.
- Ran end-to-end script runs (surrogate SIDM backend) for both:
  - `docs/cluster_ensemble_config.yaml`
  - `docs/dwarf_ensemble_config.yaml`
- Reconfirmed acceptance checks from output tables:
  - cluster median mass in expected range (near `3e14 Msun`),
  - dwarf median mass in expected range (near `1e10 Msun`),
  - radial bins match configuration,
  - reproducibility and single-halo limit pass.

## Key Technical Decisions and Rationale

1. **Log-density blending over linear blending**
   - Chosen to suppress stitching-induced oscillations in projected `DeltaSigma`.

2. **Tier-2 optional by config, not forced in core API**
   - Avoids behavior changes for existing Tier-1 users and scripts.

3. **Separate regime defaults**
   - Cluster defaults emphasize outskirts leverage (`fraction_r200m`).
   - Dwarf defaults remain conservative (`fraction_r200c`) while still allowing completeness checks.

4. **Runtime fallback for optional external models**
   - Ensures constrained environments still produce complete outputs.

## Final Artifact Checklist (Phase 3)

- [x] `outer_profiles.py`
- [x] `stitch.py`
- [x] Tier-2 config support in YAML loader and config defaults
- [x] Tier-2 integration in ensemble forecast runner
- [x] Single-halo Tier-2 diagnostics (cluster + dwarf)
- [x] Stacked Tier-1 vs Tier-2 comparison figures
- [x] Tier-aware `Delta chi^2` summary tables
- [x] New tests for Tier-2 components and parser behavior
- [x] Updated `SPEC.md`, `PLAN.md`, `README.md`, `docs/todo.md`, `docs/lessons.md`
- [x] Updated `outputs/figures/CAPTION.md` and `outputs/INVENTORY.md`
