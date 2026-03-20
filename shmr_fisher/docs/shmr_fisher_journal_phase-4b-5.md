# Development Journal: Post-Phase 4 Enhancements & Phase 5

**Date:** 2026-03-20
**Branches:** `main`, `feature/fixed-params`, `feature/phase5-smf-sigma-obs`
**Scope:** YAML configuration system, fixed parameter support, dwarf galaxy forecasts, document alignment, stellar mass uncertainty, and SMF data vector.

This journal covers a full day of intensive development that spanned six major work items, roughly in chronological order.

---

## 1. Mass-Dependent Scatter, 2-Halo Term, and Figure Reorganization

**Commit:** `fb91837` (03:20)

This commit consolidated several previously developed features into a clean merge on `main`:

- **Mass-dependent SHMR scatter (Cao & Tinker 2020):** The scatter at fixed halo mass is no longer assumed constant. Instead, it follows a tanh transition:

  ```
  sigma(Mh) = sigma_high + sigma_rise * [1 - tanh((log_Mh - log_Mh_break) / delta)]
  ```

  with fiducial values sigma_high = 0.18 dex (clusters), sigma_rise = 0.10 dex (so dwarfs see ~0.38 dex). When `SHMRParams.use_mass_dependent_scatter = True`, the Fisher matrix varies `scatter_sigma_high` and `scatter_sigma_rise` instead of the single `sigma_logMs`.

- **Analytic 2-halo lensing term:** The lensing prediction now includes large-scale correlated structure:

  ```
  DeltaSigma_total = DeltaSigma_1h(NFW) + b_eff * DeltaSigma_matter(R, z)
  ```

  The matter correlation function is projected numerically via `scipy.integrate.quad` and cached per redshift. This gives physically correct behavior at R > 2 Mpc without the 77x slowdown that the colossus `diemer23` profile would incur.

- **Figure reorganization:** Output figures moved into per-phase subdirectories (`outputs/phase1/`, `outputs/phase2/`, etc.) and scripts updated to generate figures in the correct locations.

---

## 2. YAML-Driven Survey Configuration

**Commit:** `9da1cb9` (07:44)

Replaced the hard-coded Python survey definitions with a flexible YAML configuration system.

### New files:
- `configs/default.yaml` — 4-survey comparison (Stage-III through Stage-V)
- `configs/stage4_vs_stage5.yaml` — 2-survey head-to-head
- `shmr_fisher/config_io.py` — `RunConfig` dataclass and `load_run_config()` parser

### Key design decisions:

1. **Per-run output directories:** Each YAML config produces outputs in `outputs/{run_name}/` where `run_name` is the config filename stem. A copy of the input YAML is saved for reproducibility.

2. **Z-dependent completeness via YAML:** Instead of `eval()` on arbitrary Python, the schema uses `log_Mstar_min_z_dep: {base, slope}` which constructs `f(z) = base + slope * z` safely.

3. **Backward compatibility:** Legacy CLI flags (`--surveys`, `--systematics`) still work. The YAML path is preferred via `--config`.

4. **Dynamic field filtering in `config_io.py`:** The loader filters YAML keys to match valid `ForecastConfig` fields, so adding new config fields (like `sigma_log_Mstar_obs` later) requires zero changes to the loader.

### .gitignore update (commit `89cae0c`):
- Track PNG figures (they are the main scientific results)
- Ignore data files (JSON, NPZ) since they can be regenerated

---

## 3. Fixed Parameters for Focused Science Cases

**Branch:** `feature/fixed-params`
**Commits:** `b8c2074` (08:21), `a25e596` (13:23)
**Merged:** `0787006` (13:24)

### Motivation

For dwarf galaxy forecasts, the high-mass SHMR parameters (`gamma_0`, `log_M1_0`) are irrelevant because dwarf halos (Mh ~ 10^{10-11}) sit deep in the low-mass power-law regime where the SHMR reduces to:

```
M*/Mh ≈ 2N * (Mh/M1)^beta
```

Varying gamma_0 and log_M1_0 wastes Fisher information and can cause numerical issues.

### Implementation

1. **New field in `ForecastConfig`:**
   ```python
   fixed_params: list[str] = field(default_factory=list)
   ```

2. **Updated `get_varied_params()` in `fisher.py`:** Filters out parameters listed in `fixed_params`, with validation that warns on unknown parameter names.

3. **New config: `configs/dwarf_regime.yaml`:**
   ```yaml
   forecast:
     vary_z_evolution: false
     fixed_params: [gamma_0, log_M1_0]
   surveys:
     stage5_dwarf:
       name: "Stage-V Dwarf Focus"
       z_min: 0.05
       z_max: 0.4
       log_Mstar_min: 7.5
   ```

### Dwarf galaxy forecast results

Created a 4-survey dwarf comparison with realistic survey parameters:

| Survey | Area | z-range | N_gal | log M*_min | log M*_max |
|--------|------|---------|-------|-----------|-----------|
| DESI | 14,000 deg^2 | 0.01-0.08 | 100k | 8.0 | 9.0 |
| DESI-II | 5,000 deg^2 | 0.01-0.15 | 1M | 7.5 | 9.0 |
| Spec-S5 | 5,000 deg^2 | 0.01-0.20 | 12M | 7.0 | 9.0 |
| Stage-V Dwarf Mania | 8,000 deg^2 | 0.01-0.30 | 40M | 7.0 | 9.0 |

**Key finding — scatter constraints from dwarfs alone are weak:**

With stellar mass capped at 10^9 (dwarf-only), the best survey achieves only ~82% fractional error on sigma_high and ~74% on sigma_rise. The reason: all dwarf halos sit on the same side of the SHMR peak, so scatter just smears things within one regime without creating distinctive signatures.

**Full-range comparison (removing the M* cap, 3x galaxies):**

Extending to the full stellar mass range (log M* = 7–12) transforms scatter constraints from ~82% to ~0.6% — a 130x improvement. The lever arm across both sides of the SHMR peak (M_h ~ 10^12) is essential for constraining scatter.

Output figures:
- `outputs/dwarf_regime/dwarf_error_comparison.png`
- `outputs/dwarf_regime/dwarf_delta_sigma.png`
- `outputs/dwarf_regime/fullrange_error_comparison.png`

---

## 4. Document Alignment Check

**Commit:** `3cac96e` (13:58) — on `main`

Before implementing Phase 5, all four project documents were carefully updated to reflect the current codebase and plan the new features:

### PLAN.md
- Added **Phase 5** (Tasks 13–16): stellar mass uncertainty + SMF data vector
- Marked completed stretch goals: nuisance params, mass-dependent scatter, YAML configs, fixed_params, 2-halo term
- Updated known simplifications list

### DESIGN.md
- Added **§9: Stellar Mass Measurement Uncertainty** — mathematical framework for sigma_eff = sqrt(sigma_SHMR^2 + sigma_obs^2), magnitude table, implementation plan
- Added **§10: Stellar Mass Function as a Data Vector** — SMF covariance model (Poisson + cosmic variance), relationship to existing n_gal, Fisher contribution
- Updated data flow diagram (§7) to show sigma_eff and SMF
- Updated limitations (§11) with new items

### SPEC.md
- Added `sigma_log_Mstar_obs`, `include_smf`, `fixed_params` to `ForecastConfig` spec
- Added `smf_covariance()` function specification
- Updated `get_varied_params()` and `compute_fisher_matrix()` docstrings
- Added acceptance tests AT-9 (sigma_obs effect) and AT-10 (SMF covariance consistency)

### README.md
- Updated feature list with SMF, stellar mass uncertainty, fixed parameters
- Updated systematic error options section

---

## 5. Phase 5 Implementation: Stellar Mass Uncertainty & SMF Data Vector

**Branch:** `feature/phase5-smf-sigma-obs`
**Commit:** `bf38b41` (14:21)

### 5.1 Stellar mass measurement uncertainty (sigma_obs)

**The insight:** When we bin galaxies by *observed* (estimated) stellar mass rather than true stellar mass, the HOD per observed bin is equivalent to using an effective scatter:

```
sigma_eff(Mh) = sqrt(sigma_SHMR(Mh)^2 + sigma_obs^2)
```

This is exact because both the SHMR scatter and the measurement error are Gaussian in log-space, and the convolution of two Gaussians is another Gaussian.

**Changes:**
- `config.py`: Added `sigma_log_Mstar_obs: float = 0.0` to `ForecastConfig`
- `halo_model.py`: All HOD functions (`n_cen`, `n_sat`, `delta_sigma_bin`, `effective_bias`, `galaxy_number_density`) accept `sigma_obs` parameter. When > 0, scatter is replaced with sigma_eff before computing erf integrals.
- `fisher.py`: Passes `sigma_obs = forecast_config.sigma_log_Mstar_obs` to all halo model calls and derivative lambda functions.

**Validation result:** For the Cao & Tinker scatter with sigma_obs = 0.1 dex:
- At high Mh (sigma_SHMR = 0.18): sigma_eff = 0.206, +15% increase
- At low Mh (sigma_SHMR = 0.38): sigma_eff = 0.393, +3% increase
- Impact on dwarf forecasts: negligible (<2% change in constraints)

### 5.2 SMF data vector with cosmic variance

**The problem:** The old pipeline used Poisson-only noise for the n_gal clustering observable: `var(n) = n / V`. This was artificially optimistic for large surveys where cosmic variance dominates.

**The fix:** New `smf_covariance()` function in `covariance.py`:

```
var(Phi) = Phi / V_survey + b_eff^2 * sigma_m^2(V, z) * Phi^2
```

where `sigma_m` is the matter RMS fluctuation on the survey scale, computed from colossus:

```python
R_survey = (3 * V / (4 * pi))^{1/3}  # effective radius [Mpc]
sigma_m = cosmo.sigma(R_survey * h, z)  # colossus matter variance
```

**Changes:**
- `config.py`: Added `include_smf: bool = True` to `ForecastConfig`
- `covariance.py`: New `smf_covariance()` function
- `fisher.py`: When `include_smf=True`, uses `smf_covariance()` instead of Poisson-only `var_n` for the n_gal Fisher contribution

### 5.3 Validation results — disentangled effects

Ran four configurations for the Stage-V Dwarf Mania survey to isolate each effect:

| Configuration | N_0 | beta_0 | sigma_high | sigma_rise |
|--------------|-----|--------|------------|------------|
| Baseline (old pipeline) | 0.3% | 0.1% | 82% | 74% |
| sigma_obs = 0.1 only | 0.3% | 0.1% | 81% | 72% |
| SMF covariance only | 3.1% | 0.7% | 352% | 317% |
| Both (new default) | 3.2% | 0.7% | 341% | 306% |

**Key findings:**

1. **sigma_obs = 0.1 dex** has almost no effect on dwarf forecasts. This is physically expected: for dwarf halos, sigma_SHMR ≈ 0.38 dex already, so adding 0.1 dex in quadrature gives 0.393 — barely different. The effect would be larger for massive halos where sigma_SHMR ≈ 0.18 dex.

2. **SMF cosmic variance** is the dominant correction. The old pipeline was too optimistic because it treated n_gal as a Poisson-counting observable. In reality, cosmic variance adds a noise floor that doesn't decrease with more galaxies. For the Stage-V Dwarf Mania survey (8,000 deg^2, z=0.01-0.30):
   - Shape parameter errors (N_0, beta_0) increase ~10x
   - Scatter errors increase ~4x

3. **This is the physically correct behavior.** The n_gal observable was providing unrealistically tight constraints. With cosmic variance included, the SMF constraining power is properly bounded by the survey volume and tracer bias.

4. **The lensing and clustering (b_eff) observables are unaffected** by these changes. They remain the dominant source of Fisher information for the SHMR shape parameters.

---

## 6. Current State and Open Items

### Branch topology:
```
main (3cac96e)  ←  Phase 5 docs committed
  └── feature/phase5-smf-sigma-obs (bf38b41)  ←  Implementation, ready to merge
```

### Files modified in Phase 5 (relative to main):
| File | Lines changed | Description |
|------|--------------|-------------|
| `shmr_fisher/config.py` | +12 | New ForecastConfig fields |
| `shmr_fisher/halo_model.py` | +78/-15 | sigma_obs propagation through HOD |
| `shmr_fisher/covariance.py` | +65 | New smf_covariance() function |
| `shmr_fisher/fisher.py` | +19/-2 | sigma_obs passthrough + SMF switch |
| `scripts/run_forecast.py` | +2 | Output new settings |
| `scripts/validate_phase5.py` | +183 (new) | End-to-end validation |

### Backward compatibility:
- `sigma_log_Mstar_obs = 0.0` (default) → identical to previous behavior
- `include_smf = True` (default) → **changes behavior** from Poisson-only to Poisson + cosmic variance. This is intentional: the new default is more physically correct. Set `include_smf = False` to recover old behavior.

### Open items for future work:
- Re-run the full-range forecast (all M* bins, no cap) with SMF covariance to see if scatter constraints remain strong
- Investigate whether the SMF provides complementary information to lensing for scatter constraints when the full M* range is used
- Consider adding sigma_obs as a *varied* parameter in the Fisher matrix (rather than fixed) to forecast how well surveys can self-calibrate their stellar mass uncertainty
- Update dwarf_regime.yaml to include `sigma_log_Mstar_obs: 0.1` as the default for spectroscopic surveys
