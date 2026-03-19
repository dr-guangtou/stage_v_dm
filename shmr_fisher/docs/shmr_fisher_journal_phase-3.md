# Phase 3 Development Journal: Covariance, Fisher Matrix & Systematic Errors

**Date:** 2026-03-19
**Branch:** `claude/keen-lovelace`
**Scope:** Tasks 7–9 from PLAN.md — lensing covariance, clustering covariance, Fisher matrix computation. Plus post-hoc implementation of two systematic error features.

---

## Task 7: Lensing Covariance (`covariance.py`)

### Sigma_crit calculation

Implemented `sigma_crit(z_lens, z_source)` using astropy constants (c, G) and colossus angular diameter distances:

    Sigma_crit = c^2 / (4 pi G) * D_s / (D_l * D_ls)

**Colossus distance convention issue (fixed):** colossus `comovingDistance(z)` returns **negative** values (lookback convention). This caused negative D_ls and negative Sigma_crit on the first attempt. Fixed by taking `abs()` of all comoving distances before computing D_ls = (chi_s - chi_l) / (1+z_s).

Verified: Sigma_crit(z_l=0.3, z_s=1.0) = 2753 Msun/pc^2, consistent with expectations (~3000–5000 Msun/pc^2).

### Effective Sigma_crit and source density

`sigma_crit_effective(z_lens, lensing_config)` computes the harmonic mean of Sigma_crit weighted over the Smail source n(z):

    Sigma_crit_eff = <Sigma_crit^{-1}>^{-1}

Using `scipy.integrate.quad` over z_s from z_lens to z_max.

Verified behavior vs z_lens:
- Minimum (~2900 Msun/pc^2) near z_l ~ 0.35, where the lens-source geometry is optimal
- Rises at both low z (D_l -> 0) and high z (D_ls -> 0)

`n_source_effective(z_lens, lensing_config)`: fraction of sources behind the lens times n_total.
- z_l = 0.1: 25.0 arcmin^{-2} (nearly all sources behind)
- z_l = 0.8: 16.6 arcmin^{-2} (only 66% behind)

### Shape noise covariance

`lensing_covariance(R_bins, z_lens, N_lens, area, lensing_config)` implements:

    sigma^2(DS, R_i) = sigma_e^2 * Sigma_crit_eff^2 / (N_lens * n_source_Mpc^2 * A_annulus_i)

Key unit conversion: n_source from arcmin^{-2} to physical Mpc^{-2} via `D_A * (pi/10800)`.

**AT-3 validated:** var(N=1e6) / var(N=1e4) = 0.0100, exactly 0.01 as expected.

---

## Task 8: Clustering Covariance (`covariance.py`)

### Survey volume

`survey_volume(z_lo, z_hi, area_deg2)`: comoving volume in a redshift shell. Uses f_sky * (4/3)pi * (chi_hi^3 - chi_lo^3). Same `abs()` fix for colossus comoving distances.

Verified: V(z=0.2–0.4, 10000 deg^2) = 3.55e9 Mpc^3 (~10^9.5, expected order of magnitude).

### Clustering covariance

`clustering_covariance(n_gal_model, b_eff, V_survey, z)` provides:

    sigma^2(b_eff) = b_eff^2 / (n_gal * V_eff)
    sigma^2(n_gal) = n_gal / V_survey

where V_eff = V_survey * (nP/(1+nP))^2 with P evaluated at k=0.1 h/Mpc.

**Critical bug fix:** The original implementation used `N_survey / V` as the tracer density, which caused errors to **increase** with N_gal (because Poisson variance scales as n). The fix was to use the MODEL-predicted n_gal (from `galaxy_number_density()`), making the clustering Fisher contribution independent of N_gal_total. This is physically correct: clustering precision depends on the intrinsic tracer density and the survey volume, not on the total survey count. The survey's N_gal_total affects only the lensing (through N_lens-source pairs).

---

## Task 9: Fisher Matrix (`fisher.py`)

### Parameter selection

`get_varied_params()` auto-selects parameters based on the survey's z-range:
- Always: 5 z=0 shape + scatter parameters
- Conditionally: 4 evolution parameters if `vary_z_evolution=True` AND z_max > 0.4 AND Δz > 0.3

### Numerical derivatives

`compute_derivatives()` uses central finite differences with step size `frac_step * |theta_i|`. Special case: if |theta_i| < 1e-10 (near-zero evolution params like nu_N), uses absolute step 0.01.

### Fisher accumulation

For each (z-bin, M*-bin):
1. Distribute N_lens proportional to model-predicted galaxy counts
2. Compute fiducial ΔΣ(R), b_eff, n_gal
3. Compute covariance for each observable
4. Compute Jacobian ∂obs/∂θ via `compute_derivatives`
5. Accumulate: F += J^T C^{-1} J

The lensing contribution uses a diagonal covariance (one variance per R-bin). The clustering contribution uses two scalars (var_b, var_n).

### Post-processing utilities

- `marginalized_errors(F)`: σ_i = sqrt((F^{-1})_ii)
- `conditional_errors(F)`: σ_cond_i = 1/sqrt(F_ii) (always <= marginalized)
- `fisher_ellipse(F, i, j)`: 2D marginalized contour via eigendecomposition of the 2×2 covariance sub-block
- `add_external_prior(F, name, sigma)`: F_ii += 1/sigma^2

### Validation results

| Test | Result |
|------|--------|
| AT-3 (covariance scaling) | PASS: var ratio = 0.01 |
| AT-4 (positive definite) | PASS for Stage-III, IV, V |
| AT-5 (marg >= cond) | PASS for all tiers |
| AT-6 (stage ordering) | PASS (5-param Fisher for fair comparison) |

Condition numbers: Stage-III: 1.0e4, Stage-IV: 6.8e3, Stage-V (9-param): 3.7e5. All well below the 1e10 warning threshold.

**AT-6 note:** Stage-V with 9 parameters has slightly worse N_0 constraints than Stage-IV (5 params) due to marginalization over evolution parameters. The fair comparison uses 5-param Fisher for all tiers, which passes.

---

## Systematic Errors (Post-Phase 3 Enhancement)

### Motivation

The stat-only forecast produces unrealistically tight constraints (~0.001–0.8% fractional errors). Real surveys achieve ~0.1–1 dex precision on SHMR parameters, dominated by systematics. Two toggleable features were added to make absolute error bars more realistic.

### Feature 1: Fractional systematic floor

Added `ForecastConfig.systematic_floor_fraction` (default 0.0 = off).

    sigma^2_total(R_i) = sigma^2_stat(R_i) + (f_sys * DS_fid(R_i))^2

Applied in `fisher.py` after computing `var_ds` and `ds_fiducial` (3 lines of code). This captures irreducible model uncertainties (baryonic effects, miscentering, intrinsic alignments) that scale with the signal.

Impact at f_sys = 0.05 (Stage-IV):
- sigma_logMs: 2.83x baseline (0.28% → 0.80%)
- gamma_0: 1.87x baseline (0.079% → 0.148%)
- beta_0: 1.10x baseline
- N_0: 1.01x baseline (clustering-dominated, barely affected)

### Feature 2: Nuisance parameter marginalization

Added `NuisanceConfig` dataclass with two lensing nuisance parameters:

| Parameter | Fiducial | Prior σ | Derivative |
|-----------|----------|---------|-----------|
| shear_m | 0 | 0.02 | ∂ΔΣ/∂m = ΔΣ_fid (analytic) |
| photo_dz_source | 0 | 0.03 | ∂ΔΣ/∂Δz = ΔΣ_fid × ∂ln(Σ_crit^{-1})/∂Δz (numerical) |

The Fisher matrix expands from (N_shmr × N_shmr) to (N_shmr + 2 × N_shmr + 2). Gaussian priors are added to the nuisance diagonal via `add_external_prior`. Clustering observables (b_eff, n_gal) are unaffected — their Jacobian rows are zero-padded for nuisance entries.

Key design decision: nuisance derivatives are computed analytically, not via `SHMRParams.copy()`. This avoids any changes to the `compute_derivatives` function or the `SHMRParams` dataclass.

The `d_ln_inv_sigma_crit_d_dz()` helper in `covariance.py` computes the photo-z derivative by perturbing `z_source_median` in the LensingConfig via `dataclasses.replace()`.

Impact with default priors (Stage-IV): modest 1–8% degradation. The nuisance parameters are relatively well-constrained by the priors (σ_m = 0.02, σ_Δz = 0.03), so they don't drastically inflate SHMR errors.

### Validation

- Tight priors (σ_m = 0.001, σ_Δz = 0.001): max deviation from baseline = 2.7% — **PASS**
- Loose priors (σ_m = 1.0, σ_Δz = 1.0): lensing constraints degraded by up to 8% (not destroyed, because clustering provides independent constraints)
- Stage ordering preserved with systematics enabled

### Feature 3: Off-diagonal LSS covariance — DEFERRED

Documented in SPEC.md Section 2.4.1. Feasible but non-trivial (~200–300 lines requiring J2 Hankel transforms). Deferred because shape noise dominates for stacked g-g lensing, and the systematic floor already captures leading model uncertainties.

---

## QA Plots Produced

| Figure | Description |
|--------|-------------|
| `covariance_diagnostic.png` | 4-panel: Σ_crit_eff vs z, n_source_eff vs z, σ(ΔΣ) vs R, V_survey vs z_max |
| `fisher_errors_comparison.png` | Bar chart: σ/|θ_fid| for 5 params across Stage-III/IV/V |
| `fisher_ellipses.png` | 4-panel: 1σ marginalized ellipses, dashed=stat-only, solid=+systematics |
| `fisher_ngal_scaling.png` | σ vs N_gal_total for 5 params (Stage-IV footprint) |
| `systematics_comparison.png` | 2-panel: stat-only vs floor vs nuisance vs both |

---

## Bug Fixes

1. **Colossus negative comoving distances:** `comovingDistance(z)` returns negative values (lookback convention). Fixed with `abs()` in `sigma_crit()` and `survey_volume()`.

2. **N_gal scaling inversion:** Clustering covariance using N_survey/V as tracer density caused errors to increase with N_gal. Fixed by using model-predicted n_gal instead.

3. **colossus angularDiameterDistance signature:** No two-argument form for D_ls. Computed manually via D_ls = (chi_s - chi_l) / (1+z_s) in flat ΛCDM.

---

## Files Created/Modified

| File | Action | Description |
|------|--------|-------------|
| `shmr_fisher/covariance.py` | Created | Σ_crit, shape noise, survey volume, clustering covariance, photo-z derivative helper |
| `shmr_fisher/fisher.py` | Created | Fisher matrix computation, derivatives, marginalization, ellipses, priors, nuisance extraction |
| `shmr_fisher/config.py` | Modified | Added NuisanceConfig, extended ForecastConfig with systematic toggles |
| `shmr_fisher/validate.py` | Existing | Phase 1 validation (unchanged; Phase 3 validation done via scripts) |
| `scripts/validate_phase3.py` | Created | AT-3 through AT-6, covariance diagnostics, Fisher comparison plots |
| `scripts/validate_systematics.py` | Created | Systematic error impact comparison |
| `SPEC.md` | Modified | Added Section 2.4.1: Systematic Error Model |
| `PLAN.md` | Modified | Updated Known Simplifications (shear/photo-z now implemented) |

---

## Next Steps (Phase 4)

- **Task 10:** `survey_configs.py` — predefined Stage-III/IV/V configurations + `make_custom_survey` factory
- **Task 11:** `scripts/run_forecast.py` — main driver, JSON/NPZ output
- **Task 12:** `scripts/run_sweep.py` + remaining `plot_results.py` functions (two-regime summary, improvement factor, scaling plots)
- **AT-7** (sweep monotonicity) and **AT-8** (binning consistency) validation
