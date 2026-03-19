# Phase 1 Development Journal: Foundation

**Date:** 2026-03-19
**Branch:** `claude/keen-lovelace`
**Scope:** Tasks 1–3 from PLAN.md — environment setup, configuration dataclasses, and SHMR model implementation.

---

## Task 1: Environment & Package Structure

### What was done

Created the `shmr_fisher/` Python package with the directory layout specified in CLAUDE.md:

```
shmr_fisher/
├── shmr_fisher/
│   ├── __init__.py
│   ├── config.py
│   ├── shmr_model.py
│   ├── validate.py
│   └── plot_results.py
├── scripts/          (empty, for Phase 4)
├── outputs/
│   └── shmr_validation.pdf
└── pyproject.toml
```

### Dependencies

Installed via `uv add`:
- **colossus** — cosmology, HMF, halo profiles (Diemer 2018)
- **numpy** (v2.x — note: `np.trapz` renamed to `np.trapezoid`)
- **scipy** — root-finding (`brentq`), special functions
- **matplotlib** — validation plots
- **astropy** — units/constants (needed in later phases for Sigma_crit)

### Cosmology initialization

`shmr_fisher/__init__.py` calls `cosmology.setCosmology('planck18')` at import time. Verified that colossus loads Planck18 parameters correctly.

---

## Task 2: Configuration Module (`config.py`)

### Four dataclasses implemented

#### `SHMRParams`
Stores the 9-parameter Moster+2013 SHMR model:
- 4 z=0 shape parameters: `log_M1_0` (11.59), `N_0` (0.0351), `beta_0` (1.376), `gamma_0` (0.608)
- 4 evolution parameters: `nu_M1` (1.195), `nu_N` (−0.0247), `nu_beta` (−0.826), `nu_gamma` (0.329)
- 1 scatter parameter: `sigma_logMs` (0.15 dex)

All fiducial values from Moster+2013 Table 1.

Utility methods:
- `to_dict()` → flat dictionary of all 9 parameters. Verified: returns correct keys/values.
- `from_dict(d)` → reconstructs from dictionary. Verified: round-trips perfectly (`from_dict(to_dict()) == original`).
- `copy(**overrides)` → new instance with selective parameter replacement via `dataclasses.replace()`. Verified: `copy(N_0=0.05)` changes only `N_0`, original is unmodified. This is essential for finite-difference derivatives in the Fisher matrix.

#### `SurveyConfig`
Generic spectroscopic survey specification with:
- `area_deg2`, `z_min`, `z_max`, `n_gal_total` — survey footprint
- `log_Mstar_min`, `log_Mstar_max`, `dlog_Mstar` — stellar mass binning
- `log_Mstar_min_func` — optional callable for z-dependent completeness

Key methods:
- `get_stellar_mass_bins(z)` — returns list of `(lo, hi)` tuples in log10(Msun). Verified behavior:
  - Constant completeness: 6 bins from 9.0 to 12.0 in steps of 0.5 dex.
  - z-dependent completeness (`8.5 + 1.0*z`): 7 bins at z=0 (from 8.5), 5 bins at z=1 (from 9.5).
  - Raises `ValueError` if completeness limit >= upper bound.
- `get_n_gal_per_deg2()` — simple derived quantity.

#### `LensingConfig`
Fixed LSST-like parameters: `n_source_per_arcmin2=25.0`, `sigma_e=0.26`, `z_source_median=1.0`.

#### `ForecastConfig`
Controls the Fisher computation: radial bins (0.1–30 Mpc, 10 bins), redshift binning (dz=0.2), derivative step size (1%), halo mass integration grid (10^10–10^15.5, 200 points), and the `vary_z_evolution` flag for automatically fixing evolution parameters for low-z surveys.

---

## Task 3: SHMR Model (`shmr_model.py`)

### Functions implemented

#### `shmr_params_at_z(params, z)`
Evaluates z-dependent parameters via Moster+2013 Eq. 11–14:
```
log10(M1(z)) = log_M1_0 + nu_M1 * z/(1+z)
N(z)         = N_0      + nu_N  * z/(1+z)
beta(z)      = beta_0   + nu_beta * z/(1+z)
gamma(z)     = gamma_0  + nu_gamma * z/(1+z)
```

Verified parameter evolution:

| z   | log M1 | N      | beta  | gamma |
|-----|--------|--------|-------|-------|
| 0.0 | 11.590 | 0.0351 | 1.376 | 0.608 |
| 0.5 | 11.988 | 0.0269 | 1.101 | 0.718 |
| 1.0 | 12.188 | 0.0227 | 0.963 | 0.772 |
| 2.0 | 12.387 | 0.0186 | 0.825 | 0.827 |

Key trends: M1 increases (characteristic halo mass grows), N decreases (star formation efficiency drops), beta decreases (low-mass slope flattens), gamma increases (high-mass slope steepens). All consistent with Moster+2013.

#### `mean_log_Mstar(log_Mh, params, z)`
Forward SHMR implementing Moster+2013 Eq. 2:
```
M*/Mh = 2N * [(Mh/M1)^(-beta) + (Mh/M1)^gamma]^(-1)
```

Accepts both scalar and array inputs (via `np.asarray`).

Verified at z=0:

| log Mh | log M* | M*/Mh  |
|--------|--------|--------|
| 10.0   | 6.658  | 0.0005 |
| 11.0   | 9.006  | 0.0101 |
| 11.59  | 10.135 | 0.0351 |
| 12.0   | 10.535 | 0.0343 |
| 13.0   | 10.988 | 0.0097 |
| 14.0   | 11.381 | 0.0024 |
| 15.0   | 11.773 | 0.0006 |

Note: at log Mh = 11.59 (= M1), the ratio M*/Mh = 2N/(1+1) = N = 0.0351, confirming the normalization is exact.

#### `mean_log_Mh(log_Mstar, params, z)`
Inverse SHMR via `scipy.optimize.brentq` on the interval [9, 16]. The SHMR is monotonically increasing, so the root is unique.

Includes bracket-checking with a descriptive error message if the target M* lies outside the range spanned by the search interval.

#### `phi_Mstar_given_Mh(log_Mstar, log_Mh, params, z)`
Log-normal scatter distribution. Verified:
- PDF peak = 2.658 (matches analytic 1/(sqrt(2pi)*0.15) = 2.660).
- Integral over +/- 0.5 dex (3.3 sigma) = 0.9991 (expected 0.9993).

---

## Validation Results

### AT-1: Round-trip test (PASS)

`mean_log_Mh(mean_log_Mstar(12.0, fiducial, 0), fiducial, 0) = 12.000000`

Error: < 1e-7 dex, well within the 0.01 dex threshold. The brentq solver with `xtol=1e-6` provides machine-level accuracy for this smooth monotonic function.

Round-trip verified at 16 combinations of (z, log Mh): all errors < 1e-7 dex.

| z   | log Mh = 11 | 12      | 13      | 14      |
|-----|-------------|---------|---------|---------|
| 0.0 | 3.2e-10     | 9.8e-08 | 3.6e-08 | 2.5e-09 |
| 0.5 | 7.3e-11     | 8.8e-08 | 2.2e-08 | 1.0e-09 |
| 1.0 | 1.6e-09     | 6.8e-08 | 5.6e-12 | 2.1e-10 |
| 2.0 | 1.2e-09     | 6.2e-09 | 1.4e-09 | 5.6e-11 |

### Peak M*/Mh (PASS)

| z    | Peak M*/Mh | at log Mh |
|------|-----------|-----------|
| 0.00 | 0.0379    | 11.77     |
| 0.25 | 0.0314    | 11.96     |
| 0.50 | 0.0275    | 12.09     |
| 1.00 | 0.0229    | 12.24     |
| 1.50 | 0.0203    | 12.33     |
| 2.00 | 0.0186    | 12.39     |

At z=0: peak = 3.8% at log Mh = 11.77, within the required range (2–6% at log Mh 11.5–12.5).

The peak amplitude decreases monotonically with redshift (star formation efficiency was lower in the past), and the peak location shifts to higher halo mass (M1 grows with z). Both trends are physically expected and consistent with Moster+2013 Fig. 4.

### Validation plot

Saved to `outputs/shmr_validation.pdf`. Shows M*/Mh vs Mh at z = 0, 0.5, 1.0, 2.0.

Key features visible in the plot:
1. Classic double power-law shape at all redshifts.
2. z=0 curve peaks at ~3.8% near log Mh ~ 11.8.
3. Peak amplitude decreases and shifts rightward with increasing z.
4. Low-mass slope (left side) is steepest at z=0 (beta=1.376) and flattens at high z.
5. All curves converge at the high-mass end (log Mh > 14), consistent with the high-mass slope gamma being relatively similar across redshifts.

---

## Design Decisions

1. **`dataclasses.replace()` for `copy()`** — avoids manual field copying and stays in sync if new fields are added. The `copy(**overrides)` pattern will be used extensively for finite-difference derivatives.

2. **`np.asarray` wrapping in `mean_log_Mstar`** — allows the function to accept both scalars and arrays transparently, returning the same type. This is important because the Fisher matrix code will call it with grids, while the inverse function calls it with scalars.

3. **`brentq` for SHMR inversion** — the SHMR is monotonic (guaranteed by the functional form), so bracketing + bisection is robust and fast. The wide bracket [9, 16] covers all physically relevant halo masses.

4. **Validation before proceeding** — ran AT-1 and all SPEC.md criteria before moving to Phase 2. This prevents error accumulation in downstream modules.

---

## Numpy 2.x Compatibility Note

NumPy 2.x (installed via uv) renamed `np.trapz` to `np.trapezoid`. This will need to be accounted for in all integration code going forward. The `validate.py` module and any future halo model integrals should use `np.trapezoid`.

---

## Files Created

| File | Purpose | Lines |
|------|---------|-------|
| `shmr_fisher/__init__.py` | Package init, sets Planck18 cosmology | 18 |
| `shmr_fisher/config.py` | 4 dataclasses (SHMRParams, SurveyConfig, LensingConfig, ForecastConfig) | 196 |
| `shmr_fisher/shmr_model.py` | 4 core functions (shmr_params_at_z, mean_log_Mstar, mean_log_Mh, phi_Mstar_given_Mh) | 173 |
| `shmr_fisher/validate.py` | validate_shmr() with AT-1 + peak checks | 102 |
| `shmr_fisher/plot_results.py` | plot_shmr_validation() | 82 |
| `outputs/shmr_validation.pdf` | Validation figure: M*/Mh vs Mh at 4 redshifts | — |

---

## Next Steps (Phase 2)

Phase 2 will implement `halo_model.py`:
- **Task 4:** HOD functions `N_cen` and `N_sat` from the SHMR + scatter.
- **Task 5:** Lensing signal `delta_sigma_nfw` and `delta_sigma_bin` using colossus NFW profiles. This is where unit conversion bugs are most likely — will follow the fallback strategy in CLAUDE.md.
- **Task 6:** Clustering summary statistics `effective_bias` and `galaxy_number_density`.
