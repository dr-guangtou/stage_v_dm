# DESIGN.md — Forecast Pipeline Design

This document describes the end-to-end logic of the `shmr_fisher` Fisher
forecast pipeline: how the parameterized SHMR model is converted into
observable predictions, how noise is modeled, and how the Fisher matrix
is assembled.  Every equation maps to a specific function in the codebase;
cross-references are given as `module.function()`.

---

## 1. The SHMR Model

**Module:** `shmr_model.py`
**Configuration:** `config.SHMRParams`

### 1.1 Mean Relation (Moster+2013)

The stellar-to-halo mass ratio is a double power-law in halo mass:

```
M*(Mh, z) / Mh = 2 N(z) / [(Mh/M1(z))^(-β(z)) + (Mh/M1(z))^γ(z)]
```

so that `log₁₀(M*) = log₁₀(Mh) + log₁₀(M*/Mh)`.

Each shape parameter evolves linearly in the variable `z/(1+z)`:

| Parameter | z=0 value | Evolution | Formula |
|-----------|-----------|-----------|---------|
| `log₁₀ M1` | `log_M1_0 = 11.59` | `nu_M1 = 1.195` | `log_M1_0 + nu_M1 · z/(1+z)` |
| `N`       | `N_0 = 0.0351`     | `nu_N = −0.0247` | `N_0 + nu_N · z/(1+z)` |
| `β`       | `beta_0 = 1.376`   | `nu_beta = −0.826` | `beta_0 + nu_beta · z/(1+z)` |
| `γ`       | `gamma_0 = 0.608`  | `nu_gamma = 0.329` | `gamma_0 + nu_gamma · z/(1+z)` |

Fiducial values are from Moster+2013 Table 1.

**Code path:**
- `shmr_params_at_z(params, z)` → returns `(log_M1, N, β, γ)` at redshift z.
- `mean_log_Mstar(log_Mh, params, z)` → forward relation, returns `log₁₀(M*)`.
- `mean_log_Mh(log_Mstar, params, z)` → inverse relation via `scipy.optimize.brentq`
  root-finding on the interval `log Mh ∈ [9, 16]`.  The SHMR is monotonic,
  so the root is unique.


### 1.2 Scatter at Fixed Halo Mass

At fixed halo mass, the stellar mass is log-normally distributed:

```
P(log M* | Mh) = (2π σ²)^(-1/2) exp[−(log M* − ⟨log M*⟩)² / (2σ²)]
```

The pipeline supports two scatter models, selected by the boolean flag
`SHMRParams.use_mass_dependent_scatter`.

#### Constant scatter (default, `use_mass_dependent_scatter = False`)

```
σ(Mh) = sigma_logMs = 0.15 dex    (constant for all Mh)
```

A single parameter `sigma_logMs` is varied in the Fisher matrix.

#### Mass-dependent scatter (Cao & Tinker 2020)

```
σ(Mh) = σ_high + σ_rise · [1 − tanh((log Mh − log Mh_break) / δ)]
```

| Parameter | Fiducial | Meaning |
|-----------|----------|---------|
| `scatter_sigma_high` | 0.18 dex | Asymptotic scatter at high Mh (Mh ≫ Mh_break) |
| `scatter_sigma_rise` | 0.10 dex | Half-amplitude of the rise toward low Mh |
| `scatter_log_Mh_break` | 12.0 | Transition halo mass [log₁₀(M☉)] |
| `scatter_delta` | 0.4 dex | Width of the tanh transition |

Asymptotic behavior:
- High-mass halos (log Mh ≫ 12): σ → σ_high = 0.18 dex
- Low-mass halos (log Mh ≪ 12): σ → σ_high + 2 σ_rise = 0.38 dex

When this mode is active, the Fisher matrix varies `scatter_sigma_high` and
`scatter_sigma_rise` (the two parameters the forecast is most sensitive to).
The shape parameters `scatter_log_Mh_break` and `scatter_delta` are held
fixed.

No redshift evolution of the scatter is assumed (consistent with
Behroozi+2019, who find only weak evolution).

**Code path:**
- `scatter_at_Mh(log_Mh, params)` → returns σ as a float or array.
- `phi_Mstar_given_Mh(log_Mstar, log_Mh, params, z)` → the full PDF,
  calls `scatter_at_Mh` internally.


---

## 2. From SHMR to Halo Occupation

**Module:** `halo_model.py`, functions `n_cen()` and `n_sat()`

The SHMR does not directly predict "how many galaxies are in a stellar mass
bin."  That requires converting the continuous SHMR + scatter into a **halo
occupation distribution (HOD)** — the mean number of galaxies per halo that
fall in a given stellar mass bin `[log M*_lo, log M*_hi]`.

### 2.1 Central Occupation

At each halo mass Mh, the probability that the central galaxy falls in the
bin is the integral of the log-normal scatter over the bin edges:

```
N_cen(Mh) = ½ [erf(x_hi) − erf(x_lo)]
```

where

```
x = (log M*_edge − ⟨log M*⟩(Mh, z)) / (√2 · σ(Mh))
```

With mass-dependent scatter, `σ(Mh)` is an array that varies across the
halo mass grid.  The erf expressions are fully vectorized — no loops.

**Physical effect of mass-dependent scatter on N_cen:**
With constant scatter (σ = 0.15), the occupation transitions sharply from
0 to 1 as Mh crosses the characteristic mass for that stellar mass bin.
With mass-dependent scatter, low-mass halos (where σ ~ 0.38 dex) have a
much broader transition — more low-mass halos contribute to a given stellar
mass bin, and the occupation never reaches exactly 1 at intermediate masses.
At high Mh (σ ~ 0.18 dex), the transition is sharper than the constant-
scatter case.


### 2.2 Satellite Occupation

Satellites follow a simplified Leauthaud+2011 / Zheng+2005 power-law
prescription:

```
N_sat(Mh) = N_cen_thresh(Mh) · (Mh / M1_sat)^α_sat
```

where:
- `N_cen_thresh(Mh) = ½ [1 + erf((⟨log M*⟩ − log M*_lo) / (√2 · σ(Mh)))]`
  is the threshold central occupation (all galaxies above `log M*_lo`).
  This also uses mass-dependent σ(Mh).
- `M_min` is the halo mass where `N_cen_thresh = 0.5`, found via the
  inverse SHMR: `M_min = 10^mean_log_Mh(log M*_lo)`.
- `M1_sat = 17 · M_min` (empirical ratio from Leauthaud+2011 Eq. 5).
- `α_sat = 1.0`.
- Satellites are suppressed below `M_min`.

To convert from a threshold sample to a bin sample, the satellite count
is scaled by the bin fraction: `N_cen(bin) / N_cen(threshold)`.

**The satellite parameters (α_sat, M1/Mmin ratio) are fixed — they are
NOT varied in the Fisher matrix.**  This is acceptable for a relative
survey comparison where the goal is to compare constraining power across
surveys, not to produce absolute error bars.


---

## 3. Observable Predictions

### 3.1 Galaxy-Galaxy Lensing: ΔΣ(R)

**Module:** `halo_model.py`, functions `delta_sigma_nfw()` and `delta_sigma_bin()`

The excess surface mass density around galaxies in a stellar mass bin is
the HOD-weighted average of single-halo NFW profiles:

```
ΔΣ(R | bin, z) = [∫ dMh (dn/dMh) · N_total(Mh) · ΔΣ_NFW(R, Mh, z)] / n_gal
```

where:
- `dn/dMh` is the Tinker+2008 halo mass function from colossus, converted
  from its native units `(dn/d ln M) [h³ Mpc⁻³]` to `[M☉⁻¹ Mpc⁻³]`.
- `N_total(Mh) = N_cen(Mh) + N_sat(Mh)`.
- `ΔΣ_NFW(R, Mh, z)` is the NFW excess surface mass density for a single
  halo, computed via colossus with the Diemer+2019 concentration-mass
  relation and the `200m` mass definition.
- `n_gal = ∫ dMh (dn/dMh) · N_total(Mh)` is the galaxy number density
  (the denominator normalizes to a per-galaxy signal).

The integral is evaluated via `np.trapz` over a log-spaced grid of
`n_Mh_bins` (default 200) halo masses from `log_Mh_min` to `log_Mh_max`
(default 10.0 to 15.5 in log₁₀ M☉).  The integration variable is
`d(log Mh)`, with the Jacobian `Mh · ln 10` included.

#### Unit conversion for colossus NFW profiles

Colossus profiles work internally in comoving kpc/h.  The conversion chain
is documented in `delta_sigma_nfw()`:

```
Input:   R [physical Mpc]
  → R_colossus = R × 10³ × h × (1+z)    [comoving kpc/h]

Output:  ΔΣ_colossus [h M☉ / comoving kpc²]
  → ΔΣ_physical = ΔΣ_colossus × h × (10³)² × (1+z)²    [M☉ / physical Mpc²]
```

The `h` factor removes the h in the numerator; `(10³)²` converts kpc² →
Mpc² in the denominator; `(1+z)²` converts comoving area → physical area.

**Effect of mass-dependent scatter on ΔΣ:**
The occupation function N_total(Mh) acts as a weighting kernel that selects
which halos contribute to a given stellar mass bin.  Broader scatter at low
masses mixes in halos from a wider mass range, which typically **dilutes**
the mean lensing signal (more low-mass halos with weak ΔΣ contribute) and
**shifts** the effective mean halo mass for each bin.


### 3.2 Galaxy Clustering: b_eff and n_gal

**Module:** `halo_model.py`, functions `effective_bias()` and `galaxy_number_density()`

Rather than predicting the full projected correlation function w_p(r_p),
the pipeline uses two summary statistics per (z, M*) bin:

**Effective bias:**
```
b_eff = [∫ dMh (dn/dMh) · N_total(Mh) · b(Mh, z)] / n_gal
```

where `b(Mh, z)` is the Tinker+2010 linear halo bias for the `200m`
definition, queried from colossus.

**Galaxy number density:**
```
n_gal = ∫ dMh (dn/dMh) · N_total(Mh)    [comoving Mpc⁻³]
```

Both integrals use the same log-spaced halo mass grid and `np.trapz` as
the lensing signal.

**Effect of mass-dependent scatter on clustering:**
Broader scatter changes the halo mass distribution contributing to each
stellar mass bin, which shifts b_eff.  For low-mass bins, the broader
scatter pulls in lower-mass, lower-bias halos, reducing b_eff.


---

## 4. Covariance Models

**Module:** `covariance.py`

### 4.1 Lensing Covariance (Shape Noise)

The lensing covariance is **diagonal** (no off-diagonal radial bin
correlations) and dominated by intrinsic galaxy shape noise:

```
Var(ΔΣ, Rᵢ) = σ_e² · Σ_crit_eff² / (N_lens · n_source_eff · A_annulus_i)
```

where:
- `σ_e = 0.26` — shape noise per ellipticity component (from `LensingConfig`).
- `Σ_crit_eff` — effective critical surface mass density, computed by
  averaging `Σ_crit⁻¹` over a Smail-type source redshift distribution
  `n(z) ∝ z² exp(−(z/z₀)^1.5)` with `z₀ = z_median / 1.41`
  (function `sigma_crit_effective()`).
- `N_lens` — number of spectroscopic lens galaxies in this (z, M*) bin
  (see §5.2 for how this is distributed).
- `n_source_eff` — effective source galaxy density [Mpc⁻²] behind the
  lens redshift.  Computed from `n_source_per_arcmin2` × fraction of
  sources behind `z_lens`, converted to physical Mpc⁻² using the angular
  diameter distance (function `n_source_effective()`).
- `A_annulus_i = π(R_hi² − R_lo²)` — projected annulus area at the lens
  [physical Mpc²].

#### Systematic floor (optional)

When `ForecastConfig.systematic_floor_fraction` > 0, an irreducible
systematic variance is added in quadrature:

```
Var_total = Var_shape + (f_sys · ΔΣ_fid)²
```

This represents systematic uncertainties from baryonic effects,
miscentering, and photometric calibration that cannot be reduced by
adding more lens-source pairs.  Typical values: `f_sys = 0.05–0.10`.


### 4.2 Clustering Covariance

**Effective bias variance** (cosmic variance + shot noise):
```
Var(b_eff) = b_eff² / (n_gal · V_eff)
```

where `V_eff = V_survey · [nP·b² / (1 + nP·b²)]²` is the FKP-weighted
effective volume, with the linear matter power spectrum P evaluated at
k = 0.1 h/Mpc (the scale carrying most bias information).

**Number density variance** (Poisson):
```
Var(n_gal) = n_gal / V_survey
```

The survey volume `V_survey` is computed in `survey_volume()` as the
fraction of the full spherical shell `(4/3)π(χ_hi³ − χ_lo³)` covered
by the survey area.


---

## 5. Survey Design and Its Effect on S/N

**Configuration:** `config.SurveyConfig`, predefined surveys in `survey_configs.py`

### 5.1 Survey Parameters

| Parameter | Description | Unit |
|-----------|-------------|------|
| `area_deg2` | Sky footprint | deg² |
| `z_min`, `z_max` | Redshift range | — |
| `n_gal_total` | Total spectroscopic galaxy count | — |
| `log_Mstar_min` | Stellar mass completeness floor | log₁₀(M☉) |
| `log_Mstar_min_func` | Optional z-dependent completeness | f(z) → log₁₀(M☉) |
| `dlog_Mstar` | Stellar mass bin width | dex (default 0.5) |
| `log_Mstar_max` | Upper edge of highest bin | log₁₀(M☉) (default 12.0) |

Stellar mass bins are generated by `get_stellar_mass_bins(z)`: edges run
from `get_log_Mstar_min(z)` to `log_Mstar_max` in steps of `dlog_Mstar`.

### 5.2 How N_lens Is Distributed

The total `n_gal_total` spectroscopic galaxies are distributed across
(z, M*) bins **proportional to the model-predicted galaxy count** in each
bin:

```
N_lens(z, M*) = n_gal_total × [n_gal_model(bin) · V(z-bin)] / Σ_all [n_gal_model · V]
```

This ensures that bins with more predicted galaxies get proportionally more
lens galaxies, matching the expected survey yield.

### 5.3 How Each Survey Parameter Affects Constraints

**Survey area (`area_deg2`):**
- Lensing: more lens-source pairs → lower shape noise variance
  (Var ∝ 1/N_lens, and N_lens ∝ area).
- Clustering: larger volume → more Fourier modes → lower cosmic variance
  (Var(b) ∝ 1/V_eff).
- Scaling: errors decrease roughly as area^(−1/2).

**Redshift range (`z_min`, `z_max`):**
- More redshift bins → more independent (z, M*) measurements → more
  Fisher information.
- Wider z-range → enables constraining redshift evolution parameters
  (ν_M1, ν_N, ν_β, ν_γ).  The pipeline auto-enables evolution
  parameters when `z_max > 0.4` and `z_max − z_min > 0.3`.
- Higher z_max → Σ_crit_eff decreases (more lensing efficiency from
  the larger lens-source separation), but source density behind the
  lens also decreases.  These partially cancel.

**Total galaxy count (`n_gal_total`):**
- Directly sets N_lens per bin → controls the lensing S/N.
- Does NOT affect clustering covariance, which depends on the
  model-predicted n_gal and the survey volume (not the observed count).

**Mass completeness (`log_Mstar_min`):**
- Determines how many stellar mass bins exist.
- Deeper completeness adds low-mass bins (M* < 10^9.5) that probe
  halos where σ(Mh) is largest (0.30–0.38 dex).
- These bins provide the strongest leverage on `scatter_sigma_rise`
  and the low-mass slope `beta_0`.

### 5.4 Predefined Surveys

| Name | Area | z-range | N_gal | log M*_min | Archetype |
|------|------|---------|-------|------------|-----------|
| `stage3_shallow_wide` | 7,500 deg² | 0.02–0.2 | 700k | 9.5 | SDSS/GAMA-like |
| `stage4_low_z` | 14,000 deg² | 0.05–0.4 | 10M | 9.0 | DESI BGS-like |
| `stage4_high_z` | 14,000 deg² | 0.4–1.0 | 5M | 10.8 | DESI LRG-like |
| `stage5_wide` | 10,000 deg² | 0.05–1.0 | 50M | 8.5 (+1.0·z) | MUST-like wide |
| `stage5_deep` | 3,000 deg² | 0.1–1.5 | 20M | 8.0 (+1.0·z) | MUST-like deep |

The imaging lensing survey is held fixed at Stage-IV levels across all
comparisons (`LensingConfig` defaults: n_source = 25 arcmin⁻², σ_e = 0.26,
z_source_median = 1.0).


---

## 6. Fisher Matrix Assembly

**Module:** `fisher.py`, main function `compute_fisher_matrix()`

### 6.1 Parameter Selection

The function `get_varied_params()` determines which SHMR parameters to
vary based on the scatter mode and survey redshift range:

**Always varied (z=0 shape):**
- `log_M1_0`, `N_0`, `beta_0`, `gamma_0`

**Scatter parameters (one of):**
- Constant mode: `sigma_logMs` (1 parameter)
- Mass-dependent mode: `scatter_sigma_high`, `scatter_sigma_rise` (2 parameters)

**Conditionally varied (if `vary_z_evolution=True` AND `z_max > 0.4`
AND `z_max − z_min > 0.3`):**
- `nu_M1`, `nu_N`, `nu_beta`, `nu_gamma`

**Optional nuisance parameters (if `include_nuisance_params=True`):**
- `shear_m` — multiplicative shear calibration bias (fiducial = 0)
- `photo_dz_source` — source photo-z bias (fiducial = 0)

Total parameter count ranges from 5 (constant scatter, low-z only) to
11 (mass-dependent scatter, wide-z, with nuisance).


### 6.2 Binning

1. **Redshift bins:** edges from `z_min` to `z_max` in steps of `dz`
   (default 0.2).  Each bin is evaluated at `z_mid = (z_lo + z_hi) / 2`.

2. **Stellar mass bins:** at each `z_mid`, `SurveyConfig.get_stellar_mass_bins(z_mid)`
   returns bins from the (possibly z-dependent) completeness limit to
   `log_Mstar_max = 12.0` in steps of `dlog_Mstar = 0.5 dex`.

3. **Radial bins for lensing:** `n_R_bins` (default 10) log-spaced bins
   from `R_min_Mpc = 0.1` to `R_max_Mpc = 30.0` physical Mpc.  Bin
   centers are geometric means: `R_center = √(R_lo · R_hi)`.


### 6.3 Numerical Derivatives

**Function:** `compute_derivatives()`

Derivatives of observables with respect to SHMR parameters are computed
via **central finite differences**:

```
∂f/∂θᵢ ≈ [f(θᵢ + δ) − f(θᵢ − δ)] / (2δ)
```

where `δ = frac_step · |θᵢ|` (default `frac_step = 0.01`, i.e., 1%).
For parameters near zero (|θᵢ| < 10⁻¹⁰, e.g., `nu_N`), an absolute
step of 0.01 is used instead.

Each derivative evaluation requires two full observable computations
(one at θ + δ, one at θ − δ), each of which involves the HMF integral
over the halo mass grid.  For P parameters and B bins, this costs
`2 × P × B` halo model evaluations.

The result is a **Jacobian matrix** `J` of shape `(N_params, N_obs)`.

For **nuisance parameters**, the derivatives are analytic:
- `∂(ΔΣ_obs)/∂m = ΔΣ_fid` (shear calibration: ΔΣ_obs = (1+m) · ΔΣ_true)
- `∂(ΔΣ_obs)/∂(δz_s) = ΔΣ_fid · d(ln Σ_crit_eff⁻¹)/d(δz_s)` (photo-z bias)

Nuisance parameters do not affect clustering observables, so their
clustering Jacobian rows are zero.


### 6.4 Fisher Accumulation

For each (z-bin, M*-bin), the Fisher matrix receives three additive
contributions:

**Lensing (ΔΣ):**
```
F += J_ΔΣᵀ · C_ΔΣ⁻¹ · J_ΔΣ
```

where `C_ΔΣ` is diagonal (shape noise + optional systematic floor),
`J_ΔΣ` has shape `(N_params, N_R_bins)`, and only radial bins with
valid (positive, finite) variance and signal are included.

**Effective bias (b_eff):**
```
F += J_bᵀ · J_b / Var(b_eff)
```

where `J_b` has shape `(N_params, 1)`.

**Galaxy number density (n_gal):**
```
F += J_nᵀ · J_n / Var(n_gal)
```

where `J_n` has shape `(N_params, 1)`.

The total Fisher matrix is the sum over all (z, M*) bins.

### 6.5 Priors on Nuisance Parameters

When nuisance parameters are included, Gaussian priors are added to the
Fisher diagonal after accumulation:

```
F[i_shear, i_shear]   += 1 / σ_m²
F[i_photoz, i_photoz]  += 1 / σ_δz²
```

Default prior widths from `NuisanceConfig`: σ_m = 0.02, σ_δz = 0.03.


### 6.6 Post-Processing

**Marginalized errors** (the primary output):
```
σ(θᵢ) = √[(F⁻¹)ᵢᵢ]
```

These account for all parameter degeneracies.  When nuisance parameters
are included, their contribution is automatically marginalized out in
the matrix inversion.

**Conditional (unmarginalized) errors:**
```
σ_cond(θᵢ) = 1 / √(Fᵢᵢ)
```

Always smaller than marginalized errors.

**Fisher ellipses** for parameter pairs `(i, j)`:
Eigendecomposition of the 2×2 sub-covariance from `(F⁻¹)`, plotted as
rotated ellipses at n-σ contour levels.

**Condition number check:**
If `cond(F) > 10¹⁰`, a warning is issued.  Common causes: evolution
parameters unconstrained by a low-z-only survey, or strong parameter
degeneracies (e.g., `log_M1_0` and `N_0`).


---

## 7. Pipeline Data Flow Summary

```
SHMRParams (config.py)
  │
  ├─ mean_log_Mstar(log_Mh, z)          # shmr_model.py
  ├─ scatter_at_Mh(log_Mh)              # shmr_model.py
  │
  ▼
HOD: N_cen(Mh), N_sat(Mh)               # halo_model.py
  │
  ├─── × dn/dMh (Tinker08 HMF)
  ├─── × ΔΣ_NFW(R, Mh, z)  ─────────── ΔΣ(R | bin, z)     # halo_model.py
  ├─── × b(Mh, z) (Tinker10 bias) ───── b_eff(bin, z)      # halo_model.py
  └─── (integrate) ──────────────────── n_gal(bin, z)       # halo_model.py
                                              │
                                              ▼
                                    Covariance (covariance.py)
                                    ├── Var(ΔΣ):  shape noise + sys. floor
                                    ├── Var(b):   cosmic var. + shot noise
                                    └── Var(n):   Poisson
                                              │
                                              ▼
                                    Fisher matrix (fisher.py)
                                    ├── J = ∂(observables)/∂(params)
                                    ├── F = Σ_bins  Jᵀ C⁻¹ J
                                    ├── + nuisance priors
                                    └── σ(θᵢ) = √[(F⁻¹)ᵢᵢ]
```


---

## 8. Lensing Profile: 1-Halo NFW + 2-Halo Term

### 8.1 Current Implementation

The lensing prediction combines two physical components:

```
ΔΣ_total(R) = ΔΣ_1h(R) + ΔΣ_2h(R)
```

**1-halo term (`delta_sigma_nfw`):**
The inner halo profile uses a standard NFW model via colossus, with the
Diemer+2019 concentration-mass relation and the `200m` mass definition.
Colossus provides an analytic Wright & Brainerd (2000) closed-form
expression for ΔΣ_NFW, making this evaluation fast (~5 ms per halo).

**2-halo term (`delta_sigma_2halo`):**
At large projected radii (R ≳ 2–5 Mpc), the correlated large-scale
structure around halos contributes additional lensing signal.  This is
modeled as:

```
ΔΣ_2h(R, Mh, z) = b_h(Mh, z) · ΔΣ_matter(R, z)
```

where:
- `b_h(Mh, z)` is the Tinker+2010 linear halo bias (same as used for
  the clustering observable).
- `ΔΣ_matter(R, z)` is the excess surface mass density of the mean
  matter field, computed from the matter correlation function ξ_mm(r)
  by numerical line-of-sight projection and azimuthal averaging.

The matter correlation function is obtained from colossus
`correlation.xi(r, z)`, which evaluates the nonlinear matter power
spectrum (Eisenstein & Hu 1998 transfer function + Smith+2003 halofit
corrections) and Fourier-transforms to configuration space.

**Projection integral for ΔΣ_matter:**

```
Σ_matter(R, z) = ρ_m(z) · ∫_{-π_max}^{+π_max} [1 + ξ_mm(√(R² + π²), z)] dπ
```

The line-of-sight integral is truncated at `pi_max` (default 100 Mpc)
and computed via `scipy.integrate.quad`.  The excess surface density is
then:

```
ΔΣ_matter(R) = Σ̄_matter(<R) − Σ_matter(R)
```

where `Σ̄(<R)` is the mean surface density within radius R, computed by
integrating `Σ_matter(R') · 2πR' dR'` from 0 to R and dividing by πR².

The 2-halo ΔΣ is precomputed once per redshift on the radial grid
(`_delta_sigma_matter_cached`) and cached, since it does not depend on
halo mass.  The halo-mass-dependent bias factor is applied at integration
time in `delta_sigma_bin()`.

**Bin-averaged signal including 2-halo:**

```
ΔΣ(R | bin, z) = [∫ dMh (dn/dMh) · N_total(Mh) · [ΔΣ_NFW(R,Mh) + b_h(Mh) · ΔΣ_mm(R)]] / n_gal
```

This is equivalent to:

```
ΔΣ_total = ΔΣ_1h + b_eff · ΔΣ_mm
```

where `b_eff` is the HOD-weighted effective bias — the same quantity
used for the clustering observable.  However, the integral form is used
to correctly weight by halo mass.


### 8.2 Investigation of Colossus `diemer23` Profile

Before implementing the analytic 2-halo approach, we investigated using
the Diemer (2023) profile from colossus as an alternative to NFW.  This
profile models the splashback-truncated orbiting component of halos and
can be combined with infalling and 2-halo outer terms via
`profile_composite.compositeProfile()`.

#### Available outer terms in colossus

| Shortcode | Class | Description |
|-----------|-------|-------------|
| `'infalling'` | `OuterTermInfalling` | Power-law infalling matter with smooth transition. Parameters: `pl_delta_1` (normalization), `pl_s` (slope). |
| `'cf'` | `OuterTermCorrelationFunction` | True 2-halo term: `ρ_m · b_h · ξ_mm(r)`. Derives bias self-consistently from R200m. |
| `'mean'` | `OuterTermMeanDensity` | Constant mean matter density. Causes `surfaceDensity` to diverge; not suitable for lensing. |

#### Benchmark results (M = 10¹³ M☉, z = 0.3, 20 radial points)

| Profile variant | Time per halo | Notes |
|-----------------|--------------|-------|
| NFW (analytic W&B) | 0.005 s | Closed-form ΔΣ |
| D23 orbiting only | 0.22 s | Numerical line-of-sight integration |
| D23 + infalling | 0.38 s | Additional outer term integration |

**Slowdown factor: ~77×** compared to NFW.

For the full Fisher pipeline (200 halo masses × 255 derivative
evaluations), this translates to:

| Profile | Estimated wall time |
|---------|-------------------|
| NFW | ~4 minutes |
| D23 + infalling | ~5.5 hours |

#### Profile comparison (ratio to NFW at z = 0.3)

| Region | D23 orbiting / NFW | D23 + infalling / NFW |
|--------|-------------------|----------------------|
| R ~ 0.1 Mpc (inner) | ~1.03 | ~0.94 |
| R ~ 1 Mpc (1h–2h transition) | ~0.98 | ~1.00 |
| R ~ 5 Mpc (outer) | ~0.50 | ~1.27 |
| R ~ 20 Mpc (far outer) | ~0.34 | ~5.73 |

The D23 orbiting-only profile truncates faster than NFW at large R
(by design — it models the splashback truncation).  The infalling term
overcompensates at R > 5 Mpc.  At R < 2 Mpc, where most of the Fisher
information resides, the profiles agree to within ~3%.

#### Why the `cf` (2-halo) term fails

The correlation function outer term (`OuterTermCorrelationFunction`)
attempts to add `ρ_m · b_h · ξ_mm(r)` to the 3D density profile.
Because ξ_mm goes negative at large separations, the total density
becomes negative, which causes:

1. With `interpolate=True` (default): "Found negative value in density,
   cannot create interpolation table" error.
2. With `interpolate=False`: numerical integration failures in
   `surfaceDensityInner()` when `quad` encounters zero-size arrays from
   filtering out negative density values.

This is a known limitation of the colossus composite profile framework
for lensing projections.

#### Decision rationale

Given the 77× speed penalty and the fact that the inner profile (where
most Fisher information lives) differs by only ~3%, we chose the
**NFW + analytic 2-halo** approach:

- Keep the fast analytic NFW for the 1-halo term.
- Add `ΔΣ_2h = b_h · ΔΣ_matter` using colossus's matter correlation
  function, with line-of-sight projection computed once per redshift.
- This gives physically correct large-R behavior with negligible
  speed overhead (ΔΣ_matter is cached per redshift).

The comparison plot is saved as
`outputs/profile_comparison_nfw_vs_diemer23.png`.


---

## 9. Key Assumptions and Limitations

1. **NFW + linear 2-halo term.**  The lensing signal uses an NFW profile
   for the 1-halo term plus a linear-bias 2-halo term.  The transition
   region (R ~ 1–3 Mpc) is not explicitly modeled — the NFW and 2-halo
   terms are simply summed.  This double-counts some matter near the
   virial radius but has negligible effect on the Fisher matrix, which
   is driven by derivatives rather than absolute signal levels.

2. **Satellite profile = host NFW.**  Satellites are assumed to trace the
   host halo's NFW profile.  No subhalo lensing or off-centering is
   modeled.

3. **Diagonal lensing covariance.**  No off-diagonal correlations between
   radial bins from large-scale structure or intrinsic alignments.  This
   underestimates the true covariance and leads to slightly optimistic
   constraints.

4. **Clustering reduced to two numbers.**  The full w_p(r_p) information
   content is compressed to (b_eff, n_gal) per bin.  This discards
   scale-dependent clustering information (e.g., the transition from
   1-halo to 2-halo) but simplifies the covariance model.

5. **Fixed satellite prescription.**  The satellite parameters (α_sat,
   M1/Mmin) are not varied.  This is acceptable for relative survey
   comparisons but means the absolute error bars should not be taken at
   face value.

6. **No scatter redshift evolution.**  σ(Mh) is assumed constant with
   redshift.  This is broadly consistent with current empirical
   constraints but is an area for future improvement.

7. **No covariance between stellar mass bins.**  Each (z, M*) bin
   contributes independently to the Fisher matrix.  In reality,
   overlapping halo populations create correlations between adjacent
   M* bins.


---

## References

- Moster, Naab & White (2013), MNRAS, 428, 3121 — SHMR parameterization
- Cao & Tinker (2020), MNRAS, 498, 5080 — mass-dependent SHMR scatter
- Leauthaud et al. (2012), ApJ, 744, 159 — HOD from SHMR
- Tinker et al. (2008), ApJ, 688, 709 — halo mass function
- Tinker et al. (2010), ApJ, 724, 878 — halo bias
- Diemer & Kravtsov (2015); Diemer (2019) — concentration-mass relation
- Diemer (2023), ApJ, 947, 59 — splashback-truncated density profiles
- Wright & Brainerd (2000), ApJ, 534, 34 — analytic NFW lensing signal
- Eisenstein & Hu (1998), ApJ, 496, 605 — matter transfer function
- Smith et al. (2003), MNRAS, 341, 1311 — halofit nonlinear power spectrum
- Tegmark, Taylor & Heavens (1997), ApJ, 480, 22 — Fisher matrix formalism
- Oguri & Takada (2011), PRD, 83, 023008 — Fisher forecast for lensing
- Behroozi et al. (2019), MNRAS, 488, 3143 — SHMR scatter evolution
