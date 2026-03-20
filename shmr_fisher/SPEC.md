# SPEC.md — SHMR Forecast Tool: Technical Specification

## 1. Scope

A Python package that computes Fisher matrix forecasts for SHMR parameter constraints from joint galaxy-galaxy lensing + galaxy clustering, given arbitrary spectroscopic survey configurations and fixed Stage-IV lensing assumptions. The tool is first-order and analytic (no simulations, no MCMC), designed for rapid survey comparison and parameter sweeps.

---

## 2. Module Specifications

### 2.1 `config.py` — Configuration Dataclasses

#### `SHMRParams`
Stores the Moster+2013 SHMR model with optional mass-dependent scatter (Cao & Tinker 2020).

**Core parameters (9 — always present):**

| Parameter | Symbol | Fiducial | Physical meaning |
|-----------|--------|----------|------------------|
| `log_M1_0` | log₁₀(M₁/M☉) at z=0 | 11.59 | Characteristic halo mass at z=0 |
| `N_0` | N(z=0) | 0.0351 | Peak stellar-to-halo mass ratio at z=0 |
| `beta_0` | β(z=0) | 1.376 | Low-mass power-law slope at z=0 |
| `gamma_0` | γ(z=0) | 0.608 | High-mass power-law slope at z=0 |
| `nu_M1` | νM₁ | 1.195 | Redshift evolution of log M₁ |
| `nu_N` | νN | −0.0247 | Redshift evolution of N |
| `nu_beta` | νβ | −0.826 | Redshift evolution of β |
| `nu_gamma` | νγ | 0.329 | Redshift evolution of γ |
| `sigma_logMs` | σ_log M* | 0.15 | Log-normal scatter at fixed Mh [dex], constant with z (used when `use_mass_dependent_scatter=False`) |

**Mass-dependent scatter parameters (4 — active when `use_mass_dependent_scatter=True`):**

Following Cao & Tinker (2020, MNRAS, 498, 5080), the scatter varies with halo mass via a tanh transition:

    σ(log M* | Mh) = σ_high + σ_rise × [1 − tanh((log Mh − log Mh,break) / δ)]

| Parameter | Symbol | Fiducial | Physical meaning |
|-----------|--------|----------|------------------|
| `use_mass_dependent_scatter` | — | False | Toggle mass-dependent scatter model |
| `scatter_sigma_high` | σ_high | 0.18 | Asymptotic scatter at high Mh [dex] |
| `scatter_sigma_rise` | σ_rise | 0.10 | Half-amplitude of rise toward low Mh [dex] |
| `scatter_log_Mh_break` | log Mh,break | 12.0 | Transition halo mass [log₁₀(M☉)] |
| `scatter_delta` | δ | 0.4 | Transition width [dex] |

When mass-dependent scatter is enabled, `sigma_logMs` is replaced by these 4 parameters. The Fisher matrix varies `scatter_sigma_high` and `scatter_sigma_rise` (the other two are held fixed). At high Mh (≫ break), σ → σ_high ≈ 0.18 dex. At low Mh (≪ break), σ → σ_high + 2σ_rise ≈ 0.38 dex.

Must support:
- `to_dict()` → dictionary of {name: value} for all parameters
- `from_dict(d)` → construct from dictionary
- `copy(**overrides)` → return a new instance with specified parameters replaced (for finite differences)

#### `SurveyConfig`
| Field | Type | Units | Description |
|-------|------|-------|-------------|
| `name` | str | — | Human-readable label |
| `area_deg2` | float | deg² | Sky area |
| `z_min` | float | — | Minimum redshift |
| `z_max` | float | — | Maximum redshift |
| `n_gal_total` | float | — | Total number of spec-z galaxies |
| `log_Mstar_min` | float | log₁₀(M☉) | Global mass completeness floor |
| `log_Mstar_max` | float | log₁₀(M☉) | Upper bin edge (default 12.0) |
| `dlog_Mstar` | float | dex | Stellar mass bin width (default 0.5) |
| `log_Mstar_min_func` | callable or None | — | Optional z-dependent completeness: f(z) → log₁₀(M☉) |

Methods:
- `get_log_Mstar_min(z)` → effective mass limit at redshift z
- `get_stellar_mass_bins(z)` → list of (lo, hi) tuples
- `get_n_gal_per_deg2()` → average surface density

#### `LensingConfig`
Fixed for all forecasts. Fields: `n_source_per_arcmin2` (25.0), `sigma_e` (0.26), `z_source_median` (1.0).

#### `NuisanceConfig`
Lensing nuisance parameter priors for marginalization. Fields: `sigma_m` (0.02, shear calibration prior), `sigma_dz_source` (0.03, source photo-z bias prior).

#### `ForecastConfig`
Controls forecast behavior. Fields: `R_min_Mpc`, `R_max_Mpc`, `n_R_bins`, `dz`, `frac_step`, `vary_z_evolution`, `log_Mh_min`, `log_Mh_max`, `n_Mh_bins`, `systematic_floor_fraction` (0.0), `include_nuisance_params` (False).

---

### 2.2 `shmr_model.py` — SHMR Parameterization

#### Core functions

```python
def shmr_params_at_z(params: SHMRParams, z: float) -> tuple[float, float, float, float]:
    """
    Evaluate z-dependent SHMR parameters.
    
    Moster+2013 Eq. 11–14:
        log₁₀(M₁(z))  = log_M1_0  + nu_M1  × z/(1+z)
        N(z)           = N_0       + nu_N   × z/(1+z)
        β(z)           = beta_0    + nu_beta × z/(1+z)
        γ(z)           = gamma_0   + nu_gamma × z/(1+z)

    Parameters
    ----------
    params : SHMRParams
    z : float
        Redshift.

    Returns
    -------
    log_M1, N, beta, gamma : float
        SHMR parameters evaluated at redshift z.
    """
```

```python
def mean_log_Mstar(log_Mh: np.ndarray | float, params: SHMRParams, z: float) -> np.ndarray | float:
    """
    Mean log₁₀(M*/M☉) as a function of halo mass and redshift.

    Moster+2013 Eq. 2:
        M*/Mh = 2N × [(Mh/M₁)^(−β) + (Mh/M₁)^γ]^(−1)

    Parameters
    ----------
    log_Mh : array-like
        log₁₀(Mh/M☉). Must be > 0.
    params : SHMRParams
    z : float

    Returns
    -------
    log_Mstar : same shape as log_Mh
        log₁₀(M*/M☉).
    """
```

```python
def mean_log_Mh(log_Mstar: float, params: SHMRParams, z: float) -> float:
    """
    Inverse SHMR via root-finding.

    Uses scipy.optimize.brentq on [log_Mh_min=9, log_Mh_max=16].
    Raises ValueError if no root found.
    """
```

```python
def phi_Mstar_given_Mh(log_Mstar: np.ndarray, log_Mh: float,
                        params: SHMRParams, z: float) -> np.ndarray:
    """
    P(log M* | Mh, z): log-normal PDF.

    P = (1 / (√(2π) σ)) × exp[−(log M* − ⟨log M*⟩)² / (2σ²)]

    Returns probability density in units of per dex.
    """
```

#### Validation criteria (Task 3)
- At z=0: peak of M*/Mh should be between 0.02 and 0.06, occurring at log Mh between 11.5 and 12.5.
- `mean_log_Mstar(12.0, fiducial, z=0)` should return ~10.4–10.8.
- `mean_log_Mh(mean_log_Mstar(12.0, ...))` should return 12.0 (round-trip test).
- Plot: M*/Mh vs Mh at z=0, 0.5, 1.0, 2.0 → save as `outputs/phase1/shmr_validation.png`.

---

### 2.3 `halo_model.py` — Observable Predictions

#### HOD functions

```python
def N_cen(log_Mh_grid: np.ndarray, log_Mstar_lo: float, log_Mstar_hi: float,
          params: SHMRParams, z: float) -> np.ndarray:
    """
    Mean central occupation for a stellar mass bin.

    N_cen(Mh) = 0.5 × [erf(x_hi) − erf(x_lo)]
    where x = (log M* − ⟨log M*(Mh,z)⟩) / (√2 × σ_logMs)

    Parameters
    ----------
    log_Mh_grid : array, shape (N_Mh,)
        Halo mass grid in log₁₀(M☉).
    log_Mstar_lo, log_Mstar_hi : float
        Stellar mass bin edges in log₁₀(M☉).

    Returns
    -------
    n_cen : array, shape (N_Mh,)
        Mean number of centrals per halo.
    """
```

```python
def N_sat(log_Mh_grid: np.ndarray, log_Mstar_lo: float, log_Mstar_hi: float,
          params: SHMRParams, z: float) -> np.ndarray:
    """
    Mean satellite occupation (fixed prescription, not varied in Fisher).

    Uses Leauthaud+2011 form with fixed alpha_sat=1.0, M1_sat/M_min=17.
    """
```

#### Lensing signal

```python
def delta_sigma_nfw(R_Mpc: np.ndarray, Mh: float, z: float) -> np.ndarray:
    """
    ΔΣ(R) for a single NFW halo using colossus.

    Parameters
    ----------
    R_Mpc : array
        Projected radii in physical Mpc.
    Mh : float
        Halo mass M200m in M☉ (NOT M☉/h).
    z : float

    Returns
    -------
    delta_sigma : array
        ΔΣ in M☉/Mpc². Same shape as R_Mpc.

    Notes
    -----
    Unit conversion from colossus (comoving kpc/h) is handled internally.
    See CLAUDE.md for fallback strategy if units are problematic.
    """
```

```python
def delta_sigma_bin(R_Mpc: np.ndarray, log_Mstar_lo: float, log_Mstar_hi: float,
                    params: SHMRParams, z: float,
                    log_Mh_grid: np.ndarray = None) -> np.ndarray:
    """
    Halo-model prediction for ΔΣ(R) averaged over a stellar mass bin.

    Total signal = 1-halo + 2-halo:
        ΔΣ_1h(R|bin,z) = ∫ dMh (dn/dMh) [N_cen × ΔΣ_NFW + N_sat × ΔΣ_sat] / n_gal
        ΔΣ_2h(R|bin,z) = b_eff(bin,z) × ΔΣ_matter(R, z)

    where ΔΣ_matter is the matter-matter excess surface density computed from
    the projected matter correlation function ξ_mm (via colossus), cached per
    redshift for performance.

    Integration is done via np.trapezoid over log_Mh_grid.

    Parameters
    ----------
    R_Mpc : array, shape (N_R,)
    log_Mstar_lo, log_Mstar_hi : float
    params : SHMRParams
    z : float
    log_Mh_grid : array, optional
        Halo mass grid. Default: np.linspace(10, 15.5, 200).

    Returns
    -------
    ds : array, shape (N_R,)
        ΔΣ in M☉/Mpc² (1-halo + 2-halo).
    """
```

#### Clustering summary statistics

```python
def effective_bias(log_Mstar_lo: float, log_Mstar_hi: float,
                   params: SHMRParams, z: float,
                   log_Mh_grid: np.ndarray = None) -> float:
    """
    HOD-weighted effective halo bias for a stellar mass bin.

    b_eff = ∫ dMh (dn/dMh) [N_cen + N_sat] × b(Mh,z) / n_gal

    Uses colossus.lss.bias.haloBias(Mh, z, mdef='200m', model='tinker10').
    """
```

```python
def galaxy_number_density(log_Mstar_lo: float, log_Mstar_hi: float,
                          params: SHMRParams, z: float,
                          log_Mh_grid: np.ndarray = None) -> float:
    """
    Predicted comoving galaxy number density [Mpc^{-3}] in a stellar mass bin.

    n_gal = ∫ dMh (dn/dMh) [N_cen + N_sat]
    """
```

#### Validation criteria (Tasks 4–6)
- `delta_sigma_nfw(R, 1e13, 0.3)` at R=0.5 Mpc should give ΔΣ ~ 10–100 M☉/pc² (i.e., 10^13–10^14 M☉/Mpc²).
- `effective_bias` for a bin around M* ~ 10^{10.5} at z=0.3 should be ~1.0–2.0.
- `galaxy_number_density` for M* ∈ [10.0, 10.5] at z=0.3 should be ~10^{-3}–10^{-2} Mpc^{-3}.
- Plot: ΔΣ(R) for bins [9.5,10], [10,10.5], [10.5,11], [11,11.5] at z=0.3 → save as `outputs/phase2/delta_sigma_model.png`.

---

### 2.4 `covariance.py` — Noise Model

```python
def sigma_crit_effective(z_lens: float, lensing_config: LensingConfig) -> float:
    """
    Effective critical surface mass density for a lens at z_lens,
    averaged over the source redshift distribution.

    ⟨Σ_crit^{-1}⟩ = ∫ dz_s  n(z_s)  Σ_crit^{-1}(z_l, z_s)   [for z_s > z_l]

    Source n(z_s): Smail distribution n(z) ∝ z² exp(−(z/z₀)^{1.5})
    with z₀ = z_source_median / 1.41.

    Returns Σ_crit,eff = ⟨Σ_crit^{-1}⟩^{-1} in M☉/Mpc².
    """
```

```python
def n_source_effective(z_lens: float, lensing_config: LensingConfig) -> float:
    """
    Effective source density [per arcmin²] behind the lens redshift.

    n_eff(z_l) = n_total × ∫_{z_l}^∞ n(z_s) dz_s / ∫_0^∞ n(z_s) dz_s
    """
```

```python
def lensing_covariance(R_bins_Mpc: np.ndarray, z_lens: float, N_lens: float,
                       survey_area_deg2: float,
                       lensing_config: LensingConfig) -> np.ndarray:
    """
    Diagonal covariance for ΔΣ in radial bins (shape noise only).

    σ²(ΔΣ, R_i) = σ_e² × Σ_crit_eff² / (N_lens × n_source_eff × A_annulus_i)

    where A_annulus_i = π(R_{i+1}² − R_i²) in physical Mpc² per lens,
    converted to sky area using the angular diameter distance.

    Parameters
    ----------
    R_bins_Mpc : array, shape (N_bins+1,)
        Radial bin edges in physical Mpc.
    z_lens : float
    N_lens : float
        Number of lens galaxies in this (z, M*) bin.
    survey_area_deg2 : float
    lensing_config : LensingConfig

    Returns
    -------
    var_ds : array, shape (N_bins,)
        Variance σ²(ΔΣ) at each radial bin center, in (M☉/Mpc²)².
    """
```

```python
def clustering_covariance(N_gal: float, b_eff: float,
                          survey_volume_Mpc3: float, z: float) -> tuple[float, float]:
    """
    Variance on (b_eff, n_gal) from galaxy clustering.

    σ²(b_eff) = b_eff² / (n_gal × V_eff)   [cosmic variance limit]
    σ²(n_gal) = n_gal / V_survey             [Poisson]

    For the bias, we use V_eff ≈ V_survey × [n_gal × P_lin(k=0.1, z) × b_eff²]
                                             / [1 + n_gal × P_lin × b_eff²]
    where P_lin is evaluated at k ~ 0.1 h/Mpc (the scale carrying most of the bias info).

    Parameters
    ----------
    N_gal : float
        Total galaxies in this (z, M*) bin.
    b_eff : float
        Effective halo bias.
    survey_volume_Mpc3 : float
        Comoving survey volume in this z-bin [Mpc³].
    z : float

    Returns
    -------
    var_b : float
        Variance on b_eff.
    var_n : float
        Variance on n_gal [Mpc^{-6}].
    """
```

```python
def survey_volume(z_lo: float, z_hi: float, area_deg2: float) -> float:
    """
    Comoving survey volume in a redshift shell.
    Uses colossus cosmology for comoving volume element.

    Returns volume in Mpc³.
    """
```

#### Validation criteria (Tasks 7–8)
- For N_lens = 10^5, n_s = 25/arcmin², z_l = 0.3: σ(ΔΣ) at R ~ 1 Mpc should be ~1–10 M☉/pc².
- σ(ΔΣ) should scale as 1/√N_lens. Test by comparing N_lens = 10^4 vs 10^6.
- `survey_volume(0.2, 0.4, 10000)` should be ~10^9 Mpc³ (order of magnitude).

#### 2.4.1 Systematic Error Model (Optional Features)

Two toggleable systematic error features are available to produce more realistic absolute constraints. The statistical-only forecast is sufficient for relative survey comparisons; these features are needed when absolute error bars matter.

**Feature 1: Fractional systematic floor on ΔΣ**

Adds an irreducible fractional uncertainty to each ΔΣ radial bin:

    σ²_total(R_i) = σ²_stat(R_i) + (f_sys × ΔΣ_fid(R_i))²

This captures baryonic effects on halo profiles, miscentering, intrinsic alignments, and other model uncertainties that scale with the signal amplitude. Toggle: `ForecastConfig.systematic_floor_fraction` (default 0.0 = off; typical 0.05–0.10).

**Feature 2: Nuisance parameter marginalization**

Two lensing nuisance parameters are added to the Fisher matrix:

| Parameter | Symbol | Fiducial | Prior σ | Physical effect |
|-----------|--------|----------|---------|-----------------|
| `shear_m` | m | 0 | 0.02 | Multiplicative shear calibration bias: ΔΣ_obs = (1+m) ΔΣ_true |
| `photo_dz_source` | Δz_s | 0 | 0.03 | Source photo-z bias: shifts n(z_s), modifying Σ_crit_eff |

Derivatives are computed analytically (not via finite differences on SHMRParams):
- ∂ΔΣ/∂m = ΔΣ_fid (exact at fiducial m=0)
- ∂ΔΣ/∂Δz_s = ΔΣ_fid × ∂ln(Σ_crit_eff⁻¹)/∂Δz_s (numerical via `d_ln_inv_sigma_crit_d_dz`)

Clustering observables (b_eff, n_gal) are unaffected by lensing nuisance parameters.

Toggle: `ForecastConfig.include_nuisance_params` (default False). Configuration: `NuisanceConfig` dataclass with `sigma_m` and `sigma_dz_source`.

**Feature 3: Off-diagonal LSS covariance — DEFERRED**

Not implemented. The off-diagonal covariance from large-scale structure requires Hankel transforms (J2 Bessel functions) of the ΔΣ window function against P(k), plus survey window treatment (~200–300 lines). Deferred because: (1) shape noise dominates for stacked g-g lensing, (2) the systematic floor already captures leading model uncertainties, (3) off-diagonal terms matter mainly at large R where S/N is low. A future parametric approximation (correlation coefficient r ~ 0.1–0.3 between adjacent bins) would be simpler.

---

### 2.5 `fisher.py` — Fisher Matrix

```python
def get_varied_params(params: SHMRParams, forecast_config: ForecastConfig,
                      z_min: float, z_max: float) -> list[tuple[str, float]]:
    """
    Auto-select parameters to vary based on survey z-range.

    Rules:
    - Always vary: log_M1_0, N_0, beta_0, gamma_0, sigma_logMs (5 params)
    - Also vary nu_M1, nu_N, nu_beta, nu_gamma IF:
        forecast_config.vary_z_evolution is True
        AND z_max > 0.4
        AND (z_max - z_min) > 0.3

    Returns list of (param_name, fiducial_value).
    """
```

```python
def compute_derivatives(observable_func: callable, params: SHMRParams,
                        varied_params: list[tuple[str, float]],
                        frac_step: float, **kwargs) -> np.ndarray:
    """
    Numerical derivatives of an observable w.r.t. SHMR parameters.

    Uses central finite differences:
        ∂f/∂θ_i ≈ [f(θ_i + δθ) − f(θ_i − δθ)] / (2δθ)
    where δθ = frac_step × |θ_i|.

    Special case: if |θ_i| < 1e-10 (parameter near zero, e.g. nu_N),
    use absolute step of 0.01 instead.

    Parameters
    ----------
    observable_func : callable
        f(params, **kwargs) → array of shape (N_obs,).
    params : SHMRParams
        Fiducial parameters.
    varied_params : list of (name, fiducial_value)
    frac_step : float

    Returns
    -------
    derivs : array, shape (N_params, N_obs)
        Jacobian matrix.
    """
```

```python
def compute_fisher_matrix(
    shmr_params: SHMRParams,
    survey_config: SurveyConfig,
    lensing_config: LensingConfig,
    forecast_config: ForecastConfig,
    nuisance_config: NuisanceConfig | None = None,
) -> tuple[np.ndarray, list[str], dict]:
    """
    Main Fisher matrix computation.

    Algorithm:
    1. Build redshift bins: edges from z_min to z_max in steps of dz.
    2. For each z-bin (evaluated at z_mid):
       a. Get stellar mass bins via survey_config.get_stellar_mass_bins(z_mid).
       b. For each M*-bin:
          - Compute N_lens in this bin (from survey's total N_gal,
            distributed across bins proportional to predicted n_gal × volume).
          - Compute fiducial observables: ΔΣ(R), b_eff, n_gal.
          - Compute covariance: lensing noise, clustering noise.
          - Compute Jacobian ∂obs/∂θ.
          - Accumulate Fisher contribution:
            F += J^T × C^{-1} × J
    3. Check condition number; warn if > 1e10.
    4. Return (fisher_matrix, param_names, metadata).

    The metadata dict should contain:
        - z_bins: list of (z_lo, z_hi)
        - Mstar_bins_per_z: dict mapping z_mid → list of (M*_lo, M*_hi)
        - N_lens_per_bin: dict mapping (z_mid, M*_lo) → N_lens
        - fiducial_observables: stored for plotting
        - condition_number: float
    """
```

```python
def marginalized_errors(fisher_matrix: np.ndarray) -> np.ndarray:
    """σ(θ_i) = √(F^{-1}_{ii})."""
```

```python
def conditional_errors(fisher_matrix: np.ndarray) -> np.ndarray:
    """σ_cond(θ_i) = 1/√(F_{ii}). Unmarginalized, for comparison."""
```

```python
def fisher_ellipse(fisher_matrix: np.ndarray, i: int, j: int,
                   n_sigma: float = 1.0, n_points: int = 100) -> tuple[np.ndarray, np.ndarray]:
    """
    Return (x, y) coordinates of the n_sigma Fisher ellipse
    for parameter pair (i, j), marginalized over all other parameters.
    """
```

```python
def extract_shmr_constraints(fisher: np.ndarray, param_names: list[str]
                              ) -> tuple[np.ndarray, list[str]]:
    """
    Extract SHMR-only constraints by marginalizing over nuisance parameters.

    If the Fisher matrix includes nuisance parameters (shear_m, photo_dz_source),
    this function marginalizes over them and returns errors on only the
    SHMR parameters.
    """
```

```python
def add_external_prior(fisher: np.ndarray, param_names: list[str],
                       param_name: str, sigma_prior: float) -> np.ndarray:
    """Add Gaussian prior: F[i,i] += 1/σ²."""
```

#### Validation criteria (Task 9)
- Fisher matrix must be symmetric and positive definite.
- Condition number < 10^12 for 5-param model; < 10^14 for 9-param model (otherwise warn).
- Doubling N_gal_total should reduce marginalized errors by ~√2 (check within 20%).
- `marginalized_errors` ≥ `conditional_errors` element-wise (always true by construction).

---

### 2.6 `survey_configs.py` — Predefined Surveys

Provide 5 predefined configurations (Stage-III through Stage-V) as documented in the plan. Also provide:

```python
def make_custom_survey(name, area_deg2, z_min, z_max, n_gal_total,
                       log_Mstar_min, dlog_Mstar=0.5, **kwargs) -> SurveyConfig:
    """Factory for custom surveys."""
```

---

### 2.7 `validate.py` — Validation Suite

```python
def validate_shmr(params: SHMRParams) -> bool:
    """Check SHMR shape at z=0 and z=1. Return True if all checks pass."""

def validate_delta_sigma() -> bool:
    """Check ΔΣ for single NFW halo. Return True if order-of-magnitude correct."""

def validate_covariance() -> bool:
    """Check σ(ΔΣ) scaling with N_lens. Return True if consistent."""

def validate_fisher(fisher: np.ndarray, param_names: list[str]) -> bool:
    """Check positive-definiteness, condition number, and error ordering."""

def run_all_validations() -> dict[str, bool]:
    """Run all checks, print results, return dict of pass/fail."""
```

---

### 2.8 `plot_results.py` — Visualization

All plot functions take a `save_path` argument (default None = show interactively).

```python
def plot_shmr_validation(params: SHMRParams, save_path=None):
    """M*/Mh vs Mh at z=0, 0.5, 1.0, 2.0. Four curves, one panel."""

def plot_delta_sigma_model(params: SHMRParams, z: float, save_path=None):
    """ΔΣ(R) for multiple M* bins at given z."""

def plot_delta_sigma_with_errors(params, z, survey_configs, lensing_config, save_path=None):
    """ΔΣ(R) with error bars for different survey tiers. One panel per M* bin."""

def plot_fisher_comparison(results: dict, save_path=None):
    """Bar chart of σ(θ_i) across surveys."""

def plot_fisher_ellipses(fisher, param_names, pairs, save_path=None):
    """2D Fisher ellipses for specified parameter pairs."""

def plot_scaling(sweep_results: dict, param_name: str, save_path=None):
    """σ(param) vs sweep variable."""

def plot_two_regime_summary(results: dict, save_path=None):
    """Side-by-side: z=0 shape constraints (left) vs z-evolution constraints (right)."""

def plot_improvement_factor(results_baseline: dict, results_target: dict, save_path=None):
    """Ratio of σ(Stage-V) / σ(Stage-IV) per parameter."""
```

---

### 2.9 `shmr_fisher/config_io.py` — YAML Configuration Loader

Parses a YAML config file into a `RunConfig` dataclass that drives the entire forecast pipeline.

```python
@dataclass
class RunConfig:
    run_name: str                           # From filename stem or explicit override
    shmr_params: SHMRParams
    forecast_config: ForecastConfig
    lensing_config: LensingConfig
    nuisance_config: NuisanceConfig | None
    surveys: dict[str, SurveyConfig]        # 1–4 surveys
    output_dir: Path                        # outputs/{run_name}/

def load_run_config(yaml_path: str | Path) -> RunConfig:
    """Load and validate a YAML config file."""
```

The YAML schema supports z-dependent completeness via `log_Mstar_min_z_dep: {base, slope}`, which constructs a lambda `f(z) = base + slope * z` — avoiding `eval()` while covering all current use cases.

Validation rules:
- `surveys` section must have 1–4 entries.
- Each survey must have all required fields: `name`, `area_deg2`, `z_min`, `z_max`, `n_gal_total`, `log_Mstar_min`.
- Warns if `include_nuisance_params` is True but no `nuisance` section is provided.

### 2.10 `scripts/run_forecast.py` — Main Driver

```
Usage (YAML mode — preferred):
    uv run python scripts/run_forecast.py --config configs/default.yaml

Usage (legacy mode):
    uv run python scripts/run_forecast.py --surveys stage4_low_z stage5_wide
    uv run python scripts/run_forecast.py --systematics

Outputs (YAML mode):
    outputs/{run_name}/config.yaml           — copy of input config
    outputs/{run_name}/forecast_results.json  — all constraints
    outputs/{run_name}/forecast_results.npz   — Fisher matrices
    outputs/{run_name}/*.png                  — figures (300 dpi)
```

### 2.11 `scripts/run_sweep.py` — Parameter Sweep Driver

```
Usage (YAML mode):
    uv run python scripts/run_sweep.py --config configs/default.yaml --survey-key stage4_low_z --param area_deg2 --values 1000,5000,10000,14000

Usage (legacy mode):
    uv run python scripts/run_sweep.py --base stage5_wide --param n_gal_total --values 1e6,5e6,1e7,5e7,1e8

Outputs:
    outputs/{run_name}/sweep_{param}_{base}.png
    outputs/{run_name}/sweep_{param}_{base}.json
```

### 2.12 `scripts/generate_science_figures.py` — Science Figure Generator

```
Usage:
    uv run python scripts/generate_science_figures.py --config configs/default.yaml
    uv run python scripts/generate_science_figures.py --config configs/default.yaml --skip-sweeps

Produces all science figures (ΔΣ with errors, two-regime summary, improvement
factor, parameter sweeps) in outputs/{run_name}/.
```

---

## 3. Data Flow

```
YAML config (or CLI args)
        │
        ▼
    config_io.load_run_config()  ──► RunConfig
        │
        ├── SHMRParams (with mass-dependent scatter option)
        ├── SurveyConfig(s) (1–4 surveys)
        ├── ForecastConfig (systematic floor, nuisance toggle)
        ├── LensingConfig, NuisanceConfig
        │
        ▼
    fisher.compute_fisher_matrix()
        │
        ├──► For each (z_bin, M*_bin):
        │       │
        │       ├── halo_model.delta_sigma_bin()  ──► ΔΣ(R) = 1-halo NFW + 2-halo
        │       ├── halo_model.effective_bias()    ──► b_eff scalar
        │       ├── halo_model.galaxy_number_density() ──► n_gal scalar
        │       │
        │       ├── covariance.lensing_covariance()  ──► σ²(ΔΣ) per R-bin
        │       │     + systematic floor: σ² += (f_sys × ΔΣ_fid)²
        │       ├── covariance.clustering_covariance() ──► σ²(b_eff), σ²(n_gal)
        │       │
        │       ├── fisher.compute_derivatives()   ──► Jacobian ∂obs/∂θ
        │       │     (SHMR params via finite diff + nuisance analytically)
        │       │
        │       └──► F += J^T C^{-1} J
        │
        ├──► Add nuisance priors: F[i,i] += 1/σ² for (shear_m, photo_dz)
        │
        ▼
    Fisher matrix F_{ij}
        │
        ├──► fisher.extract_shmr_constraints()  ──► σ(θ_i) SHMR-only
        ├──► fisher.marginalized_errors()       ──► σ(θ_i) all params
        ├──► fisher.fisher_ellipse()            ──► 2D contours
        └──► fisher.add_external_prior()        ──► F with external priors
```

---

## 4. Acceptance Tests

### AT-1: Round-trip SHMR
`mean_log_Mh(mean_log_Mstar(12.0, fiducial, 0), fiducial, 0) == 12.0` within 0.01 dex.

### AT-2: ΔΣ order of magnitude
`delta_sigma_nfw([0.5], 1e13, 0.3)` in range [10^13, 10^15] M☉/Mpc² (i.e., 10–1000 M☉/pc²).

### AT-3: Covariance scaling
`lensing_covariance(..., N_lens=1e6) / lensing_covariance(..., N_lens=1e4) ≈ 0.01` (ratio of variances).

### AT-4: Fisher positive-definiteness
All eigenvalues of the Fisher matrix must be > 0.

### AT-5: Error ordering
For every survey: `marginalized_errors >= conditional_errors` element-wise.

### AT-6: Stage ordering
For the shared 5 z=0 parameters: `σ(Stage-V) < σ(Stage-IV) < σ(Stage-III)`.

### AT-7: Sweep monotonicity
`σ(param)` should decrease monotonically as `n_gal_total` increases (all else fixed).

### AT-8: Binning consistency
Changing `dlog_Mstar` from 0.5 to 0.25 should not change total Fisher information by more than ~30% (finer bins add some information from better mass resolution, but should not change drastically).

---

## 5. Non-Functional Requirements

- **Performance:** A single survey forecast (5 z-bins × 6 M*-bins × 10 R-bins × 9 params × 2 derivative evaluations) should complete in < 5 minutes on a laptop. If slower, precompute ΔΣ on a (Mh, z) grid and interpolate.
- **Reproducibility:** All random seeds (if any) must be fixed. All results must be deterministic given the same config.
- **Portability:** Pure Python + numpy/scipy/colossus/matplotlib/astropy. No compiled extensions.
