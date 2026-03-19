# SHMR Forecast: Implementation Plan for Claude Code

## Goal

Build a Python-based Fisher forecast tool that predicts how well a spectroscopic survey (combined with Stage-IV imaging lensing) can constrain the stellar-halo mass relation (SHMR) — both the mean relation and its scatter — as a function of survey parameters (galaxy number density, redshift range, survey area).

This is a **first-order, relative forecast** — we want to compare survey configurations, not produce absolute error bars. We fix lensing assumptions (LSST-like) and vary only the spectroscopic survey properties. The tool should be **generic**: not tuned to any specific survey, but parameterized so that any Stage-III through Stage-V spectroscopic survey can be specified via a simple configuration.

---

## Architecture Overview

```
shmr_fisher/
├── config.py              # Survey configs, cosmology, SHMR fiducial params
├── shmr_model.py          # Parameterized SHMR (Moster+2013 form) with z-evolution
├── halo_model.py          # Halo model predictions: DeltaSigma, b_eff, n_gal
├── covariance.py          # Analytic covariance matrices for lensing + clustering
├── fisher.py              # Fisher matrix computation
├── survey_configs.py      # Predefined survey configurations to compare
├── run_forecast.py        # Main driver script
├── plot_results.py        # Visualization of Fisher constraints
└── README.md              # Documentation
```

---

## Phase 1: Foundation (Tasks 1–3)

### Task 1: Environment & Dependencies

**Install:**
```bash
pip install colossus numpy scipy matplotlib emcee astropy
```

**Key dependency: `colossus`** (Diemer 2018). This is a pure-Python package, no compilation needed, and provides:
- Cosmology (distances, volumes, critical density)
- Halo mass function (Tinker08)
- Halo bias (Tinker10)
- NFW profiles (density, enclosed mass, surface density, ΔΣ)
- Concentration–mass relation (Diemer19)

We do **NOT** use `halotools` (requires downloading N-body halo catalogs, ~GB-scale; overkill for an analytic Fisher forecast). We do **NOT** use `halomod` (less mature for g-g lensing, harder to customize the SHMR). We build the halo model prediction layer ourselves on top of `colossus` primitives.

**Rationale:** `colossus` gives us all the building blocks (HMF, bias, NFW profile, Σ, ΔΣ) without requiring simulation data. For a first-order analytic forecast, we do not need to populate mock catalogs. We compute everything via integrals over the halo mass function.

### Task 2: Configuration Module (`config.py`)

Define a dataclass or dictionary structure for:

```python
@dataclass
class Cosmology:
    """Fiducial cosmology. Use Planck18."""
    # Just set colossus cosmology: cosmology.setCosmology('planck18')

@dataclass
class SHMRParams:
    """
    Parameterize the mean SHMR using the Moster+2013 style double power-law
    for M*/Mh as a function of Mh, WITH redshift evolution:

        M*(Mh, z) = 2 * N(z) * Mh * [(Mh/M1(z))^(-beta(z)) + (Mh/M1(z))^(gamma(z))]^(-1)

    Redshift evolution follows Moster+2013 Eq. 11–14:
        log10(M1(z))    = log_M1_0    + nu_M1    * z/(1+z)
        N(z)            = N_0          + nu_N     * z/(1+z)
        beta(z)         = beta_0       + nu_beta  * z/(1+z)
        gamma(z)        = gamma_0      + nu_gamma * z/(1+z)

    The scatter sigma_logMs is kept constant with z (standard assumption;
    Behroozi+2019 find only weak evolution).

    z=0 parameters (fiducial from Moster+2013 Table 1):
        log_M1_0:   Characteristic halo mass at z=0           ~ 11.59
        N_0:        Normalization (peak M*/Mh) at z=0          ~ 0.0351
        beta_0:     Low-mass slope at z=0                      ~ 1.376
        gamma_0:    High-mass slope at z=0                     ~ 0.608

    Redshift evolution parameters:
        nu_M1:      Evolution of log_M1                        ~ 1.195
        nu_N:       Evolution of normalization                 ~ -0.0247
        nu_beta:    Evolution of low-mass slope                ~ -0.826
        nu_gamma:   Evolution of high-mass slope               ~ 0.329

    Scatter:
        sigma_logMs: Log-normal scatter in M* at fixed Mh [dex] ~ 0.15
    """
    # z=0 SHMR shape
    log_M1_0: float = 11.59
    N_0: float = 0.0351
    beta_0: float = 1.376
    gamma_0: float = 0.608
    # Redshift evolution
    nu_M1: float = 1.195
    nu_N: float = -0.0247
    nu_beta: float = -0.826
    nu_gamma: float = 0.329
    # Scatter (constant with z)
    sigma_logMs: float = 0.15

@dataclass
class SurveyConfig:
    """
    Spectroscopic survey parameters — fully generic.

    A survey is defined by its area, redshift range, and a description
    of the galaxy sample (number density + mass completeness) which CAN
    vary with redshift. For simplicity, we specify n_gal and log_Mstar_min
    per redshift bin. For a uniform survey, these are constant across bins.

    The stellar mass binning scheme is also configurable:
        - dlog_Mstar: width of stellar mass bins [dex]
        - log_Mstar_max: upper edge of highest bin [log10(Msun)]
    """
    name: str
    area_deg2: float                # Survey area [deg^2]
    z_min: float                    # Minimum redshift
    z_max: float                    # Maximum redshift
    n_gal_total: float              # Total galaxies in survey (for reference)
    log_Mstar_min: float            # Global stellar mass completeness limit [log10(Msun)]
    log_Mstar_max: float = 12.0     # Upper edge of highest bin
    dlog_Mstar: float = 0.5         # Stellar mass bin width [dex] (user-configurable)

    # Optional: per-z-bin overrides for mass completeness
    # If provided, log_Mstar_min_func(z) returns the mass limit at redshift z
    # Default: constant log_Mstar_min at all z
    log_Mstar_min_func: callable = None  # callable(z) -> log_Mstar_min(z)

    def get_log_Mstar_min(self, z):
        """Return the stellar mass completeness limit at redshift z."""
        if self.log_Mstar_min_func is not None:
            return self.log_Mstar_min_func(z)
        return self.log_Mstar_min

    def get_n_gal_per_deg2(self):
        """Derived: average number density per deg^2."""
        return self.n_gal_total / self.area_deg2

    def get_stellar_mass_bins(self, z=None):
        """
        Return stellar mass bin edges for a given redshift.
        Bin edges run from get_log_Mstar_min(z) to log_Mstar_max
        in steps of dlog_Mstar.

        Returns: list of (lo, hi) tuples in log10(Msun)
        """
        mmin = self.get_log_Mstar_min(z) if z is not None else self.log_Mstar_min
        edges = np.arange(mmin, self.log_Mstar_max + 0.01, self.dlog_Mstar)
        return [(edges[i], edges[i+1]) for i in range(len(edges)-1)]


@dataclass
class LensingConfig:
    """
    Fixed lensing survey parameters (LSST-like / Euclid-like).
    These are held constant across all forecast comparisons.
    """
    n_source_per_arcmin2: float = 25.0   # Effective source density
    sigma_e: float = 0.26               # Shape noise per component
    z_source_median: float = 1.0         # Median source redshift
    # Source redshift distribution: use Smail-type n(z) ~ z^2 exp(-(z/z0)^1.5)
    # with z0 ≈ z_median / 1.41

@dataclass 
class ForecastConfig:
    """
    Master configuration controlling the forecast behavior.
    """
    # Radial bins for ΔΣ measurement
    R_min_Mpc: float = 0.1          # Minimum projected radius [physical Mpc]
    R_max_Mpc: float = 30.0         # Maximum projected radius [physical Mpc]
    n_R_bins: int = 10              # Number of log-spaced radial bins

    # Redshift binning
    dz: float = 0.2                 # Default redshift bin width
    # (can be overridden per survey)

    # Fisher matrix settings
    frac_step: float = 0.01         # Fractional step for numerical derivatives
    
    # Which SHMR parameters to vary in the Fisher matrix
    # (allows freezing evolution params for low-z-only surveys)
    vary_z_evolution: bool = True    # If False, fix nu_* params

    # Halo mass integration range
    log_Mh_min: float = 10.0        # Minimum halo mass [log10(Msun)]
    log_Mh_max: float = 15.5        # Maximum halo mass [log10(Msun)]
    n_Mh_bins: int = 200            # Number of log-spaced halo mass bins
```

**Key design principle:** All binning choices (stellar mass bin width, redshift bin width, radial bins) are configurable parameters, not hard-coded. The default `dlog_Mstar = 0.5` is sensible for a first pass, but the user can set `dlog_Mstar = 0.25` for finer exploration of the low-mass regime, or `dlog_Mstar = 1.0` for a quick coarse run.

### Task 3: SHMR Model (`shmr_model.py`)

Implement the mean SHMR with redshift evolution:

```python
def shmr_params_at_z(params: SHMRParams, z: float):
    """
    Evaluate the redshift-dependent SHMR parameters at redshift z.
    Following Moster+2013 Eq. 11–14:
        log10(M1(z))  = log_M1_0  + nu_M1  * z/(1+z)
        N(z)          = N_0       + nu_N   * z/(1+z)
        beta(z)       = beta_0    + nu_beta * z/(1+z)
        gamma(z)      = gamma_0   + nu_gamma * z/(1+z)

    Returns: (log_M1, N, beta, gamma) at the given redshift
    """
    zfac = z / (1.0 + z)
    log_M1 = params.log_M1_0 + params.nu_M1 * zfac
    N      = params.N_0      + params.nu_N   * zfac
    beta   = params.beta_0   + params.nu_beta * zfac
    gamma  = params.gamma_0  + params.nu_gamma * zfac
    return log_M1, N, beta, gamma


def mean_Mstar_over_Mh(log_Mh, params: SHMRParams, z: float):
    """
    Mean M*/Mh ratio at given Mh and redshift.

    Moster+2013 Eq. 2:
        M*/Mh = 2N * [(Mh/M1)^(-beta) + (Mh/M1)^(gamma)]^(-1)

    Returns: M*/Mh (dimensionless ratio)
    """
    log_M1, N, beta, gamma = shmr_params_at_z(params, z)
    x = 10**(log_Mh - log_M1)  # Mh / M1
    return 2.0 * N * (x**(-beta) + x**gamma)**(-1)


def mean_log_Mstar(log_Mh, params: SHMRParams, z: float):
    """
    Mean log10(M*) at given log10(Mh) and redshift.
    Returns: log10(M*) in solar masses
    """
    ratio = mean_Mstar_over_Mh(log_Mh, params, z)
    return log_Mh + np.log10(ratio)


def mean_log_Mh(log_Mstar, params: SHMRParams, z: float):
    """
    Inverse SHMR: mean log10(Mh) at given log10(M*) and redshift.
    Use scipy root-finding (brentq) on the forward relation.
    """

def phi_Mstar_given_Mh(log_Mstar, log_Mh, params: SHMRParams, z: float):
    """
    P(log M* | Mh, z): log-normal distribution centered on
    mean_log_Mstar(Mh, z) with dispersion sigma_logMs.
    """
```

**Validation step:** Plot the SHMR (M*/Mh vs Mh) at z=0, 0.5, 1.0, 2.0 and verify:
- Peak efficiency ~3–5% at Mh ~ 10^12 Msun at z=0
- Mild evolution: peak shifts to slightly higher Mh at higher z
- Low-mass slope steepens with z (beta increases)
- Compare visually with Moster+2013 Fig. 4 or Behroozi+2019 Fig. 8

---

## Phase 2: Halo Model Predictions (Tasks 4–6)

### Task 4: HOD from SHMR (`halo_model.py` — Part A)

Given the SHMR, derive the Halo Occupation Distribution:

```python
def N_cen(log_Mh, log_Mstar_lo, log_Mstar_hi, params: SHMRParams, z: float):
    """
    Mean number of central galaxies in stellar mass bin [M*_lo, M*_hi].

    N_cen(Mh | M*_lo, M*_hi) = integral of P(log M* | Mh) from M*_lo to M*_hi
                              = 0.5 * [erf(x_hi) - erf(x_lo)]
    where:
        x = (log_Mstar - mean_log_Mstar(Mh, z)) / (sqrt(2) * sigma_logMs)
    """

def N_sat(log_Mh, log_Mstar_lo, log_Mstar_hi, params: SHMRParams, z: float):
    """
    Mean number of satellites in stellar mass bin.

    Simplest approach: use Leauthaud+2011 / Zheng+2005 style:
        N_sat = N_cen_thresh * ((Mh - M0) / M1_sat)^alpha_sat

    For a first-order forecast, the satellite parameters are FIXED
    (not varied in the Fisher matrix):
        M0 ≈ 0 (simplification)
        M1_sat ≈ 17 * M_min  (empirical ratio from Leauthaud+2011)
        alpha_sat ≈ 1.0
    where M_min is the halo mass where N_cen(threshold) = 0.5.

    The satellite fraction is typically 5–15% and contributes mainly
    at small R in the lensing signal. For a relative forecast, fixing
    this is acceptable.
    """
```

**Important simplification:** For a relative forecast, we do not need to perfectly model satellites. The key signal is from centrals. Satellites affect the small-scale 1-halo term, but at first order, fixing the satellite prescription and varying only the SHMR parameters is sufficient.

### Task 5: Model Observables — ΔΣ(R) (`halo_model.py` — Part B)

The galaxy-galaxy lensing signal for a stellar mass bin [M*_lo, M*_hi] at redshift z:

```
ΔΣ(R | bin, z) = ∫ dMh  dn/dMh  [N_cen(Mh|bin,z) * ΔΣ_NFW(R, Mh, c(Mh,z), z)
                                   + N_sat(Mh|bin,z) * ΔΣ_sat(R, Mh, z)]  /  n_gal(bin, z)
```

where:
- `dn/dMh` = halo mass function → use `colossus.lss.mass_function`
- `ΔΣ_NFW(R, Mh, c)` = excess surface mass density of an NFW profile
  → `colossus` provides `profile_nfw.py` with `surfaceDensity()` and `deltaSigma()` methods
- `c(Mh, z)` = concentration–mass relation → use `colossus.halo.concentration`
- `n_gal(bin, z)` = total galaxy number density in this bin = ∫ dMh (dn/dMh) [N_cen + N_sat]

**For `colossus` NFW ΔΣ:**
```python
from colossus.halo import profile_nfw
from colossus.halo import concentration

def delta_sigma_nfw(R_Mpc, Mh, z):
    """
    Compute ΔΣ(R) for a single NFW halo.

    Args:
        R_Mpc: projected radii in physical Mpc (array)
        Mh: halo mass in Msun (M200m)
        z: redshift

    Returns:
        ΔΣ in Msun/Mpc^2 (convert to Msun/pc^2 later if needed)
    """
    c = concentration.concentration(Mh, '200m', z, model='diemer19')
    p = profile_nfw.NFWProfile(M=Mh, c=c, z=z, mdef='200m')
    # colossus uses comoving kpc/h internally — be very careful with units
    # surfaceDensity returns Sigma(R) in Msun/kpc^2 (comoving)
    # deltaSigma returns ΔΣ(R) = Sigma_mean(<R) - Sigma(R)
    R_kpc_h = R_Mpc * 1e3 * cosmo.h  # convert physical Mpc to comoving kpc/h
    ds = p.deltaSigma(R_kpc_h)  # Msun / kpc^2 (comoving, h units)
    # Convert to physical Msun/Mpc^2:
    ds_physical = ds * (1e3)**2 * cosmo.h**2 * (1 + z)**2
    return ds_physical
```

**⚠️ UNIT WARNING:** `colossus` profiles work in comoving kpc/h. You MUST convert carefully. The `deltaSigma()` method of the NFW profile class returns ΔΣ in comoving (Msun h / kpc^2). Converting to physical Msun/pc^2 requires multiplying by h * (1+z)^2 / (1e3)^2 ... actually, check the colossus docs carefully. This is the most common source of bugs. **Write a unit test comparing the output against a known analytic NFW ΔΣ for M=10^14, c=5, z=0.**

### Task 6: Model Observables — Clustering Summary Statistics (`halo_model.py` — Part C)

For the first-order forecast, we use a **simplified clustering description** rather than computing full w_p(r_p). Per stellar mass bin at redshift z, we compute two summary statistics:

**Effective halo bias:**
```
b_eff(bin, z) = ∫ dMh (dn/dMh) * [N_cen(Mh|bin,z) + N_sat(Mh|bin,z)] * b(Mh, z) / n_gal(bin, z)
```
where b(Mh, z) is the halo bias → `colossus.lss.bias.haloBias()`.

**Galaxy number density:**
```
n_gal(bin, z) = ∫ dMh (dn/dMh) * [N_cen(Mh|bin,z) + N_sat(Mh|bin,z)]
```

These two quantities per bin capture the dominant clustering information:
- b_eff directly constrains the mean halo mass of the bin
- n_gal constrains the integral of the HOD (i.e., the normalization/shape of the SHMR)

**Rationale for simplification:** Computing full w_p(r_p) via the halo model requires implementing the 1-halo convolution of NFW profiles in Fourier space, then Hankel-transforming back. This is doable but fiddly and bug-prone. For a relative forecast comparing survey configurations, (b_eff, n_gal) captures the dominant constraining power. Full w_p is a stretch goal.

---

## Phase 3: Covariance & Fisher Matrix (Tasks 7–9)

### Task 7: Lensing Covariance (`covariance.py`)

The ΔΣ measurement error in radial bin i for a lens sample in stellar mass bin α at redshift z:

```
σ²(ΔΣ_i) = σ_SN² + σ_LSS²
```

**Shape noise term (dominates at small R):**
```
σ_SN²(ΔΣ, R_i) = σ_e² * <Σ_crit^{-2}>^{-2} / (2π R_i ΔR_i * n_source_eff * N_lens / A_survey)
```

For the first-order forecast, the key scaling is:

```
σ(ΔΣ) ∝ σ_e * Σ_crit_eff / sqrt(N_lens * n_source_eff * A_annulus)
```

The number of source galaxies behind each lens in the annulus is:
```
N_source_per_lens_per_annulus = n_source_eff(z_l) * A_annulus
```

where n_source_eff is the effective source density at redshifts > z_lens.

**Σ_crit calculation:**
```python
def sigma_crit(z_lens, z_source):
    """
    Critical surface mass density.
    Σ_crit = c² / (4πG) * D_s / (D_l * D_ls)

    Use colossus or astropy for angular diameter distances.
    Returns in physical Msun/Mpc^2 (or Msun/pc^2).
    """
```

For a source redshift distribution n(z_s), use the effective inverse critical surface density:
```
<Σ_crit^{-1}> = ∫ dz_s  n(z_s)  Σ_crit^{-1}(z_l, z_s)  [for z_s > z_l]
```

**For the diagonal covariance (sufficient for first-order Fisher):**
```python
def cov_delta_sigma(R_bins, z_lens, N_lens, lensing_config, survey_area_sr):
    """
    Diagonal covariance matrix for ΔΣ in radial bins.

    Returns: array of shape (N_R_bins,) giving σ²(ΔΣ) at each R.
    """
```

**Key scaling to verify:** For SDSS-like (N_lens ~ 10^5, n_s ~ 1/arcmin^2), σ(ΔΣ) ~ few Msun/pc^2 at R ~ 1 Mpc. For LSST lensing + DESI-like spec-z (N_lens ~ 10^6, n_s ~ 25/arcmin^2), should be ~10x better.

### Task 8: Clustering Covariance (`covariance.py` continued)

For the effective bias b_eff, the error is set by cosmic variance + shot noise on large-scale clustering:

```
σ²(b_eff) / b_eff² ≈ 1 / (n_gal * V_eff) * (1 + 1/(n_gal * P_gg(k_eff)))
```

**Practical simplified form:**
```
σ(b_eff) / b_eff ≈ 1 / sqrt(n_gal * V_eff)  [for n*P >> 1, cosmic variance limit]
σ(b_eff) / b_eff ≈ 1 / sqrt(N_gal)  [for shot-noise limit]
```

where V_eff is the effective survey volume in the redshift bin, computed from:
```
V_eff = A_survey * ∫_{z_lo}^{z_hi} dV/dz dz    [comoving volume in the shell]
```

For the galaxy number density constraint:
```
σ(n_gal) / n_gal ≈ 1 / sqrt(N_gal)  [Poisson]
```

These two "observables" per bin (b_eff, n_gal) provide the clustering constraints.

### Task 9: Fisher Matrix (`fisher.py`)

**Two-regime approach for redshift evolution parameters:**

The SHMR has 9 parameters (4 z=0 shape + 4 evolution + 1 scatter). Not all surveys constrain all parameters equally. The forecast should handle this via the `ForecastConfig.vary_z_evolution` flag and by intelligently choosing which parameters to vary based on the survey's redshift range:

- **Low-z surveys (z < 0.4):** The z/(1+z) factor is < 0.29, so evolution parameters are weakly constrained. For such surveys, fix the evolution parameters (nu_M1, nu_N, nu_beta, nu_gamma) at their fiducial values and forecast only the z=0 shape + scatter (5 parameters). This avoids ill-conditioned Fisher matrices.

- **Wide-z surveys (spanning z = 0 to z > 0.5):** Vary all 9 parameters. The z-evolution is constrained by the leverage between low-z and high-z bins.

- **High-z only surveys (z > 0.5):** The z=0 parameters and evolution parameters are degenerate (they can trade off). Either fix z=0 parameters to external priors (from e.g., SDSS) or reparameterize to "parameters at z_pivot" + evolution.

```python
def get_varied_params(params: SHMRParams, forecast_config: ForecastConfig,
                      z_min: float, z_max: float) -> list:
    """
    Determine which SHMR parameters to vary in the Fisher matrix,
    based on the survey's redshift range and user configuration.

    Returns: list of (param_name, fiducial_value) tuples
    """
    # Always vary the z=0 shape and scatter
    varied = [
        ('log_M1_0', params.log_M1_0),
        ('N_0', params.N_0),
        ('beta_0', params.beta_0),
        ('gamma_0', params.gamma_0),
        ('sigma_logMs', params.sigma_logMs),
    ]

    # Vary evolution params only if survey spans enough redshift range
    # AND user has enabled it
    if forecast_config.vary_z_evolution and (z_max - z_min) > 0.3 and z_max > 0.4:
        varied += [
            ('nu_M1', params.nu_M1),
            ('nu_N', params.nu_N),
            ('nu_beta', params.nu_beta),
            ('nu_gamma', params.nu_gamma),
        ]

    return varied


def compute_fisher_matrix(
    shmr_params: SHMRParams,
    survey_config: SurveyConfig,
    lensing_config: LensingConfig,
    forecast_config: ForecastConfig,
):
    """
    Compute the Fisher matrix for SHMR parameters.

    Steps:
        1. Determine redshift bins from survey z range and forecast_config.dz
        2. For each z-bin, determine stellar mass bins from survey config
           (using survey_config.get_stellar_mass_bins(z) which respects
            the z-dependent mass completeness limit and dlog_Mstar)
        3. Determine which SHMR parameters to vary
        4. For each (z-bin, M*-bin) combination:
           a. Compute fiducial ΔΣ(R), b_eff, n_gal
           b. Compute covariance (lensing noise, clustering noise)
           c. Compute numerical derivatives wrt each varied parameter
           d. Accumulate into Fisher matrix
        5. Return Fisher matrix and metadata

    The Fisher matrix accumulates contributions from all bins:

    F_ij = Σ_{z-bins} Σ_{M*-bins} [
        Σ_{R-bins} (∂ΔΣ/∂θ_i) * (1/σ²(ΔΣ)) * (∂ΔΣ/∂θ_j)
      + (∂b_eff/∂θ_i) * (1/σ²(b_eff)) * (∂b_eff/∂θ_j)
      + (∂n_gal/∂θ_i) * (1/σ²(n_gal)) * (∂n_gal/∂θ_j)
    ]

    Derivatives are computed numerically via central finite differences:
        ∂ΔΣ/∂θ_i ≈ [ΔΣ(θ_i + δθ) - ΔΣ(θ_i - δθ)] / (2*δθ)

    Step sizes: use forecast_config.frac_step * |fiducial_value|.

    Returns:
        fisher_matrix: (N_params, N_params) array
        param_names: list of parameter names
        metadata: dict with bin info, fiducial values, etc.
    """
```

**Marginalized constraints:**
```python
def marginalized_errors(fisher_matrix):
    """
    Compute 1σ marginalized errors from Fisher matrix.
    σ(θ_i) = sqrt(F^{-1}_{ii})
    """
    cov = np.linalg.inv(fisher_matrix)
    return np.sqrt(np.diag(cov))
```

**Adding priors:** For nuisance parameters (if included), add to diagonal:
```
F_ii += 1/σ_prior²
```

For evolution parameters, one can optionally add external priors (e.g., from SDSS constraints on z=0 SHMR):
```python
def add_external_prior(fisher, param_names, param_name, sigma_prior):
    """Add a Gaussian prior on a single parameter."""
    idx = param_names.index(param_name)
    fisher[idx, idx] += 1.0 / sigma_prior**2
```

---

## Phase 4: Survey Comparison & Visualization (Tasks 10–12)

### Task 10: Survey Configurations (`survey_configs.py`)

Define several representative survey configurations spanning Stage-III through Stage-V. These are **generic archetypes**, not tied to specific survey names. The user can modify any parameter.

```python
surveys = {
    # --- Stage III: existing spectroscopic surveys ---
    "stage3_shallow_wide": SurveyConfig(
        name="Stage-III Shallow Wide",
        area_deg2=7500,
        z_min=0.02, z_max=0.2,
        n_gal_total=700_000,
        log_Mstar_min=9.5,
    ),

    # --- Stage IV: current-generation surveys ---
    "stage4_low_z": SurveyConfig(
        name="Stage-IV Low-z",
        area_deg2=14000,
        z_min=0.05, z_max=0.4,
        n_gal_total=10_000_000,
        log_Mstar_min=9.0,
    ),
    "stage4_high_z": SurveyConfig(
        name="Stage-IV High-z",
        area_deg2=14000,
        z_min=0.4, z_max=1.0,
        n_gal_total=5_000_000,
        log_Mstar_min=10.8,
    ),

    # --- Stage V: next-generation wide-field spectroscopic surveys ---
    "stage5_wide": SurveyConfig(
        name="Stage-V Wide",
        area_deg2=10000,
        z_min=0.05, z_max=1.0,
        n_gal_total=50_000_000,
        log_Mstar_min=8.5,
        # Mass completeness degrades with z:
        log_Mstar_min_func=lambda z: 8.5 + 1.0 * z,  # 8.5 at z=0, 9.5 at z=1
    ),
    "stage5_deep": SurveyConfig(
        name="Stage-V Deep",
        area_deg2=3000,
        z_min=0.1, z_max=1.5,
        n_gal_total=20_000_000,
        log_Mstar_min=8.0,
        log_Mstar_min_func=lambda z: 8.0 + 1.0 * z,  # 8.0 at z=0, 9.5 at z=1.5
    ),
}

# --- Convenience function for custom survey exploration ---
def make_custom_survey(name, area_deg2, z_min, z_max, n_gal_total,
                       log_Mstar_min, dlog_Mstar=0.5, **kwargs):
    """
    Quick factory for creating a custom survey config.
    Allows easy parameter sweeps.
    """
    return SurveyConfig(
        name=name,
        area_deg2=area_deg2,
        z_min=z_min, z_max=z_max,
        n_gal_total=n_gal_total,
        log_Mstar_min=log_Mstar_min,
        dlog_Mstar=dlog_Mstar,
        **kwargs,
    )
```

### Task 11: Run Forecast (`run_forecast.py`)

Main driver:

```python
def main():
    """
    For each survey config:
        1. Determine stellar mass bins (configurable width, z-dependent floor)
        2. Determine redshift bins
        3. Auto-select which SHMR params to vary based on z-range
        4. Compute Fisher matrix
        5. Extract marginalized errors on SHMR params
        6. Store results

    Then compare surveys.
    """
    forecast_config = ForecastConfig()
    lensing_config = LensingConfig()
    shmr_params = SHMRParams()

    results = {}
    for name, survey in surveys.items():
        fisher, param_names, metadata = compute_fisher_matrix(
            shmr_params=shmr_params,
            survey_config=survey,
            lensing_config=lensing_config,
            forecast_config=forecast_config,
        )

        errors = marginalized_errors(fisher)
        results[name] = {
            'fisher': fisher,
            'param_names': param_names,
            'errors': dict(zip(param_names, errors)),
            'metadata': metadata,
        }

    # Save results
    save_results(results, 'forecast_results.json')

    return results


def parameter_sweep(base_survey, sweep_param, sweep_values,
                    forecast_config=None, lensing_config=None, 
                    shmr_params=None):
    """
    Utility for exploring how constraints scale with a single survey parameter.

    Example usage:
        # How do constraints improve with total galaxy count?
        results = parameter_sweep(
            base_survey=surveys['stage5_wide'],
            sweep_param='n_gal_total',
            sweep_values=[1e6, 5e6, 1e7, 5e7, 1e8],
        )

        # How do constraints improve with mass completeness depth?
        results = parameter_sweep(
            base_survey=surveys['stage5_wide'],
            sweep_param='log_Mstar_min',
            sweep_values=[9.5, 9.0, 8.5, 8.0, 7.5],
        )

        # How do constraints improve with finer stellar mass bins?
        results = parameter_sweep(
            base_survey=surveys['stage5_wide'],
            sweep_param='dlog_Mstar',
            sweep_values=[1.0, 0.5, 0.25],
        )
    """
```

### Task 12: Visualization (`plot_results.py`)

Produce the following diagnostic and science plots:

1. **SHMR validation plot:** M*/Mh vs Mh at z=0, 0.5, 1.0 with the fiducial model. Show the redshift evolution clearly.

2. **Signal plot:** ΔΣ(R) for a few representative stellar mass bins at z=0.3. Show the model prediction and the expected error bars for each survey tier.

3. **Fisher constraint comparison (main result):**
   - Bar chart or table showing σ(θ_i) for each survey, for the shared parameters (z=0 shape + scatter). For surveys that also constrain evolution, show those separately.
   - 2D Fisher ellipses for key parameter pairs (e.g., log_M1_0 vs sigma_logMs, beta_0 vs gamma_0).

4. **Scaling plots:** How do constraints on key parameters scale with:
   - N_gal_total (at fixed area, depth)
   - Survey area (at fixed n_gal/deg^2)
   - log_Mstar_min (depth / mass completeness limit)
   - z_max (redshift reach)

5. **Two-regime summary figure:** Side-by-side panels showing:
   - Left: constraints on z=0 SHMR shape + scatter from low-z surveys (dwarf-galaxy focused; all surveys provide these)
   - Right: constraints on z-evolution parameters from wide-z surveys (only surveys with z_max > 0.4 and Δz > 0.3)
   This highlights the complementarity: deep low-z surveys nail the z=0 SHMR for dwarfs, while wide+high-z surveys constrain how it evolves.

6. **Improvement factor plot:** Ratio of Stage-V constraints to Stage-IV constraints, as a function of SHMR parameter. This is the main "science case" figure.

---

## Critical Implementation Notes

### Unit Conventions
- **Halo masses:** M200m in Msun (NOT Msun/h). Colossus can work in both; be explicit.
- **Distances:** Physical Mpc for ΔΣ radial bins. Colossus profiles use comoving kpc/h internally.
- **ΔΣ units:** Msun/pc^2 is the standard observational unit. Compute in Msun/Mpc^2 internally, then convert.
- **Set `colossus` cosmology once at the top:** `cosmology.setCosmology('planck18')`

### Numerical Stability
- The HMF integral is over ~5 orders of magnitude in mass. Use log-spaced integration (e.g., `np.logspace(10, 15.5, 200)` in Msun).
- NFW ΔΣ has an analytic form — use `colossus`'s implementation, don't re-derive.
- Finite difference derivatives: use step sizes of ~1–2% of the parameter value. Check convergence by comparing 1% vs 2% steps.
- For the Fisher matrix with 9 parameters, check the condition number. If > 10^10, the matrix is ill-conditioned — likely means evolution parameters are unconstrained. Fall back to fixing them.

### Validation Checkpoints
After each phase, verify:

1. **After Task 3:** Plot SHMR at z=0 and z=1. Peak efficiency should be ~3–5% at Mh ~ 10^12 at z=0, shifting mildly to higher mass at higher z.
2. **After Task 5:** Plot ΔΣ(R) for a single halo of M=10^{13} Msun at z=0.3. Should peak at ~10–50 Msun/pc^2 at R ~ 0.5 Mpc and decline as ~R^{-1} at large R.
3. **After Task 7:** Verify ΔΣ error bars scale as expected: σ ∝ 1/√(N_lens * n_source).
4. **After Task 9:** Check that Fisher matrix is positive definite. Check that doubling N_gal roughly halves the errors (in the shape-noise-dominated regime).
5. **After Task 10:** Sanity check: Stage-V should give significantly tighter constraints than Stage-III. If not, something is wrong.

### Known Simplifications (acceptable for first order)
- No cross-bin covariance (stellar mass bins treated as independent)
- No off-diagonal ΔΣ covariance (radial bins treated as independent; see SPEC.md Section 2.4.1 Feature 3 for feasibility assessment)
- Satellite prescription fixed, not varied
- ~~No shear calibration or photo-z systematics~~ **IMPLEMENTED:** optional nuisance parameter marginalization for shear calibration bias (m) and source photo-z bias (Δz_s) with configurable Gaussian priors. Toggle: `ForecastConfig.include_nuisance_params`.
- ~~No systematic error floor~~ **IMPLEMENTED:** optional fractional systematic floor on ΔΣ covariance. Toggle: `ForecastConfig.systematic_floor_fraction`.
- Single effective source redshift per lens bin (vs. full n(z_s) integration)
- Clustering info simplified to (b_eff, n_gal) rather than full w_p(r_p)
- Scatter assumed constant with z (weak assumption; Behroozi+2019 finds little evolution)

---

## Timeline for 24-Hour Execution

| Time | Task | Deliverable |
|------|------|-------------|
| 0–2h | Tasks 1–2: Setup, config | Working environment, all config dataclasses |
| 2–5h | Task 3: SHMR model | Validated SHMR plot at multiple redshifts |
| 5–10h | Tasks 4–5: HOD + ΔΣ | Working ΔΣ prediction with validation plot |
| 10–13h | Task 6: Clustering summary | b_eff and n_gal computation, validated |
| 13–17h | Tasks 7–8: Covariance | Validated error bars on ΔΣ and clustering |
| 17–20h | Task 9: Fisher matrix | Working Fisher computation with auto param selection |
| 20–22h | Tasks 10–11: Survey comparison | Results for all survey configs + parameter sweep |
| 22–24h | Task 12: Plots | All diagnostic and science figures |

### If things go wrong (fallback plan):
- If `colossus` NFW ΔΣ units are problematic → implement analytic NFW ΔΣ from Bartelmann (1996) / Wright & Brainerd (2000). The closed-form expressions are well-known.
- If full halo model integral is too slow → precompute ΔΣ on a grid of (Mh, z, R) and interpolate.
- If Fisher matrix is singular → check condition number; fix evolution parameters and redo with 5-param model.
- If redshift evolution makes derivatives numerically unstable → increase step size to 2–5% for nu_* parameters.

---

## Stretch Goals (if time permits)
1. Add w_p(r_p) as a full observable (replace b_eff/n_gal with radial bins)
2. Add nuisance parameters (shear calibration, photo-z bias) with Gaussian priors
3. Allow mass-dependent scatter σ_logMs(Mh)
4. Compare Fisher forecast with a quick MCMC on mock data (using `emcee`)
5. Produce a Jupyter notebook version for interactive exploration
6. Allow user to add external priors (e.g., SDSS z=0 SHMR constraints) to high-z-only forecasts

---

## References
- Moster, Naab & White (2013), MNRAS, 428, 3121 — SHMR parameterization with z-evolution
- Behroozi, Conroy & Wechsler (2010), ApJ, 717, 379 — SHMR analysis
- Behroozi et al. (2019), MNRAS, 488, 3143 — UniverseMachine
- Leauthaud et al. (2012), ApJ, 744, 159 — Joint lensing+clustering framework
- Zu & Mandelbaum (2015), MNRAS, 454, 1161 — iHOD framework
- Diemer (2018), ApJS, 239, 35 — COLOSSUS
- Wright & Brainerd (2000), ApJ, 534, 34 — Analytic NFW lensing
- Wechsler & Tinker (2018), ARA&A, 56, 435 — Galaxy-halo connection review
- Oguri & Takada (2011), PRD, 83, 023008 — Fisher forecast for lensing surveys
