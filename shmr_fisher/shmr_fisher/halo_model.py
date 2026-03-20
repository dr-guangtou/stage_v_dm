"""
Halo model predictions for galaxy-galaxy lensing and clustering.

This module connects the SHMR parameterization (shmr_model.py) to
observable quantities by integrating over the halo mass function.
It provides:

1. HOD functions: N_cen, N_sat — mean halo occupation for a stellar mass bin.
2. Lensing signal:
   - delta_sigma_nfw: single NFW halo (1-halo, analytic via colossus).
   - delta_sigma_matter / _sigma_matter: 2-halo term from matter correlation
     function (projected ξ_mm), cached per redshift.
   - delta_sigma_bin: bin-averaged ΔΣ = 1-halo + 2-halo, integrated over HMF.
3. Clustering: effective_bias, galaxy_number_density per stellar mass bin.

All halo model integrals use the colossus HMF (Tinker+2008) and integrate
over a log-spaced halo mass grid via _trapezoid. Halo masses are in
Msun (NOT Msun/h) throughout.

References
----------
Moster, Naab & White (2013), MNRAS, 428, 3121 — SHMR
Leauthaud et al. (2012), ApJ, 744, 159 — HOD from SHMR, satellite prescription
Tinker et al. (2008), ApJ, 688, 709 — halo mass function
Tinker et al. (2010), ApJ, 724, 878 — halo bias
Diemer & Kravtsov (2015); Diemer (2019) — concentration-mass relation
Wright & Brainerd (2000), ApJ, 534, 34 — analytic NFW lensing signal
"""

from __future__ import annotations

import warnings

import numpy as np

# NumPy 2.0 renamed np.trapz -> _trapezoid; support both
_trapezoid = getattr(np, "trapezoid", np.trapz)
from functools import lru_cache

import scipy.integrate
from colossus.cosmology import cosmology
from colossus.halo import concentration as colossus_conc
from colossus.halo import mass_defs as colossus_mdef
from colossus.halo import profile_nfw
from colossus.lss import bias as colossus_bias
from colossus.lss import mass_function as colossus_mf
# Matter correlation function is accessed via cosmology.correlationFunction()
from scipy.special import erf

from .config import SHMRParams
from .shmr_model import mean_log_Mstar, scatter_at_Mh


# ---------------------------------------------------------------------------
# HOD functions
# ---------------------------------------------------------------------------

def n_cen(
    log_Mh_grid: np.ndarray,
    log_Mstar_lo: float,
    log_Mstar_hi: float,
    params: SHMRParams,
    z: float,
) -> np.ndarray:
    """
    Mean number of central galaxies per halo in a stellar mass bin.

    Derived from the SHMR scatter: the probability that a central galaxy
    at fixed Mh has M* in [M*_lo, M*_hi] is the integral of the log-normal
    P(log M* | Mh) over that bin.

    N_cen(Mh) = 0.5 * [erf(x_hi) - erf(x_lo)]

    where x = (log M* - <log M*(Mh,z)>) / (sqrt(2) * sigma_logMs).

    Parameters
    ----------
    log_Mh_grid : array, shape (N_Mh,)
        Halo mass grid [log10(Msun)].
    log_Mstar_lo : float
        Lower edge of stellar mass bin [log10(Msun)].
    log_Mstar_hi : float
        Upper edge of stellar mass bin [log10(Msun)].
    params : SHMRParams
        SHMR parameters.
    z : float
        Redshift.

    Returns
    -------
    ncen : array, shape (N_Mh,)
        Mean central occupation (0 to 1).
    """
    mu = mean_log_Mstar(log_Mh_grid, params, z)
    # Mass-dependent scatter: sigma is an array matching log_Mh_grid
    sigma = scatter_at_Mh(log_Mh_grid, params)
    sqrt2_sigma = np.sqrt(2.0) * sigma

    x_lo = (log_Mstar_lo - mu) / sqrt2_sigma
    x_hi = (log_Mstar_hi - mu) / sqrt2_sigma

    return 0.5 * (erf(x_hi) - erf(x_lo))


def n_sat(
    log_Mh_grid: np.ndarray,
    log_Mstar_lo: float,
    log_Mstar_hi: float,
    params: SHMRParams,
    z: float,
) -> np.ndarray:
    """
    Mean number of satellite galaxies per halo in a stellar mass bin.

    Uses a simplified Leauthaud+2011/Zheng+2005 prescription with FIXED
    satellite parameters (not varied in the Fisher matrix):

        N_sat(Mh) = N_cen_thresh(Mh) * (Mh / M1_sat)^alpha_sat

    where:
        - N_cen_thresh: central occupation for a THRESHOLD sample
          (all galaxies above log_Mstar_lo), used as the satellite modulation
        - M1_sat = 17 * M_min, where M_min is the halo mass at which
          N_cen_thresh = 0.5 (empirical ratio from Leauthaud+2011)
        - alpha_sat = 1.0

    Parameters
    ----------
    log_Mh_grid : array, shape (N_Mh,)
        Halo mass grid [log10(Msun)].
    log_Mstar_lo : float
        Lower edge of stellar mass bin [log10(Msun)].
    log_Mstar_hi : float
        Upper edge of stellar mass bin [log10(Msun)].
    params : SHMRParams
        SHMR parameters.
    z : float
        Redshift.

    Returns
    -------
    nsat : array, shape (N_Mh,)
        Mean satellite occupation (>= 0).

    Notes
    -----
    The satellite fraction is typically 5-15% and contributes mainly at
    small projected radii. For a relative forecast, fixing this prescription
    is acceptable (see PLAN.md).
    """
    alpha_sat = 1.0
    # Ratio M1_sat / M_min from Leauthaud+2011 Eq. 5
    m1_over_mmin = 17.0

    # Threshold central occupation: P(M* > M*_lo | Mh)
    mu = mean_log_Mstar(log_Mh_grid, params, z)
    # Mass-dependent scatter: sigma is an array matching log_Mh_grid
    sigma = scatter_at_Mh(log_Mh_grid, params)
    sqrt2_sigma = np.sqrt(2.0) * sigma
    # N_cen_thresh = 0.5 * [1 + erf((mu - log_Mstar_lo) / (sqrt2*sigma))]
    n_cen_thresh = 0.5 * (1.0 + erf((mu - log_Mstar_lo) / sqrt2_sigma))

    # Find M_min: halo mass where N_cen_thresh = 0.5
    # This is where mu(M_min) = log_Mstar_lo, i.e., the mean SHMR hits
    # the bin lower edge. Use the inverse SHMR.
    from .shmr_model import mean_log_Mh
    try:
        log_Mmin = mean_log_Mh(log_Mstar_lo, params, z)
    except ValueError:
        # If the stellar mass is too low/high for the SHMR range,
        # satellites are negligible
        return np.zeros_like(log_Mh_grid)

    Mmin = 10.0**log_Mmin  # [Msun]
    M1_sat = m1_over_mmin * Mmin  # [Msun]
    Mh = 10.0**log_Mh_grid  # [Msun]

    # N_sat = N_cen_thresh * (Mh / M1_sat)^alpha_sat, only for Mh > M1_sat
    ratio = Mh / M1_sat
    nsat = n_cen_thresh * ratio**alpha_sat

    # Suppress satellites at low halo masses (below M_min)
    nsat[Mh < Mmin] = 0.0

    # Satellites in a BIN (not threshold): scale by the fraction of the
    # threshold that falls in this bin. Approximate: use the central
    # occupation ratio N_cen(bin) / N_cen(threshold).
    ncen_bin = n_cen(log_Mh_grid, log_Mstar_lo, log_Mstar_hi, params, z)
    # Avoid division by zero
    safe_thresh = np.maximum(n_cen_thresh, 1e-30)
    bin_fraction = ncen_bin / safe_thresh

    return nsat * bin_fraction


# ---------------------------------------------------------------------------
# Lensing signal: single NFW halo
# ---------------------------------------------------------------------------

def delta_sigma_nfw(
    R_Mpc: np.ndarray,
    Mh: float,
    z: float,
) -> np.ndarray:
    """
    Excess surface mass density DeltaSigma(R) for a single NFW halo.

    Uses colossus NFW profile with Diemer+2019 concentration-mass relation
    and the 200m mass definition.

    Parameters
    ----------
    R_Mpc : array
        Projected radii [physical Mpc].
    Mh : float
        Halo mass M200m [Msun] (NOT Msun/h).
    z : float
        Redshift.

    Returns
    -------
    ds : array
        DeltaSigma [Msun/Mpc^2]. Same shape as R_Mpc.

    Notes
    -----
    Colossus profiles work internally in comoving kpc/h. The unit
    conversion chain is:
        R [physical Mpc] -> R [comoving kpc/h] = R * 1e3 * h * (1+z)
    Colossus deltaSigma() returns in [h Msun / kpc^2] (comoving).
    To convert to [Msun / physical Mpc^2]:
        ds_phys = ds_colossus * h * (1+z)^2 * (1e3)^2
    """
    cosmo = cosmology.getCurrent()
    h = cosmo.H0 / 100.0

    # Colossus concentration expects M in Msun/h
    Mh_h = Mh * h  # [Msun/h]

    c = colossus_conc.concentration(Mh_h, '200m', z, model='diemer19')

    # Create NFW profile — colossus expects M in Msun/h
    prof = profile_nfw.NFWProfile(M=Mh_h, c=c, z=z, mdef='200m')

    # Convert R from physical Mpc to comoving kpc/h
    # physical -> comoving: multiply by (1+z)
    # Mpc -> kpc: multiply by 1e3
    # no-h -> h: multiply by h
    R_ckpc_h = R_Mpc * 1e3 * h * (1.0 + z)  # [comoving kpc/h]

    # colossus deltaSigma returns DeltaSigma in [h Msun / comoving kpc^2]
    ds_colossus = prof.deltaSigma(R_ckpc_h)

    # Convert [h Msun / kpc_com^2] to [Msun / physical Mpc^2]:
    #   * h to remove the h in numerator
    #   * (1e3)^2 to convert kpc^2 -> Mpc^2 in denominator
    #   * (1+z)^2 to convert comoving area -> physical area
    ds_physical_Mpc2 = ds_colossus * h * (1e3)**2 * (1.0 + z)**2

    return ds_physical_Mpc2


# ---------------------------------------------------------------------------
# 2-halo term: ΔΣ from correlated large-scale structure
# ---------------------------------------------------------------------------

# Cache ΔΣ_matter per (z, R-grid hash) to avoid recomputation.
# The matter contribution depends only on redshift and radial grid,
# not on halo mass or SHMR parameters.
_ds_matter_cache: dict[tuple[float, int], np.ndarray] = {}


def _sigma_matter(R_Mpc: float, z: float, pi_max: float = 100.0) -> float:
    """
    Surface mass density of the mean matter field at projected radius R.

    Σ_matter(R) = ρ_m(z) · ∫_{-π_max}^{+π_max} [1 + ξ_mm(√(R² + π²), z)] dπ

    In practice, the constant ρ_m · 2·π_max term cancels in ΔΣ = Σ̄(<R) - Σ(R),
    so we only integrate the ξ_mm part:

        Σ_mm(R) = ρ_m(z) · ∫_{-π_max}^{+π_max} ξ_mm(√(R² + π²), z) dπ

    Parameters
    ----------
    R_Mpc : float
        Projected radius [physical Mpc].
    z : float
        Redshift.
    pi_max : float
        Line-of-sight integration limit [physical Mpc].

    Returns
    -------
    sigma_mm : float
        Projected matter surface density contribution from ξ_mm
        [Msun / physical Mpc²].

    Notes
    -----
    colossus.lss.correlation.xi() expects r in comoving Mpc/h and
    returns the nonlinear matter correlation function (dimensionless).
    We convert R from physical Mpc to comoving Mpc/h for the call.
    """
    cosmo = cosmology.getCurrent()
    h = cosmo.H0 / 100.0
    # Mean matter density at z=0 in Msun/Mpc^3 (comoving)
    # colossus rho_m returns in Msun h^2 / kpc^3 (physical at z)
    # Use Omega_m * rho_crit(z=0) for comoving density
    rho_m_comoving = cosmo.rho_m(0.0) * 1e9  # [Msun h^2 / Mpc^3] → need h conversion
    # Actually: cosmo.rho_m(0) is in Msun h^2 / kpc^3
    # Convert to Msun / Mpc^3: * (1e3)^3 / h^2
    rho_m = cosmo.rho_m(0.0) * (1e3)**3 / h**2  # [Msun / comoving Mpc^3]

    # The projection integral sums along the line of sight in physical Mpc.
    # ξ_mm needs comoving Mpc/h: r_com = r_phys * (1+z) * h
    scale_to_ckpc_h = (1.0 + z) * h  # physical Mpc -> comoving Mpc/h

    def integrand(pi):
        r_3d_phys = np.sqrt(R_Mpc**2 + pi**2)  # physical Mpc
        r_3d_cmpc_h = r_3d_phys * scale_to_ckpc_h  # comoving Mpc/h
        # colossus correlationFunction expects r in comoving Mpc/h
        xi = cosmo.correlationFunction(r_3d_cmpc_h, z)
        return xi

    # Integrate from 0 to pi_max and double (symmetric)
    result, _ = scipy.integrate.quad(
        integrand, 0.0, pi_max, limit=100, epsrel=1e-3,
    )
    sigma_mm = 2.0 * rho_m * result  # factor 2 for symmetry

    return sigma_mm


def delta_sigma_matter(R_Mpc: np.ndarray, z: float) -> np.ndarray:
    """
    Excess surface mass density of the mean matter field: ΔΣ_mm(R).

    ΔΣ_mm(R) = Σ̄_mm(<R) - Σ_mm(R)

    where Σ̄_mm(<R) = (2/R²) ∫₀ᴿ Σ_mm(R') R' dR' is the mean surface
    density within projected radius R.

    This quantity is independent of halo mass. The halo-specific 2-halo
    lensing signal is obtained by multiplying by the halo bias:

        ΔΣ_2h(R, Mh, z) = b_h(Mh, z) · ΔΣ_mm(R, z)

    Results are cached per (z, R-grid) to avoid recomputation during
    Fisher derivative evaluations (which vary SHMR params but not z
    or R).

    Parameters
    ----------
    R_Mpc : array, shape (N_R,)
        Projected radii [physical Mpc].
    z : float
        Redshift.

    Returns
    -------
    ds_mm : array, shape (N_R,)
        ΔΣ of the matter field [Msun / physical Mpc²].
    """
    # Check cache
    cache_key = (round(z, 6), hash(R_Mpc.tobytes()))
    if cache_key in _ds_matter_cache:
        return _ds_matter_cache[cache_key]

    R_Mpc = np.atleast_1d(R_Mpc)

    # Compute Σ_mm on a fine grid from R_min/10 to R_max for accurate
    # mean-interior integration.
    R_min = R_Mpc.min()
    R_max = R_Mpc.max()
    R_fine = np.logspace(np.log10(R_min * 0.1), np.log10(R_max * 1.01), 80)

    sigma_fine = np.array([_sigma_matter(r, z) for r in R_fine])

    # Compute Σ̄_mm(<R) via cumulative integration:
    # Σ̄(<R) = (1/πR²) ∫₀ᴿ Σ(R') 2πR' dR' = (2/R²) ∫₀ᴿ Σ(R') R' dR'
    # Use cumulative trapezoid on the fine grid
    cumul_integrand = sigma_fine * R_fine  # Σ(R') · R'
    cumul_integral = np.zeros_like(R_fine)
    for i in range(1, len(R_fine)):
        cumul_integral[i] = _trapezoid(
            cumul_integrand[:i + 1], R_fine[:i + 1],
        )
    sigma_bar_fine = 2.0 * cumul_integral / R_fine**2
    # Fix R_fine[0]: use L'Hôpital / extrapolation
    sigma_bar_fine[0] = sigma_fine[0]

    # Interpolate Σ_mm and Σ̄_mm to the requested R grid
    from scipy.interpolate import interp1d
    sigma_interp = interp1d(
        np.log10(R_fine), sigma_fine, kind='cubic',
        fill_value='extrapolate',
    )
    sigma_bar_interp = interp1d(
        np.log10(R_fine), sigma_bar_fine, kind='cubic',
        fill_value='extrapolate',
    )

    sigma_at_R = sigma_interp(np.log10(R_Mpc))
    sigma_bar_at_R = sigma_bar_interp(np.log10(R_Mpc))

    ds_mm = sigma_bar_at_R - sigma_at_R  # [Msun / Mpc^2]

    _ds_matter_cache[cache_key] = ds_mm
    return ds_mm


def clear_2halo_cache() -> None:
    """Clear the ΔΣ_matter cache (call when cosmology changes)."""
    _ds_matter_cache.clear()


# ---------------------------------------------------------------------------
# Lensing signal: bin-averaged
# ---------------------------------------------------------------------------

def delta_sigma_bin(
    R_Mpc: np.ndarray,
    log_Mstar_lo: float,
    log_Mstar_hi: float,
    params: SHMRParams,
    z: float,
    log_Mh_grid: np.ndarray | None = None,
    include_2halo: bool = True,
) -> np.ndarray:
    """
    Halo-model prediction for DeltaSigma(R) averaged over a stellar mass bin.

    Combines the 1-halo NFW term (per-halo lensing) with a 2-halo term
    (correlated large-scale structure):

        ΔΣ(R | bin, z) = ΔΣ_1h(R) + ΔΣ_2h(R)

    where:
        ΔΣ_1h = ∫ dMh (dn/dMh) · N_total(Mh) · ΔΣ_NFW(R,Mh) / n_gal
        ΔΣ_2h = ∫ dMh (dn/dMh) · N_total(Mh) · b_h(Mh) · ΔΣ_mm(R) / n_gal
               = b_eff · ΔΣ_mm(R)

    The satellite term uses the same NFW profile centered on the host halo
    (1-halo approximation; ignores off-centering and subhalo profiles).

    Integration is done via _trapezoid over a log-spaced halo mass grid.

    Parameters
    ----------
    R_Mpc : array, shape (N_R,)
        Projected radii [physical Mpc].
    log_Mstar_lo : float
        Lower edge of stellar mass bin [log10(Msun)].
    log_Mstar_hi : float
        Upper edge of stellar mass bin [log10(Msun)].
    params : SHMRParams
        SHMR parameters.
    z : float
        Redshift.
    log_Mh_grid : array or None
        Halo mass grid [log10(Msun)]. Default: np.linspace(10, 15.5, 200).
    include_2halo : bool
        Whether to include the 2-halo term. Default True.

    Returns
    -------
    ds : array, shape (N_R,)
        DeltaSigma [Msun/Mpc^2], including 1-halo + 2-halo.
    """
    if log_Mh_grid is None:
        log_Mh_grid = np.linspace(10.0, 15.5, 200)

    Mh_grid = 10.0**log_Mh_grid  # [Msun]

    cosmo = cosmology.getCurrent()
    h = cosmo.H0 / 100.0

    # Halo mass function: dn/dln(M) in [h^3 Mpc^{-3}]
    # colossus expects M in Msun/h
    Mh_grid_h = Mh_grid * h  # [Msun/h]
    # massFunction returns dn/dlnM in comoving (Mpc/h)^{-3}
    dndlnM = colossus_mf.massFunction(
        Mh_grid_h, z, mdef='200m', model='tinker08', q_out='dndlnM'
    )
    # Convert to dn/dMh in [Msun^{-1} Mpc^{-3}] (comoving, no h):
    # dn/dMh = dn/dlnM / Mh, and convert (Mpc/h)^{-3} -> Mpc^{-3}: * h^3
    dndMh = dndlnM * h**3 / Mh_grid  # [Msun^{-1} Mpc^{-3}]

    # HOD: centrals + satellites
    ncen = n_cen(log_Mh_grid, log_Mstar_lo, log_Mstar_hi, params, z)
    nsat = n_sat(log_Mh_grid, log_Mstar_lo, log_Mstar_hi, params, z)
    n_total = ncen + nsat

    # Galaxy number density: n_gal = integral dMh (dn/dMh) * (N_cen + N_sat)
    n_gal = _trapezoid(dndMh * n_total * Mh_grid * np.log(10.0),
                         log_Mh_grid)

    if n_gal <= 0:
        warnings.warn(
            f"Zero galaxy density in bin [{log_Mstar_lo}, {log_Mstar_hi}] at z={z}. "
            "Returning zeros for DeltaSigma.",
            stacklevel=2,
        )
        return np.zeros_like(R_Mpc)

    # --- 1-halo term: NFW ---
    R_Mpc = np.atleast_1d(R_Mpc)
    ds_grid = np.zeros((len(log_Mh_grid), len(R_Mpc)))

    for i, (lmh, mh) in enumerate(zip(log_Mh_grid, Mh_grid)):
        if n_total[i] < 1e-10:
            continue  # Skip halos with negligible occupation
        ds_grid[i, :] = delta_sigma_nfw(R_Mpc, mh, z)

    # Weighted average: integral dMh (dn/dMh) * N_total * DS_NFW / n_gal
    # Use log-space integration: integral = integral d(log Mh) * Mh * ln(10) * ...
    integrand = (dndMh * n_total * Mh_grid * np.log(10.0))[:, np.newaxis] * ds_grid
    ds_1h = _trapezoid(integrand, log_Mh_grid, axis=0) / n_gal

    if not include_2halo:
        return ds_1h

    # --- 2-halo term: b_eff * ΔΣ_matter ---
    # ΔΣ_matter depends only on (R, z), cached per redshift.
    ds_mm = delta_sigma_matter(R_Mpc, z)

    # HOD-weighted effective bias for 2-halo weighting.
    # b_eff = ∫ dMh (dn/dMh) * N_total * b_h / n_gal
    # (Same quantity as effective_bias(), computed inline to reuse
    # the HMF and HOD arrays already in memory.)
    bh = colossus_bias.haloBias(
        Mh_grid_h, model='tinker10', z=z, mdef='200m',
    )
    weight = dndMh * n_total * Mh_grid * np.log(10.0)
    b_eff_val = _trapezoid(weight * bh, log_Mh_grid) / n_gal

    ds_2h = b_eff_val * ds_mm

    return ds_1h + ds_2h


# ---------------------------------------------------------------------------
# Clustering summary statistics
# ---------------------------------------------------------------------------

def effective_bias(
    log_Mstar_lo: float,
    log_Mstar_hi: float,
    params: SHMRParams,
    z: float,
    log_Mh_grid: np.ndarray | None = None,
) -> float:
    """
    HOD-weighted effective halo bias for a stellar mass bin.

    b_eff = integral dMh (dn/dMh) * [N_cen + N_sat] * b(Mh,z) / n_gal

    Uses colossus Tinker+2010 halo bias with 200m mass definition.

    Parameters
    ----------
    log_Mstar_lo : float
        Lower edge of stellar mass bin [log10(Msun)].
    log_Mstar_hi : float
        Upper edge of stellar mass bin [log10(Msun)].
    params : SHMRParams
        SHMR parameters.
    z : float
        Redshift.
    log_Mh_grid : array or None
        Halo mass grid [log10(Msun)]. Default: np.linspace(10, 15.5, 200).

    Returns
    -------
    b_eff : float
        Effective linear halo bias [dimensionless].
    """
    if log_Mh_grid is None:
        log_Mh_grid = np.linspace(10.0, 15.5, 200)

    Mh_grid = 10.0**log_Mh_grid
    cosmo = cosmology.getCurrent()
    h = cosmo.H0 / 100.0
    Mh_grid_h = Mh_grid * h

    # HMF: dn/dlnM in (Mpc/h)^{-3}
    dndlnM = colossus_mf.massFunction(
        Mh_grid_h, z, mdef='200m', model='tinker08', q_out='dndlnM'
    )
    dndMh = dndlnM * h**3 / Mh_grid  # [Msun^{-1} Mpc^{-3}]

    # HOD
    ncen = n_cen(log_Mh_grid, log_Mstar_lo, log_Mstar_hi, params, z)
    nsat = n_sat(log_Mh_grid, log_Mstar_lo, log_Mstar_hi, params, z)
    n_total = ncen + nsat

    # Halo bias: colossus expects M in Msun/h
    # Tinker+2010 bias for 200m definition
    bh = colossus_bias.haloBias(Mh_grid_h, model='tinker10', z=z, mdef='200m')

    # Integration weight: dMh = Mh * ln(10) * d(log Mh)
    weight = dndMh * n_total * Mh_grid * np.log(10.0)

    n_gal = _trapezoid(weight, log_Mh_grid)
    if n_gal <= 0:
        warnings.warn(
            f"Zero galaxy density in bias calc for [{log_Mstar_lo}, {log_Mstar_hi}] "
            f"at z={z}.",
            stacklevel=2,
        )
        return 0.0

    b_weighted = _trapezoid(weight * bh, log_Mh_grid)
    return b_weighted / n_gal


def galaxy_number_density(
    log_Mstar_lo: float,
    log_Mstar_hi: float,
    params: SHMRParams,
    z: float,
    log_Mh_grid: np.ndarray | None = None,
) -> float:
    """
    Predicted comoving galaxy number density in a stellar mass bin.

    n_gal = integral dMh (dn/dMh) * [N_cen(Mh) + N_sat(Mh)]

    Parameters
    ----------
    log_Mstar_lo : float
        Lower edge of stellar mass bin [log10(Msun)].
    log_Mstar_hi : float
        Upper edge of stellar mass bin [log10(Msun)].
    params : SHMRParams
        SHMR parameters.
    z : float
        Redshift.
    log_Mh_grid : array or None
        Halo mass grid [log10(Msun)]. Default: np.linspace(10, 15.5, 200).

    Returns
    -------
    n_gal : float
        Comoving galaxy number density [Mpc^{-3}].
    """
    if log_Mh_grid is None:
        log_Mh_grid = np.linspace(10.0, 15.5, 200)

    Mh_grid = 10.0**log_Mh_grid
    cosmo = cosmology.getCurrent()
    h = cosmo.H0 / 100.0
    Mh_grid_h = Mh_grid * h

    dndlnM = colossus_mf.massFunction(
        Mh_grid_h, z, mdef='200m', model='tinker08', q_out='dndlnM'
    )
    dndMh = dndlnM * h**3 / Mh_grid

    ncen = n_cen(log_Mh_grid, log_Mstar_lo, log_Mstar_hi, params, z)
    nsat = n_sat(log_Mh_grid, log_Mstar_lo, log_Mstar_hi, params, z)
    n_total = ncen + nsat

    weight = dndMh * n_total * Mh_grid * np.log(10.0)
    return _trapezoid(weight, log_Mh_grid)
