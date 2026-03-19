"""
Halo model predictions for galaxy-galaxy lensing and clustering.

This module connects the SHMR parameterization (shmr_model.py) to
observable quantities by integrating over the halo mass function.
It provides:

1. HOD functions: N_cen, N_sat — mean halo occupation for a stellar mass bin.
2. Lensing signal: delta_sigma_nfw (single halo), delta_sigma_bin (bin-averaged).
3. Clustering: effective_bias, galaxy_number_density per stellar mass bin.

All halo model integrals use the colossus HMF (Tinker+2008) and integrate
over a log-spaced halo mass grid via np.trapezoid. Halo masses are in
Msun (NOT Msun/h) throughout.

References
----------
Moster, Naab & White (2013), MNRAS, 428, 3121 — SHMR
Leauthaud et al. (2012), ApJ, 744, 159 — HOD from SHMR, satellite prescription
Tinker et al. (2008), ApJ, 688, 709 — halo mass function
Tinker et al. (2010), ApJ, 724, 878 — halo bias
Diemer & Kravtsov (2015); Diemer (2019) — concentration-mass relation
"""

from __future__ import annotations

import warnings

import numpy as np
from colossus.cosmology import cosmology
from colossus.halo import concentration as colossus_conc
from colossus.halo import mass_defs as colossus_mdef
from colossus.halo import profile_nfw
from colossus.lss import bias as colossus_bias
from colossus.lss import mass_function as colossus_mf
from scipy.special import erf

from .config import SHMRParams
from .shmr_model import mean_log_Mstar


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
    sigma = params.sigma_logMs
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
    sigma = params.sigma_logMs
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
# Lensing signal: bin-averaged
# ---------------------------------------------------------------------------

def delta_sigma_bin(
    R_Mpc: np.ndarray,
    log_Mstar_lo: float,
    log_Mstar_hi: float,
    params: SHMRParams,
    z: float,
    log_Mh_grid: np.ndarray | None = None,
) -> np.ndarray:
    """
    Halo-model prediction for DeltaSigma(R) averaged over a stellar mass bin.

    DeltaSigma(R | bin, z) =
        integral dMh (dn/dMh) * [N_cen(Mh) * DS_NFW(R,Mh) + N_sat(Mh) * DS_NFW(R,Mh)]
        / n_gal(bin, z)

    The satellite term uses the same NFW profile centered on the host halo
    (1-halo approximation; ignores off-centering and subhalo profiles).

    Integration is done via np.trapezoid over a log-spaced halo mass grid.

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

    Returns
    -------
    ds : array, shape (N_R,)
        DeltaSigma [Msun/Mpc^2].
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
    n_gal = np.trapezoid(dndMh * n_total * Mh_grid * np.log(10.0),
                         log_Mh_grid)

    if n_gal <= 0:
        warnings.warn(
            f"Zero galaxy density in bin [{log_Mstar_lo}, {log_Mstar_hi}] at z={z}. "
            "Returning zeros for DeltaSigma.",
            stacklevel=2,
        )
        return np.zeros_like(R_Mpc)

    # Compute DeltaSigma for each halo mass
    # ds_grid shape: (N_Mh, N_R)
    R_Mpc = np.atleast_1d(R_Mpc)
    ds_grid = np.zeros((len(log_Mh_grid), len(R_Mpc)))

    for i, (lmh, mh) in enumerate(zip(log_Mh_grid, Mh_grid)):
        if n_total[i] < 1e-10:
            continue  # Skip halos with negligible occupation
        ds_grid[i, :] = delta_sigma_nfw(R_Mpc, mh, z)

    # Weighted average: integral dMh (dn/dMh) * N_total * DS_NFW / n_gal
    # Use log-space integration: integral = integral d(log Mh) * Mh * ln(10) * ...
    integrand = (dndMh * n_total * Mh_grid * np.log(10.0))[:, np.newaxis] * ds_grid
    ds_avg = np.trapezoid(integrand, log_Mh_grid, axis=0) / n_gal

    return ds_avg


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

    n_gal = np.trapezoid(weight, log_Mh_grid)
    if n_gal <= 0:
        warnings.warn(
            f"Zero galaxy density in bias calc for [{log_Mstar_lo}, {log_Mstar_hi}] "
            f"at z={z}.",
            stacklevel=2,
        )
        return 0.0

    b_weighted = np.trapezoid(weight * bh, log_Mh_grid)
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
    return np.trapezoid(weight, log_Mh_grid)
