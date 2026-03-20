"""
Parameterized stellar-halo mass relation (SHMR) with redshift evolution.

Implements the Moster+2013 double power-law form for M*/Mh as a function
of halo mass, with linear z/(1+z) evolution of all four shape parameters.
This module provides the forward SHMR, its inverse (via root-finding),
and the log-normal scatter distribution P(log M* | Mh, z).

These functions are the foundation for the HOD-based halo model predictions
in halo_model.py. The Fisher matrix varies these parameters to forecast
survey constraining power.

References
----------
Moster, Naab & White (2013), MNRAS, 428, 3121 (hereafter M13).
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from .config import SHMRParams


def shmr_params_at_z(
    params: SHMRParams, z: float
) -> tuple[float, float, float, float]:
    """
    Evaluate the redshift-dependent SHMR parameters at redshift z.

    Following M13 Eq. 11-14, each parameter evolves linearly in z/(1+z):
        log10(M1(z)) = log_M1_0 + nu_M1 * z/(1+z)
        N(z)         = N_0      + nu_N  * z/(1+z)
        beta(z)      = beta_0   + nu_beta * z/(1+z)
        gamma(z)     = gamma_0  + nu_gamma * z/(1+z)

    Parameters
    ----------
    params : SHMRParams
        Fiducial SHMR parameters (z=0 values + evolution coefficients).
    z : float
        Redshift. Must be >= 0.

    Returns
    -------
    log_M1 : float
        Characteristic halo mass [log10(Msun)] at redshift z.
    N : float
        Peak M*/Mh ratio [dimensionless] at redshift z.
    beta : float
        Low-mass slope [dimensionless] at redshift z.
    gamma : float
        High-mass slope [dimensionless] at redshift z.

    Raises
    ------
    ValueError
        If z < 0.
    """
    if z < 0:
        raise ValueError(f"Redshift must be non-negative, got z={z}")

    # M13 Eq. 11-14: evolution factor
    zfac = z / (1.0 + z)

    log_M1 = params.log_M1_0 + params.nu_M1 * zfac
    N = params.N_0 + params.nu_N * zfac
    beta = params.beta_0 + params.nu_beta * zfac
    gamma = params.gamma_0 + params.nu_gamma * zfac

    return log_M1, N, beta, gamma


def mean_log_Mstar(
    log_Mh: np.ndarray | float, params: SHMRParams, z: float
) -> np.ndarray | float:
    """
    Mean log10(M*/Msun) as a function of halo mass and redshift.

    Implements M13 Eq. 2:
        M*/Mh = 2N * [(Mh/M1)^(-beta) + (Mh/M1)^gamma]^(-1)
    then returns log10(M*) = log10(Mh) + log10(M*/Mh).

    Parameters
    ----------
    log_Mh : array-like or float
        log10(Mh/Msun). Values must be positive (i.e., Mh > 1 Msun).
    params : SHMRParams
        SHMR parameters.
    z : float
        Redshift.

    Returns
    -------
    log_Mstar : same shape as log_Mh
        log10(M*/Msun).
    """
    log_Mh = np.asarray(log_Mh, dtype=float)
    log_M1, N, beta, gamma = shmr_params_at_z(params, z)

    # M13 Eq. 2: x = Mh / M1
    x = 10.0 ** (log_Mh - log_M1)

    # M*/Mh = 2N / (x^(-beta) + x^gamma)
    ratio = 2.0 * N / (x ** (-beta) + x**gamma)

    result = log_Mh + np.log10(ratio)

    # Return scalar if input was scalar
    return float(result) if result.ndim == 0 else result


def mean_log_Mh(log_Mstar: float, params: SHMRParams, z: float) -> float:
    """
    Inverse SHMR: mean log10(Mh/Msun) at given log10(M*/Msun) and redshift.

    Uses scipy.optimize.brentq to invert the forward relation on the
    interval [log_Mh_min=9, log_Mh_max=16]. The SHMR is monotonically
    increasing (d(log M*)/d(log Mh) > 0 for all Mh), so the root is unique.

    Parameters
    ----------
    log_Mstar : float
        log10(M*/Msun). Typical range: 7-12.
    params : SHMRParams
        SHMR parameters.
    z : float
        Redshift.

    Returns
    -------
    log_Mh : float
        log10(Mh/Msun).

    Raises
    ------
    ValueError
        If no root is found in [9, 16].
    """
    def _residual(log_Mh_trial: float) -> float:
        return mean_log_Mstar(log_Mh_trial, params, z) - log_Mstar

    # Search over a wide halo mass range [log10(Msun)]
    log_Mh_lo = 9.0
    log_Mh_hi = 16.0

    # Check that the target M* is bracketed
    f_lo = _residual(log_Mh_lo)
    f_hi = _residual(log_Mh_hi)
    if f_lo * f_hi > 0:
        raise ValueError(
            f"Cannot invert SHMR: log_Mstar={log_Mstar:.2f} at z={z:.2f} "
            f"is outside the range spanned by log_Mh=[{log_Mh_lo}, {log_Mh_hi}]. "
            f"mean_log_Mstar at bounds: [{f_lo + log_Mstar:.2f}, {f_hi + log_Mstar:.2f}]"
        )

    return brentq(_residual, log_Mh_lo, log_Mh_hi, xtol=1e-6)


def scatter_at_Mh(
    log_Mh: np.ndarray | float,
    params: SHMRParams,
) -> np.ndarray | float:
    """
    SHMR scatter sigma(log M* | Mh) as a function of halo mass.

    If params.use_mass_dependent_scatter is False, returns the constant
    params.sigma_logMs for all halo masses.

    If True, uses the Cao & Tinker (2020) parameterization:
        sigma(Mh) = sigma_high + sigma_rise * [1 - tanh((log_Mh - log_Mh_break) / delta)]

    This gives:
    - At high Mh (log_Mh >> break): sigma -> sigma_high (~0.18 dex)
    - At low  Mh (log_Mh << break): sigma -> sigma_high + 2*sigma_rise (~0.38 dex)

    Parameters
    ----------
    log_Mh : array-like or float
        log10(Mh/Msun).
    params : SHMRParams
        SHMR parameters (includes scatter configuration).

    Returns
    -------
    sigma : same shape as log_Mh
        Scatter in log M* at fixed Mh [dex].

    References
    ----------
    Cao & Tinker (2020), MNRAS, 498, 5080 — mass-dependent SHMR scatter.
    """
    if not params.use_mass_dependent_scatter:
        log_Mh = np.asarray(log_Mh, dtype=float)
        result = np.full_like(log_Mh, params.sigma_logMs)
        return float(result) if result.ndim == 0 else result

    log_Mh = np.asarray(log_Mh, dtype=float)

    # Cao & Tinker (2020) tanh transition model
    x = (log_Mh - params.scatter_log_Mh_break) / params.scatter_delta
    sigma = params.scatter_sigma_high + params.scatter_sigma_rise * (1.0 - np.tanh(x))

    return float(sigma) if sigma.ndim == 0 else sigma


def phi_Mstar_given_Mh(
    log_Mstar: np.ndarray | float,
    log_Mh: float,
    params: SHMRParams,
    z: float,
) -> np.ndarray | float:
    """
    P(log M* | Mh, z): log-normal distribution of stellar mass at fixed Mh.

    P(log M* | Mh, z) = (1 / (sqrt(2 pi) sigma)) *
                         exp[-(log M* - <log M*>)^2 / (2 sigma^2)]

    where <log M*> = mean_log_Mstar(log_Mh, params, z) and
    sigma = scatter_at_Mh(log_Mh, params).

    Parameters
    ----------
    log_Mstar : array-like or float
        log10(M*/Msun) values at which to evaluate the PDF.
    log_Mh : float
        log10(Mh/Msun) of the host halo.
    params : SHMRParams
        SHMR parameters.
    z : float
        Redshift.

    Returns
    -------
    pdf : same shape as log_Mstar
        Probability density [per dex].

    Notes
    -----
    The scatter can be constant (params.sigma_logMs) or mass-dependent
    following Cao & Tinker (2020). See scatter_at_Mh() for details.
    Redshift evolution of scatter is assumed weak (Behroozi+2019).
    """
    log_Mstar = np.asarray(log_Mstar, dtype=float)
    mu = mean_log_Mstar(log_Mh, params, z)
    sigma = scatter_at_Mh(log_Mh, params)

    if np.any(np.asarray(sigma) <= 0):
        raise ValueError(f"Scatter must be positive, got sigma={sigma}")

    # Log-normal PDF in log-space (i.e., Gaussian in log M*)
    pdf = (
        1.0
        / (np.sqrt(2.0 * np.pi) * sigma)
        * np.exp(-0.5 * ((log_Mstar - mu) / sigma) ** 2)
    )

    return float(pdf) if pdf.ndim == 0 else pdf
