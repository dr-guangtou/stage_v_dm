"""
Analytic covariance matrices for galaxy-galaxy lensing and clustering.

Provides diagonal covariance for DeltaSigma(R) from shape noise, and
variance on (b_eff, n_gal) from cosmic variance + shot noise. These
feed into the Fisher matrix computation in fisher.py.

The lensing covariance requires Sigma_crit, which depends on the
lens-source geometry. We compute an effective Sigma_crit by averaging
over a Smail-type source redshift distribution.

Units
-----
- DeltaSigma covariance: (Msun/Mpc^2)^2
- Sigma_crit: Msun/Mpc^2
- Survey volume: comoving Mpc^3
- Source density: arcmin^{-2} (converted to Mpc^{-2} internally)

References
----------
Oguri & Takada (2011), PRD, 83, 023008 — Fisher forecast for lensing.
Mandelbaum et al. (2013), MNRAS, 432, 1544 — shape noise formalism.
"""

from __future__ import annotations

import numpy as np
from astropy import constants as const
from astropy import units as u
from colossus.cosmology import cosmology
from scipy.integrate import quad

from .config import LensingConfig


# ---------------------------------------------------------------------------
# Sigma_crit and source redshift distribution
# ---------------------------------------------------------------------------

def _smail_nz(z: float, z0: float) -> float:
    """
    Smail-type source redshift distribution (unnormalized).

    n(z) ~ z^2 exp(-(z/z0)^1.5)

    Parameters
    ----------
    z : float
        Redshift.
    z0 : float
        Characteristic redshift, related to median by z0 = z_median / 1.41.

    Returns
    -------
    nz : float
        Unnormalized probability density.
    """
    if z <= 0:
        return 0.0
    return z**2 * np.exp(-(z / z0) ** 1.5)


def sigma_crit(z_lens: float, z_source: float) -> float:
    """
    Critical surface mass density for a single lens-source pair.

    Sigma_crit = c^2 / (4 pi G) * D_s / (D_l * D_ls)

    Parameters
    ----------
    z_lens : float
        Lens redshift.
    z_source : float
        Source redshift. Must be > z_lens.

    Returns
    -------
    sigma_crit : float
        Critical surface mass density [Msun/Mpc^2].

    Notes
    -----
    Uses astropy constants for c and G, and colossus for angular
    diameter distances (returned in physical Mpc).
    """
    if z_source <= z_lens:
        return np.inf  # No lensing if source is behind or at the lens

    cosmo = cosmology.getCurrent()
    h = cosmo.H0 / 100.0

    # colossus angularDiameterDistance returns D_A in Mpc/h (positive)
    D_l = cosmo.angularDiameterDistance(z_lens) / h   # [physical Mpc]
    D_s = cosmo.angularDiameterDistance(z_source) / h  # [physical Mpc]

    # D_ls: angular diameter distance from lens to source
    # In a flat universe: D_ls = (chi_s - chi_l) / (1+z_s)
    # colossus comovingDistance uses a lookback sign convention (negative),
    # so we take absolute values.
    chi_l = abs(cosmo.comovingDistance(z_lens)) / h    # [comoving Mpc]
    chi_s = abs(cosmo.comovingDistance(z_source)) / h  # [comoving Mpc]
    D_ls = (chi_s - chi_l) / (1.0 + z_source)         # [physical Mpc]

    # Sigma_crit = c^2 / (4 pi G) * D_s / (D_l * D_ls)
    # Compute the prefactor c^2/(4 pi G) in Msun/Mpc units
    c_val = const.c.to(u.Mpc / u.s).value        # speed of light [Mpc/s]
    G_val = const.G.to(u.Mpc**3 / (u.Msun * u.s**2)).value  # G [Mpc^3/(Msun s^2)]

    prefactor = c_val**2 / (4.0 * np.pi * G_val)  # [Msun/Mpc]

    return prefactor * D_s / (D_l * D_ls)  # [Msun/Mpc^2]


def sigma_crit_effective(
    z_lens: float,
    lensing_config: LensingConfig,
) -> float:
    """
    Effective critical surface mass density averaged over source n(z).

    Sigma_crit_eff = <Sigma_crit^{-1}>^{-1}

    where <Sigma_crit^{-1}> = integral dz_s n(z_s) Sigma_crit^{-1}(z_l, z_s)
    for z_s > z_l, normalized by integral dz_s n(z_s) for z_s > z_l.

    Parameters
    ----------
    z_lens : float
        Lens redshift.
    lensing_config : LensingConfig
        Lensing survey parameters (source density, shape noise, z_median).

    Returns
    -------
    sigma_crit_eff : float
        Effective Sigma_crit [Msun/Mpc^2].
    """
    z0 = lensing_config.z_source_median / 1.41  # Smail distribution parameter

    # Integrate <Sigma_crit^{-1}> over sources behind the lens
    def integrand_inv_sc(z_s):
        nz = _smail_nz(z_s, z0)
        if nz < 1e-30:
            return 0.0
        sc = sigma_crit(z_lens, z_s)
        return nz / sc

    def integrand_nz(z_s):
        return _smail_nz(z_s, z0)

    # Integration limits: from just above z_lens to z_max
    # z_max = 4.0 is sufficient for Smail distribution with z_median ~ 1
    z_max_int = max(4.0, 3.0 * lensing_config.z_source_median)

    inv_sc_avg, _ = quad(integrand_inv_sc, z_lens + 0.01, z_max_int,
                         limit=100)
    nz_norm, _ = quad(integrand_nz, z_lens + 0.01, z_max_int,
                      limit=100)

    if nz_norm < 1e-30 or inv_sc_avg < 1e-30:
        # No sources behind the lens — Sigma_crit is effectively infinite
        return np.inf

    # <Sigma_crit^{-1}> = inv_sc_avg / nz_norm
    # Sigma_crit_eff = 1 / <Sigma_crit^{-1}>
    return nz_norm / inv_sc_avg  # [Msun/Mpc^2]


def n_source_effective(
    z_lens: float,
    lensing_config: LensingConfig,
) -> float:
    """
    Effective source density behind the lens redshift.

    n_eff(z_l) = n_total * fraction_behind(z_l)

    where fraction_behind = integral_{z_l}^inf n(z_s) dz_s / integral_0^inf n(z_s) dz_s.

    Parameters
    ----------
    z_lens : float
        Lens redshift.
    lensing_config : LensingConfig
        Lensing survey parameters.

    Returns
    -------
    n_eff : float
        Effective source density [arcmin^{-2}].
    """
    z0 = lensing_config.z_source_median / 1.41
    z_max_int = max(4.0, 3.0 * lensing_config.z_source_median)

    total, _ = quad(lambda z: _smail_nz(z, z0), 0.01, z_max_int, limit=100)
    behind, _ = quad(lambda z: _smail_nz(z, z0), z_lens, z_max_int, limit=100)

    if total < 1e-30:
        return 0.0

    fraction = behind / total
    return lensing_config.n_source_per_arcmin2 * fraction


# ---------------------------------------------------------------------------
# Lensing covariance
# ---------------------------------------------------------------------------

def lensing_covariance(
    R_bins_Mpc: np.ndarray,
    z_lens: float,
    N_lens: float,
    survey_area_deg2: float,
    lensing_config: LensingConfig,
) -> np.ndarray:
    """
    Diagonal covariance for DeltaSigma in radial bins (shape noise only).

    sigma^2(DS, R_i) = sigma_e^2 * Sigma_crit_eff^2
                       / (N_lens * n_source_eff * A_annulus_i)

    where A_annulus_i is in physical Mpc^2 per lens (the projected area
    of the annulus on the sky at the lens redshift), and n_source_eff is
    converted from arcmin^{-2} to Mpc^{-2} using the angular diameter
    distance.

    Parameters
    ----------
    R_bins_Mpc : array, shape (N_bins+1,)
        Radial bin EDGES [physical Mpc].
    z_lens : float
        Lens redshift.
    N_lens : float
        Number of lens galaxies in this (z, M*) bin.
    survey_area_deg2 : float
        Survey area [deg^2].
    lensing_config : LensingConfig
        Lensing survey parameters.

    Returns
    -------
    var_ds : array, shape (N_bins,)
        Variance sigma^2(DeltaSigma) at each radial bin [Msun/Mpc^2]^2.
    """
    if N_lens <= 0:
        return np.full(len(R_bins_Mpc) - 1, np.inf)

    sc_eff = sigma_crit_effective(z_lens, lensing_config)
    n_src_arcmin2 = n_source_effective(z_lens, lensing_config)

    if n_src_arcmin2 <= 0 or np.isinf(sc_eff):
        return np.full(len(R_bins_Mpc) - 1, np.inf)

    # Convert source density from arcmin^{-2} to physical Mpc^{-2}
    # at the lens redshift.
    # 1 arcmin = D_A * (1 arcmin in radians)
    # 1 arcmin in radians = pi / (180 * 60)
    cosmo = cosmology.getCurrent()
    h = cosmo.H0 / 100.0
    D_A = cosmo.angularDiameterDistance(z_lens) / h  # [physical Mpc]
    arcmin_to_Mpc = D_A * (np.pi / (180.0 * 60.0))  # [Mpc per arcmin]

    # n_source in Mpc^{-2} = n_source_arcmin2 / arcmin_to_Mpc^2
    n_src_Mpc2 = n_src_arcmin2 / arcmin_to_Mpc**2  # [Mpc^{-2}]

    # Shape noise variance per radial bin
    sigma_e = lensing_config.sigma_e
    n_bins = len(R_bins_Mpc) - 1
    var_ds = np.zeros(n_bins)

    for i in range(n_bins):
        R_lo = R_bins_Mpc[i]
        R_hi = R_bins_Mpc[i + 1]
        # Annulus area in physical Mpc^2
        A_annulus = np.pi * (R_hi**2 - R_lo**2)

        # Total number of source galaxies in this annulus per lens
        # N_source = n_src_Mpc2 * A_annulus
        # Total source-lens pairs = N_lens * N_source
        N_pairs = N_lens * n_src_Mpc2 * A_annulus

        if N_pairs <= 0:
            var_ds[i] = np.inf
            continue

        # sigma^2(DS) = sigma_e^2 * Sigma_crit_eff^2 / N_pairs
        # The factor of 2 in sigma_e^2 accounts for both ellipticity
        # components (sigma_e is per component; total shape noise
        # variance is sigma_e^2 per component, but DS uses the
        # tangential component only, so no factor of 2 needed here).
        var_ds[i] = sigma_e**2 * sc_eff**2 / N_pairs

    return var_ds


# ---------------------------------------------------------------------------
# Clustering covariance
# ---------------------------------------------------------------------------

def survey_volume(
    z_lo: float,
    z_hi: float,
    area_deg2: float,
) -> float:
    """
    Comoving survey volume in a redshift shell.

    V = A_sr * integral_{z_lo}^{z_hi} dV/dz dz

    where dV/dz = c/H(z) * chi^2(z) * A_sr (but we use the fraction
    of the sky times the total shell volume).

    Parameters
    ----------
    z_lo : float
        Lower redshift bound.
    z_hi : float
        Upper redshift bound.
    area_deg2 : float
        Survey area [deg^2].

    Returns
    -------
    volume : float
        Comoving volume [Mpc^3].
    """
    cosmo = cosmology.getCurrent()
    h = cosmo.H0 / 100.0

    # Sky fraction
    area_sr = area_deg2 * (np.pi / 180.0)**2  # [sr]
    f_sky = area_sr / (4.0 * np.pi)

    # Comoving volume of the full shell between z_lo and z_hi
    # colossus comovingDistance uses lookback convention (negative); take abs
    chi_lo = abs(cosmo.comovingDistance(z_lo)) / h  # [comoving Mpc]
    chi_hi = abs(cosmo.comovingDistance(z_hi)) / h  # [comoving Mpc]

    # Volume of spherical shell
    V_shell = (4.0 / 3.0) * np.pi * (chi_hi**3 - chi_lo**3)  # [Mpc^3]

    return f_sky * V_shell


def smf_covariance(
    n_gal_model: float,
    b_eff: float,
    survey_volume_Mpc3: float,
    z: float,
) -> float:
    """
    Variance on the galaxy number density from the stellar mass function.

    Combines Poisson shot noise and cosmic variance contributions:

        var(n_gal) = n_gal / V_survey + b_eff^2 * sigma_m^2 * n_gal^2

    where sigma_m^2 is the matter RMS fluctuation on the survey scale,
    computed from colossus for the effective survey radius at redshift z.

    Parameters
    ----------
    n_gal_model : float
        Model-predicted comoving galaxy number density [Mpc^{-3}].
    b_eff : float
        Effective halo bias [dimensionless].
    survey_volume_Mpc3 : float
        Comoving survey volume in this z-bin [Mpc^3].
    z : float
        Central redshift of the bin.

    Returns
    -------
    var_n : float
        Variance on n_gal [Mpc^{-6}].

    Notes
    -----
    The cosmic variance term accounts for the fact that the survey
    volume is a biased tracer of the underlying matter field. The
    effective survey radius is computed as the radius of a sphere
    with the same volume, and sigma_m is the matter RMS fluctuation
    smoothed on that scale.

    References
    ----------
    Moster, Somerville & Newman (2011), ApJ, 731, 113 — Eq. 13.
    """
    if survey_volume_Mpc3 <= 0 or n_gal_model <= 0:
        return np.inf

    cosmo = cosmology.getCurrent()
    h = cosmo.H0 / 100.0

    # Effective survey radius: sphere with same volume [Mpc]
    R_survey = (3.0 * survey_volume_Mpc3 / (4.0 * np.pi)) ** (1.0 / 3.0)
    # colossus sigma() expects R in Mpc/h
    R_survey_h = R_survey * h  # [Mpc/h]
    sigma_m = cosmo.sigma(R_survey_h, z)  # matter RMS fluctuation

    # Poisson term: var = n_gal / V [Mpc^{-6}]
    var_poisson = n_gal_model / survey_volume_Mpc3

    # Cosmic variance term: var = b_eff^2 * sigma_m^2 * n_gal^2 [Mpc^{-6}]
    var_cosmic = b_eff**2 * sigma_m**2 * n_gal_model**2

    return var_poisson + var_cosmic


def clustering_covariance(
    n_gal_model: float,
    b_eff: float,
    survey_volume_Mpc3: float,
    z: float,
) -> tuple[float, float]:
    """
    Variance on (b_eff, n_gal) from galaxy clustering.

    For the effective bias:
        sigma^2(b_eff) = b_eff^2 / (n_gal * V_eff)

    where V_eff = V_survey * [nP*b^2 / (1 + nP*b^2)] with P evaluated
    at k ~ 0.1 h/Mpc (the scale carrying most of the bias information).

    For the number density (Poisson counting):
        sigma^2(n_gal) = n_gal_model / V_survey

    Parameters
    ----------
    n_gal_model : float
        Model-predicted comoving galaxy number density [Mpc^{-3}].
        This is the tracer density from the halo model, NOT N_survey/V.
    b_eff : float
        Effective halo bias.
    survey_volume_Mpc3 : float
        Comoving survey volume in this z-bin [Mpc^3].
    z : float
        Central redshift of the bin.

    Returns
    -------
    var_b : float
        Variance on b_eff.
    var_n : float
        Variance on n_gal [Mpc^{-6}].

    Notes
    -----
    The variances depend on the MODEL-predicted density (which sets the
    Poisson rate and the tracer shot noise in V_eff) and the survey
    volume. They do NOT depend on N_gal_total — that parameter affects
    only the lensing covariance (through N_lens-source pairs).
    """
    if survey_volume_Mpc3 <= 0 or n_gal_model <= 0:
        return np.inf, np.inf

    # Linear power spectrum at k ~ 0.1 h/Mpc for V_eff calculation
    cosmo = cosmology.getCurrent()
    h = cosmo.H0 / 100.0
    k_eff = 0.1  # [h/Mpc]
    P_lin = cosmo.matterPowerSpectrum(k_eff, z)  # [(Mpc/h)^3]
    P_phys = P_lin / h**3  # [Mpc^3]

    # Galaxy power spectrum at k_eff
    P_gg = b_eff**2 * P_phys

    # Effective volume (FKP weighting)
    nP = n_gal_model * P_gg
    V_eff = survey_volume_Mpc3 * (nP / (1.0 + nP))**2 if nP > 0 else 0.0

    if V_eff <= 0:
        var_b = np.inf
    else:
        # sigma(b_eff) / b_eff = 1 / sqrt(n_gal * V_eff)
        var_b = b_eff**2 / (n_gal_model * V_eff)

    # Poisson variance on n_gal: var(n_hat) = n_model / V
    # This comes from var(N) = N_model = n_model * V, then
    # n_hat = N/V, so var(n_hat) = var(N)/V^2 = n_model/V.
    var_n = n_gal_model / survey_volume_Mpc3

    return var_b, var_n


# ---------------------------------------------------------------------------
# Nuisance parameter derivatives
# ---------------------------------------------------------------------------

def d_ln_inv_sigma_crit_d_dz(
    z_lens: float,
    lensing_config: LensingConfig,
    dz_step: float = 0.005,
) -> float:
    """
    Derivative of ln(Sigma_crit_eff^{-1}) w.r.t. source photo-z bias.

    A positive dz_source shifts the effective source redshift distribution
    to higher z, increasing lensing efficiency (decreasing Sigma_crit_eff).
    This derivative quantifies how the lensing signal changes when the
    source redshift calibration is biased.

    Parameters
    ----------
    z_lens : float
        Lens redshift.
    lensing_config : LensingConfig
        Lensing survey configuration.
    dz_step : float
        Step size for central finite difference on z_source_median.

    Returns
    -------
    dln_inv_sc_ddz : float
        d(ln Sigma_crit_eff^{-1}) / d(dz_source) [dimensionless].

    Notes
    -----
    Uses central differences: perturbs z_source_median by +/- dz_step
    and computes the log-derivative of the inverse Sigma_crit_eff.
    """
    from dataclasses import replace

    config_plus = replace(
        lensing_config,
        z_source_median=lensing_config.z_source_median + dz_step,
    )
    config_minus = replace(
        lensing_config,
        z_source_median=lensing_config.z_source_median - dz_step,
    )

    sc_plus = sigma_crit_effective(z_lens, config_plus)
    sc_minus = sigma_crit_effective(z_lens, config_minus)

    if np.isinf(sc_plus) or np.isinf(sc_minus) or sc_plus <= 0 or sc_minus <= 0:
        return 0.0

    # d(ln(1/Sc))/ddz = -d(ln Sc)/ddz
    return -(np.log(sc_plus) - np.log(sc_minus)) / (2.0 * dz_step)
