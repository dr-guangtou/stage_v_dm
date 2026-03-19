"""Velocity-dependent SIDM helpers adapted from parametricSIDM scripts."""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad

from sidm_stagev_forecast.config import DEFAULT_COSMOLOGY
from sidm_stagev_forecast.cosmology import GRAVITATIONAL_CONSTANT_KPC_KMS


def sigma_vis_velocity_dependent(
    v_km_s: float,
    sigma0_over_m_cm2_g: float,
    w_km_s: float,
) -> float:
    """Return velocity-dependent cross section in cm^2/g."""
    if sigma0_over_m_cm2_g <= 0.0 or w_km_s <= 0.0 or v_km_s <= 0.0:
        return 0.0

    v2 = float(v_km_s) ** 2
    w2 = float(w_km_s) ** 2
    ratio = v2 / w2
    value = (
        3.0
        * float(sigma0_over_m_cm2_g)
        * float(w_km_s) ** 6
        * ((-4.0 * ratio) + 2.0 * (2.0 + ratio) * np.log1p(ratio))
        / float(v_km_s) ** 6
    )
    return float(max(value, 0.0))


def sigma_eff_maxwellian(
    veff_km_s: float,
    sigma0_over_m_cm2_g: float,
    w_km_s: float,
) -> float:
    """Return Maxwellian-averaged effective cross section in cm^2/g."""
    if sigma0_over_m_cm2_g <= 0.0 or w_km_s <= 0.0 or veff_km_s <= 0.0:
        return 0.0

    veff = float(veff_km_s)

    def integrand(ln_velocity: float) -> float:
        velocity = float(np.exp(ln_velocity))
        return (
            velocity**7
            * np.exp(-(velocity**2) / (4.0 * veff**2))
            * (2.0 / 3.0)
            * sigma_vis_velocity_dependent(velocity, sigma0_over_m_cm2_g, w_km_s)
            * velocity
        )

    value, _ = quad(integrand, np.log(0.01), np.log(15000.0))
    sigma_eff = value / (512.0 * veff**8)
    return float(max(sigma_eff, 0.0))


def nfw_vmax_km_s(rho_s_msun_kpc3: float, rs_kpc: float) -> float:
    """Return NFW Vmax in km/s from `(rho_s, rs)` in physical units."""
    if rho_s_msun_kpc3 <= 0.0 or rs_kpc <= 0.0:
        return 0.0
    return float(1.648 * np.sqrt(GRAVITATIONAL_CONSTANT_KPC_KMS * rho_s_msun_kpc3) * rs_kpc)


def formation_redshift_from_mass_fit(
    mass_msun: float,
    *,
    h: float = DEFAULT_COSMOLOGY.h,
) -> float:
    """Return formation-redshift proxy from the basic parametricSIDM fit."""
    mass_msun_h = float(mass_msun) * float(h)
    log10_mass = float(np.log10(max(mass_msun_h, 1.0)))
    z_form = -0.0064 * log10_mass**2 + 0.0237 * log10_mass + 1.8837
    return float(max(z_form, 0.0))


def elapsed_time_since_formation_gyr(
    *,
    mass_msun: float,
    z_observation: float,
    tlb_callable,
    formation_redshift: float | None = None,
) -> float:
    """Compute elapsed SIDM evolution time from formation to observation."""
    z_form = (
        formation_redshift_from_mass_fit(mass_msun)
        if formation_redshift is None
        else float(formation_redshift)
    )
    elapsed_time_gyr = float(tlb_callable(z_form) - tlb_callable(float(z_observation)))
    return float(max(elapsed_time_gyr, 0.0))
