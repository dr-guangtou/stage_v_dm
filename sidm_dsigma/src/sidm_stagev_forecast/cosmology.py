"""Cosmology helpers with explicit physical units.

Conventions
-----------
- Radii are in kpc.
- Masses are in Msun.
- Densities are in Msun / kpc^3.
- Mass definition defaults to M200c.
"""

import math

from sidm_stagev_forecast.config import CosmologyConfig

GRAVITATIONAL_CONSTANT_KPC_KMS = 4.30091e-6


def e_z(z: float, cosmo: CosmologyConfig) -> float:
    """Dimensionless expansion rate E(z)."""
    return math.sqrt(cosmo.omega_m * (1.0 + z) ** 3 + cosmo.omega_lambda)


def hubble_km_s_kpc(z: float, cosmo: CosmologyConfig) -> float:
    """Hubble parameter in km/s/kpc."""
    h0_km_s_mpc = 100.0 * cosmo.h
    return (h0_km_s_mpc / 1000.0) * e_z(z, cosmo)


def rho_crit_z(z: float, cosmo: CosmologyConfig) -> float:
    """Critical density in Msun / kpc^3 at redshift z."""
    hubble_value = hubble_km_s_kpc(z, cosmo)
    return 3.0 * hubble_value**2 / (8.0 * math.pi * GRAVITATIONAL_CONSTANT_KPC_KMS)


def omega_m_z(z: float, cosmo: CosmologyConfig) -> float:
    """Return the matter density parameter Omega_m(z)."""
    ez_squared = e_z(z, cosmo) ** 2
    return cosmo.omega_m * (1.0 + z) ** 3 / ez_squared


def rho_mean_matter_z(z: float, cosmo: CosmologyConfig) -> float:
    """Mean matter density in Msun / kpc^3 at redshift z."""
    return omega_m_z(z, cosmo) * rho_crit_z(z, cosmo)


def rdelta(mass_msun: float, z: float, cosmo: CosmologyConfig, definition: str = "200c") -> float:
    """Radius for spherical-overdensity mass definitions, returned in kpc."""
    definition_normalized = definition.lower()
    if definition_normalized == "200c":
        density_ref = rho_crit_z(z, cosmo)
    elif definition_normalized == "200m":
        density_ref = rho_mean_matter_z(z, cosmo)
    else:
        raise ValueError("Supported definitions are '200c' and '200m'.")
    delta = 200.0
    return (3.0 * mass_msun / (4.0 * math.pi * delta * density_ref)) ** (1.0 / 3.0)
