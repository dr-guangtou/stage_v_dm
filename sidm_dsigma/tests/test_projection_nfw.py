"""Validation tests for projection module."""

import numpy as np

from sidm_stagev_forecast.config import DEFAULT_COSMOLOGY
from sidm_stagev_forecast.profiles import nfw_parameters_from_m_c, nfw_profile_from_m_c
from sidm_stagev_forecast.projection import analytic_nfw_delta_sigma, delta_sigma_of_R


def test_numeric_projection_matches_analytic_nfw_delta_sigma() -> None:
    """Numeric projection should reproduce analytic NFW DeltaSigma within tolerance."""
    m200_msun = 1.0e14
    c200 = 4.0
    z = 0.3

    r_3d_kpc = np.geomspace(0.5, 1.0e4, 500)
    r_lensing_kpc = np.geomspace(30.0, 2.0e3, 40)

    rho_msun_kpc3 = nfw_profile_from_m_c(r_3d_kpc, m200_msun, c200, z, DEFAULT_COSMOLOGY)
    numeric = delta_sigma_of_R(r_lensing_kpc, r_3d_kpc, rho_msun_kpc3, n_z=1400)

    parameters = nfw_parameters_from_m_c(m200_msun, c200, z, DEFAULT_COSMOLOGY)
    analytic = analytic_nfw_delta_sigma(
        r_lensing_kpc, parameters.rho_s_msun_kpc3, parameters.rs_kpc
    )

    ratio = numeric["delta_sigma_msun_kpc2"] / analytic["delta_sigma_msun_kpc2"]
    max_abs_fractional_error = np.max(np.abs(ratio - 1.0))

    assert max_abs_fractional_error < 0.08
