"""Mass-conservation and profile-sanity tests."""

import numpy as np
import pytest

from sidm_stagev_forecast.config import DEFAULT_COSMOLOGY
from sidm_stagev_forecast.cosmology import rdelta
from sidm_stagev_forecast.profiles import nfw_profile_bundle, sidm_profile_from_parametric_model


def test_nfw_m200c_mass_consistency_at_r200() -> None:
    """NFW enclosed mass at r200 should match input M200 within tight tolerance."""
    m200_msun = 1.0e14
    c200 = 4.0
    z = 0.3

    r200_kpc = rdelta(m200_msun, z, DEFAULT_COSMOLOGY, definition="200c")
    r_grid_kpc = np.geomspace(1.0e-4 * r200_kpc, r200_kpc, 1500)

    profile = nfw_profile_bundle(r_grid_kpc, m200_msun, c200, z, DEFAULT_COSMOLOGY)
    measured_mass_msun = profile["m_enclosed_msun"][-1]

    relative_error = np.abs(measured_mass_msun / m200_msun - 1.0)
    assert relative_error < 1.0e-6


def test_nfw_profile_sanity_monotonic_enclosed_mass() -> None:
    """NFW profile bundle should produce physical and monotonic outputs."""
    r_grid_kpc = np.geomspace(1.0e-2, 5.0e3, 500)
    profile = nfw_profile_bundle(r_grid_kpc, 1.0e10, 15.0, 0.3, DEFAULT_COSMOLOGY)

    rho = profile["rho_msun_kpc3"]
    m_enclosed = profile["m_enclosed_msun"]
    vcirc = profile["vcirc_km_s"]

    assert np.all(np.isfinite(rho))
    assert np.all(np.isfinite(m_enclosed))
    assert np.all(np.isfinite(vcirc))
    assert np.all(rho > 0.0)
    assert np.all(np.diff(m_enclosed) >= 0.0)
    assert np.all(vcirc >= 0.0)


def test_sidm_wrapper_profile_sanity() -> None:
    """SIDM wrapper output should be finite and physically ordered for a benchmark case."""
    r_grid_kpc = np.geomspace(1.0e-2, 3.0e2, 220)

    try:
        profile = sidm_profile_from_parametric_model(
            r_kpc=r_grid_kpc,
            m200_msun=1.0e10,
            c200=15.0,
            z=0.3,
            sigma_over_m=1.0,
            model_options=None,
        )
    except (RuntimeError, FileNotFoundError) as error:
        pytest.skip(f"parametricSIDM backend unavailable in this environment: {error}")

    rho = profile["rho_msun_kpc3"]
    m_enclosed = profile["m_enclosed_msun"]
    vcirc = profile["vcirc_km_s"]

    assert np.all(np.isfinite(rho))
    assert np.all(np.isfinite(m_enclosed))
    assert np.all(np.isfinite(vcirc))
    assert np.all(rho > 0.0)
    assert np.all(np.diff(m_enclosed) >= -1.0e-8)
    assert np.all(vcirc >= 0.0)
