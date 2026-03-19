"""Tests for DK14-like outer-profile utilities."""

import numpy as np

from sidm_stagev_forecast.config import DEFAULT_COSMOLOGY
from sidm_stagev_forecast.outer_profiles import build_dk14_outer_profile
from sidm_stagev_forecast.plotting import log_slope


def test_dk14_like_profile_positive_and_finite() -> None:
    """DK14-like profile components should remain finite and positive."""
    r_kpc = np.geomspace(5.0, 8000.0, 240)
    profile = build_dk14_outer_profile(
        r_kpc=r_kpc,
        mass_msun=3.0e14,
        z=0.4,
        concentration=4.0,
        mass_def="200c",
        cosmology=DEFAULT_COSMOLOGY,
        dk14_params={"regime": "cluster"},
    )

    assert np.all(np.isfinite(profile["rho_total_msun_kpc3"]))
    assert np.all(np.isfinite(profile["rho_inner_baseline_msun_kpc3"]))
    assert np.all(np.isfinite(profile["rho_outer_term_msun_kpc3"]))
    assert np.all(profile["rho_total_msun_kpc3"] > 0.0)
    assert np.all(profile["rho_outer_term_msun_kpc3"] > 0.0)


def test_dk14_like_profile_has_outer_steepening_signature() -> None:
    """A DK14-like profile should show steepening relative to inner halo behavior."""
    r_kpc = np.geomspace(5.0, 8000.0, 320)
    profile = build_dk14_outer_profile(
        r_kpc=r_kpc,
        mass_msun=3.0e14,
        z=0.4,
        concentration=4.0,
        mass_def="200c",
        cosmology=DEFAULT_COSMOLOGY,
        dk14_params={"regime": "cluster"},
    )
    slope = log_slope(r_kpc, profile["rho_total_msun_kpc3"])

    assert np.min(slope) < -2.8
