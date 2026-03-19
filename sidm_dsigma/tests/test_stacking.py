"""Tests for Tier-1 stacking utilities."""

import numpy as np

from sidm_stagev_forecast.stacking import (
    interpolate_profile_to_common_grid,
    stack_delta_sigma_profiles,
)


def test_interpolate_profile_to_common_grid_exact_for_power_law() -> None:
    """Log-log interpolation should exactly reproduce a power-law profile."""
    r_input = np.geomspace(10.0, 1000.0, 20)
    profile_input = 5.0e4 * (r_input / 100.0) ** (-0.8)
    r_common = np.geomspace(12.0, 900.0, 25)

    interpolated = interpolate_profile_to_common_grid(r_input, profile_input, r_common)
    expected = 5.0e4 * (r_common / 100.0) ** (-0.8)

    assert np.allclose(interpolated, expected, rtol=1.0e-12, atol=0.0)


def test_stack_delta_sigma_profiles_single_halo_limit() -> None:
    """A single-halo stack should return that halo profile unchanged."""
    r_common = np.geomspace(30.0, 3000.0, 30)
    profile = 1.0e8 * (r_common / 100.0) ** (-0.6)
    projected_profiles = [
        {
            "r_projected_kpc": r_common,
            "delta_sigma_msun_kpc2": profile,
        }
    ]

    stacked = stack_delta_sigma_profiles(r_common, projected_profiles)

    assert np.allclose(stacked["delta_sigma_msun_kpc2"], profile, rtol=1.0e-12, atol=0.0)


def test_stack_delta_sigma_profiles_weighted_mean_behavior() -> None:
    """Weighted stack should reproduce the analytic two-profile weighted mean."""
    r_common = np.geomspace(30.0, 3000.0, 40)
    profile_a = 2.0e8 * (r_common / 100.0) ** (-0.7)
    profile_b = 1.2e8 * (r_common / 100.0) ** (-0.5)

    stacked = stack_delta_sigma_profiles(
        r_common_kpc=r_common,
        projected_profiles=[
            {"r_projected_kpc": r_common, "delta_sigma_msun_kpc2": profile_a},
            {"r_projected_kpc": r_common, "delta_sigma_msun_kpc2": profile_b},
        ],
        weights=np.asarray([0.25, 0.75]),
    )

    expected = 0.25 * profile_a + 0.75 * profile_b
    assert np.allclose(stacked["delta_sigma_msun_kpc2"], expected, rtol=1.0e-12, atol=0.0)
