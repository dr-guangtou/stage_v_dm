"""Tests for Tier-3 empirical outer-correction helpers."""

import numpy as np

from sidm_stagev_forecast.calibration import resolve_tier3_parameters
from sidm_stagev_forecast.config import DEFAULT_COSMOLOGY
from sidm_stagev_forecast.outer_corrections import (
    apply_sidm_outer_correction,
    modify_dk14_parameters_for_sidm,
    sigma_response,
)
from sidm_stagev_forecast.outer_profiles import default_dk14_like_parameters


def test_sigma_response_is_bounded_and_monotonic() -> None:
    """f_sigma should be in [0, 1) and increase with sigma/m."""
    values = [sigma_response(value, sigma_pivot=1.0) for value in (0.0, 0.2, 0.5, 1.0, 2.0, 10.0)]
    assert values[0] == 0.0
    assert all(0.0 <= value < 1.0 for value in values)
    assert all(second >= first for first, second in zip(values[:-1], values[1:]))


def test_rt_gamma_shift_modifies_transition_parameters() -> None:
    """rt_gamma_shift should update DK14-like transition controls for SIDM."""
    halo = {"m200_msun": 3.0e14, "z": 0.4, "c200": 4.5}
    base = default_dk14_like_parameters(regime="cluster")
    parameters = resolve_tier3_parameters({"preset": "moderate", "correction_model": "rt_gamma_shift"})

    modified = modify_dk14_parameters_for_sidm(
        halo=halo,
        sidm_params={"sigma_over_m": 2.0},
        base_dk14_params=base,
        correction_model="rt_gamma_shift",
        correction_params=parameters,
        cosmology=DEFAULT_COSMOLOGY,
    )
    assert modified["r_t_fraction_r200m"] != base["r_t_fraction_r200m"]
    assert modified["gamma"] != base["gamma"]


def test_multiplicative_outer_window_is_outskirts_localized() -> None:
    """Outer-window correction should be near unity at small radii and vary outside."""
    halo = {"m200_msun": 3.0e14, "z": 0.4, "c200": 4.5}
    r_kpc = np.geomspace(10.0, 6000.0, 220)
    rho = 1.0e6 * (r_kpc / 100.0) ** (-2.0)
    parameters = resolve_tier3_parameters(
        {"preset": "none", "correction_model": "multiplicative_outer_window", "outer_window": {"A_outer": -0.2}}
    )

    corrected = apply_sidm_outer_correction(
        r_kpc=r_kpc,
        rho_tier2_msun_kpc3=rho,
        halo=halo,
        sidm_params={"sigma_over_m": 2.0},
        correction_model="multiplicative_outer_window",
        correction_params=parameters,
        cosmology=DEFAULT_COSMOLOGY,
    )

    ratio = corrected / rho
    assert np.all(np.isfinite(corrected))
    assert np.all(corrected > 0.0)
    assert np.isclose(ratio[0], 1.0, rtol=0.05)
    assert np.abs(ratio[-1] - 1.0) > 0.01
