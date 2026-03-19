"""Tier-3 empirical SIDM outer-profile corrections.

These corrections are simulation-informed nuisance models for sensitivity tests.
They are not first-principles SIDM splashback solvers.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from sidm_stagev_forecast.config import CosmologyConfig
from sidm_stagev_forecast.cosmology import rdelta
from sidm_stagev_forecast.ensemble import concentration_from_mass_maccio


def sigma_response(sigma_over_m: float, sigma_pivot: float = 1.0) -> float:
    """Bounded monotonic SIDM response function f_sigma."""
    sigma = max(float(sigma_over_m), 0.0)
    pivot = max(float(sigma_pivot), 1.0e-8)
    return sigma / (pivot + sigma)


def concentration_residual_standardized(
    halo: dict[str, float],
    cosmology: CosmologyConfig,
    scatter_dex: float = 0.15,
) -> float:
    """Return standardized log-concentration residual relative to maccio relation."""
    c_expected = concentration_from_mass_maccio(halo["m200_msun"], halo["z"])
    residual = (np.log10(halo["c200"]) - np.log10(c_expected)) / max(scatter_dex, 1.0e-4)
    return float(np.clip(residual, -3.0, 3.0))


def modify_dk14_parameters_for_sidm(
    halo: dict[str, float],
    sidm_params: dict[str, float],
    base_dk14_params: dict[str, float],
    correction_model: str,
    correction_params: dict[str, Any],
    cosmology: CosmologyConfig,
) -> dict[str, float]:
    """Apply empirical SIDM-dependent parameter shifts to DK14-like controls."""
    model_name = str(correction_model).lower()
    output = dict(base_dk14_params)
    response = sigma_response(
        sidm_params.get("sigma_over_m", 0.0),
        sigma_pivot=float(correction_params.get("sigma_pivot", 1.0)),
    )

    concentration_factor = 1.0
    rt_shift = dict(correction_params.get("rt_shift", {}))
    if bool(rt_shift.get("use_concentration_dependence", False)):
        q_c = concentration_residual_standardized(halo=halo, cosmology=cosmology)
        concentration_factor = 1.0 + float(rt_shift.get("A_c", 0.0)) * q_c

    if model_name in {"rt_shift", "rt_gamma_shift"}:
        a_rt = float(rt_shift.get("A_rt", 0.0))
        output["r_t_fraction_r200m"] = max(
            0.4,
            float(base_dk14_params["r_t_fraction_r200m"]) * (1.0 + a_rt * response * concentration_factor),
        )

    if model_name in {"gamma_shift", "rt_gamma_shift"}:
        gamma_shift = dict(correction_params.get("gamma_shift", {}))
        beta_shift = dict(correction_params.get("beta_shift", {}))
        output["gamma"] = max(
            0.5,
            float(base_dk14_params["gamma"]) * (1.0 + float(gamma_shift.get("A_gamma", 0.0)) * response),
        )
        output["beta"] = max(
            0.5,
            float(base_dk14_params["beta"]) * (1.0 + float(beta_shift.get("A_beta", 0.0)) * response),
        )

    return output


def apply_sidm_outer_correction(
    r_kpc: np.ndarray,
    rho_tier2_msun_kpc3: np.ndarray,
    halo: dict[str, float],
    sidm_params: dict[str, float],
    correction_model: str,
    correction_params: dict[str, Any],
    cosmology: CosmologyConfig,
) -> np.ndarray:
    """Apply an optional multiplicative empirical outer-window correction."""
    model_name = str(correction_model).lower()
    if model_name != "multiplicative_outer_window":
        return np.asarray(rho_tier2_msun_kpc3, dtype=float)

    outer_window = dict(correction_params.get("outer_window", {}))
    amplitude = float(outer_window.get("A_outer", 0.0))
    response = sigma_response(
        sidm_params.get("sigma_over_m", 0.0),
        sigma_pivot=float(correction_params.get("sigma_pivot", 1.0)),
    )
    if np.isclose(amplitude * response, 0.0):
        return np.asarray(rho_tier2_msun_kpc3, dtype=float)

    turn_mode = str(outer_window.get("r_turn_mode", "fraction_r200m")).lower()
    turn_value = float(outer_window.get("r_turn_value", 1.0))
    if turn_mode == "fraction_r200m":
        r_turn_kpc = turn_value * rdelta(halo["m200_msun"], halo["z"], cosmology, definition="200m")
    elif turn_mode == "fraction_r200c":
        r_turn_kpc = turn_value * rdelta(halo["m200_msun"], halo["z"], cosmology, definition="200c")
    else:
        r_turn_kpc = turn_value
    r_turn_kpc = max(r_turn_kpc, 1.0e-6)

    width_dex = max(float(outer_window.get("width_dex", 0.2)), 1.0e-3)
    log_radius = np.log10(np.asarray(r_kpc, dtype=float))
    window = 1.0 / (1.0 + np.exp(-(log_radius - np.log10(r_turn_kpc)) / width_dex))

    log_correction = amplitude * response * window
    return np.asarray(rho_tier2_msun_kpc3, dtype=float) * np.exp(log_correction)
