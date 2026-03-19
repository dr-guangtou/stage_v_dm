"""Tier-3 empirical calibration presets for SIDM outskirts corrections.

This module provides simulation-informed but approximate parameter presets.
It does not implement first-principles SIDM splashback dynamics.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any


TIER3_PRESETS: dict[str, dict[str, Any]] = {
    "none": {
        "correction_model": "rt_gamma_shift",
        "sigma_pivot": 1.0,
        "rt_shift": {"A_rt": 0.0, "use_concentration_dependence": False, "A_c": 0.0},
        "gamma_shift": {"A_gamma": 0.0},
        "beta_shift": {"A_beta": 0.0},
        "outer_window": {"A_outer": 0.0, "r_turn_mode": "fraction_r200m", "r_turn_value": 1.0, "width_dex": 0.2},
    },
    "mild": {
        "correction_model": "rt_gamma_shift",
        "sigma_pivot": 1.0,
        "rt_shift": {"A_rt": -0.04, "use_concentration_dependence": False, "A_c": 0.0},
        "gamma_shift": {"A_gamma": -0.10},
        "beta_shift": {"A_beta": 0.0},
        "outer_window": {"A_outer": -0.05, "r_turn_mode": "fraction_r200m", "r_turn_value": 1.0, "width_dex": 0.2},
    },
    "moderate": {
        "correction_model": "rt_gamma_shift",
        "sigma_pivot": 1.0,
        "rt_shift": {"A_rt": -0.08, "use_concentration_dependence": True, "A_c": 0.35},
        "gamma_shift": {"A_gamma": -0.20},
        "beta_shift": {"A_beta": -0.05},
        "outer_window": {"A_outer": -0.10, "r_turn_mode": "fraction_r200m", "r_turn_value": 1.0, "width_dex": 0.2},
    },
    "strong": {
        "correction_model": "rt_gamma_shift",
        "sigma_pivot": 1.0,
        "rt_shift": {"A_rt": -0.15, "use_concentration_dependence": True, "A_c": 0.60},
        "gamma_shift": {"A_gamma": -0.35},
        "beta_shift": {"A_beta": -0.10},
        "outer_window": {"A_outer": -0.18, "r_turn_mode": "fraction_r200m", "r_turn_value": 1.0, "width_dex": 0.22},
    },
    "shallower_splashback": {
        "correction_model": "gamma_shift",
        "sigma_pivot": 1.0,
        "rt_shift": {"A_rt": 0.0, "use_concentration_dependence": False, "A_c": 0.0},
        "gamma_shift": {"A_gamma": -0.25},
        "beta_shift": {"A_beta": -0.08},
        "outer_window": {"A_outer": 0.0, "r_turn_mode": "fraction_r200m", "r_turn_value": 1.0, "width_dex": 0.2},
    },
    "high_c_enhanced": {
        "correction_model": "rt_gamma_shift",
        "sigma_pivot": 1.0,
        "rt_shift": {"A_rt": -0.06, "use_concentration_dependence": True, "A_c": 0.80},
        "gamma_shift": {"A_gamma": -0.15},
        "beta_shift": {"A_beta": 0.0},
        "outer_window": {"A_outer": 0.0, "r_turn_mode": "fraction_r200m", "r_turn_value": 1.0, "width_dex": 0.2},
    },
    "small_inward_shift": {
        "correction_model": "rt_shift",
        "sigma_pivot": 1.0,
        "rt_shift": {"A_rt": -0.05, "use_concentration_dependence": False, "A_c": 0.0},
        "gamma_shift": {"A_gamma": 0.0},
        "beta_shift": {"A_beta": 0.0},
        "outer_window": {"A_outer": 0.0, "r_turn_mode": "fraction_r200m", "r_turn_value": 1.0, "width_dex": 0.2},
    },
    "small_outward_shift": {
        "correction_model": "rt_shift",
        "sigma_pivot": 1.0,
        "rt_shift": {"A_rt": 0.05, "use_concentration_dependence": False, "A_c": 0.0},
        "gamma_shift": {"A_gamma": 0.0},
        "beta_shift": {"A_beta": 0.0},
        "outer_window": {"A_outer": 0.0, "r_turn_mode": "fraction_r200m", "r_turn_value": 1.0, "width_dex": 0.2},
    },
}


def resolve_tier3_parameters(tier3_config: dict[str, Any] | None) -> dict[str, Any]:
    """Resolve Tier-3 runtime parameters from config and preset selection."""
    config = {} if tier3_config is None else dict(tier3_config)
    preset_name = str(config.get("preset", "none")).lower()
    if preset_name not in TIER3_PRESETS:
        raise ValueError(f"Unknown Tier-3 preset: {preset_name}")

    parameters = deepcopy(TIER3_PRESETS[preset_name])
    for key in ("correction_model", "sigma_pivot"):
        if key in config:
            parameters[key] = config[key]

    for block in ("rt_shift", "gamma_shift", "beta_shift", "outer_window"):
        if isinstance(config.get(block), dict):
            for key, value in dict(config[block]).items():
                parameters[block][key] = value
    return parameters
