"""Profile stitching utilities for Tier-2 hybrid SIDM forecasts."""

from __future__ import annotations

from typing import Any

import numpy as np

from sidm_stagev_forecast.config import CosmologyConfig
from sidm_stagev_forecast.cosmology import rdelta


def resolve_match_radius_kpc(
    mass_msun: float,
    z: float,
    cosmo: CosmologyConfig,
    stitch_config: dict[str, Any] | None = None,
) -> float:
    """Resolve stitch match radius in kpc from config."""
    options = {} if stitch_config is None else dict(stitch_config)
    mode = str(options.get("r_match_mode", "fraction_r200m")).lower()
    value = float(options.get("r_match_value", 0.8))

    if mode == "fraction_r200m":
        return value * rdelta(mass_msun, z, cosmo, definition="200m")
    if mode == "fraction_r200c":
        return value * rdelta(mass_msun, z, cosmo, definition="200c")
    if mode == "fixed_kpc":
        if value <= 0.0:
            raise ValueError("fixed_kpc match radius must be positive.")
        return value
    raise ValueError("Unsupported r_match_mode.")


def logistic_log_radius_weight(
    r_kpc: np.ndarray,
    r_match_kpc: float,
    smooth_width_dex: float,
) -> np.ndarray:
    """Return a smooth inner-to-outer logistic weight in log-radius."""
    if r_match_kpc <= 0.0:
        raise ValueError("r_match_kpc must be positive.")
    if smooth_width_dex <= 0.0:
        raise ValueError("smooth_width_dex must be positive.")

    log_radius = np.log10(np.asarray(r_kpc, dtype=float))
    width_natural = smooth_width_dex * np.log(10.0)
    x_value = (log_radius - np.log10(r_match_kpc)) / width_natural
    return 1.0 / (1.0 + np.exp(x_value))


def stitch_inner_outer_profile(
    r_kpc: np.ndarray,
    rho_inner_sidm_msun_kpc3: np.ndarray,
    rho_outer_reference_msun_kpc3: np.ndarray,
    method: str = "logistic_logrho_blend",
    r_match_kpc: float | None = None,
    smooth_width_dex: float = 0.15,
    continuity: str = "density",
) -> np.ndarray:
    """Blend inner and outer profiles into a numerically smooth hybrid profile."""
    if continuity.lower() != "density":
        raise ValueError("Tier-2 supports continuity='density' only.")

    radius_kpc = np.asarray(r_kpc, dtype=float)
    rho_inner = np.asarray(rho_inner_sidm_msun_kpc3, dtype=float)
    rho_outer = np.asarray(rho_outer_reference_msun_kpc3, dtype=float)

    if rho_inner.shape != radius_kpc.shape or rho_outer.shape != radius_kpc.shape:
        raise ValueError("Radius and density arrays must have matching shapes.")
    if np.any(radius_kpc <= 0.0):
        raise ValueError("All radii must be positive.")
    if np.any(rho_inner <= 0.0) or np.any(rho_outer <= 0.0):
        raise ValueError("Density arrays must be strictly positive.")
    if r_match_kpc is None:
        raise ValueError("r_match_kpc must be provided.")

    method_normalized = method.lower()
    if method_normalized != "logistic_logrho_blend":
        raise ValueError("Unsupported stitch method.")

    weight_inner = logistic_log_radius_weight(
        r_kpc=radius_kpc,
        r_match_kpc=r_match_kpc,
        smooth_width_dex=smooth_width_dex,
    )

    log_rho_hybrid = weight_inner * np.log(rho_inner) + (1.0 - weight_inner) * np.log(rho_outer)
    rho_hybrid = np.exp(log_rho_hybrid)
    if np.any(~np.isfinite(rho_hybrid)):
        raise ValueError("Non-finite hybrid density encountered after stitching.")
    return rho_hybrid
