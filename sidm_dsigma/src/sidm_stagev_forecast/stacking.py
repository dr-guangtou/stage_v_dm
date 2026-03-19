"""Stacking utilities for projected lensing profiles."""

from __future__ import annotations

import numpy as np
from scipy.interpolate import interp1d


def interpolate_profile_to_common_grid(
    r_input_kpc: np.ndarray,
    profile_input: np.ndarray,
    r_common_kpc: np.ndarray,
) -> np.ndarray:
    """Interpolate a positive profile onto a common radius grid in log-log space."""
    r_input_kpc = np.asarray(r_input_kpc, dtype=float)
    profile_input = np.asarray(profile_input, dtype=float)
    r_common_kpc = np.asarray(r_common_kpc, dtype=float)

    if not np.all(np.diff(r_input_kpc) > 0.0):
        raise ValueError("r_input_kpc must be strictly increasing.")
    if not np.all(np.diff(r_common_kpc) > 0.0):
        raise ValueError("r_common_kpc must be strictly increasing.")
    if np.any(r_input_kpc <= 0.0) or np.any(r_common_kpc <= 0.0):
        raise ValueError("All radii must be positive.")
    if np.any(profile_input <= 0.0):
        raise ValueError("Profile values must be positive for log-log interpolation.")

    log_interpolator = interp1d(
        np.log(r_input_kpc),
        np.log(profile_input),
        kind="linear",
        bounds_error=False,
        fill_value=(-np.inf, -np.inf),
        assume_sorted=True,
    )

    profile_common = np.exp(log_interpolator(np.log(r_common_kpc)))
    profile_common[~np.isfinite(profile_common)] = 0.0
    return profile_common


def weighted_stack_profiles(
    profile_matrix: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    """Return weighted mean profile for shape (n_halos, n_radii)."""
    profile_matrix = np.asarray(profile_matrix, dtype=float)
    weights = np.asarray(weights, dtype=float)

    if profile_matrix.ndim != 2:
        raise ValueError("profile_matrix must have shape (n_halos, n_radii).")
    if profile_matrix.shape[0] != weights.shape[0]:
        raise ValueError("weights length must match profile_matrix first axis.")
    if np.any(weights < 0.0):
        raise ValueError("weights must be non-negative.")

    normalization = np.sum(weights)
    if normalization <= 0.0:
        raise ValueError("Sum of weights must be positive.")

    normalized_weights = weights / normalization
    return np.sum(profile_matrix * normalized_weights[:, np.newaxis], axis=0)


def stack_delta_sigma_profiles(
    r_common_kpc: np.ndarray,
    projected_profiles: list[dict[str, np.ndarray]],
    weights: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    """Interpolate halo DeltaSigma profiles to a common grid and stack them."""
    if len(projected_profiles) == 0:
        raise ValueError("projected_profiles cannot be empty.")

    n_halos = len(projected_profiles)
    if weights is None:
        weights = np.full(n_halos, 1.0 / n_halos)
    else:
        weights = np.asarray(weights, dtype=float)

    profile_rows: list[np.ndarray] = []
    for projected_profile in projected_profiles:
        interpolated = interpolate_profile_to_common_grid(
            r_input_kpc=np.asarray(projected_profile["r_projected_kpc"]),
            profile_input=np.asarray(projected_profile["delta_sigma_msun_kpc2"]),
            r_common_kpc=r_common_kpc,
        )
        profile_rows.append(interpolated)

    profile_matrix = np.vstack(profile_rows)
    stacked = weighted_stack_profiles(profile_matrix, weights)

    return {
        "r_projected_kpc": np.asarray(r_common_kpc, dtype=float),
        "delta_sigma_msun_kpc2": stacked,
        "profile_matrix_msun_kpc2": profile_matrix,
        "weights": np.asarray(weights, dtype=float),
    }
