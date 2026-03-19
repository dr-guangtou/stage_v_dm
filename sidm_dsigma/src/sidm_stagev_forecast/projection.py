"""Numerical projection from 3D density to weak-lensing observables."""

from __future__ import annotations

import numpy as np
from scipy.interpolate import interp1d


def _loglog_interpolator(r_kpc: np.ndarray, rho_msun_kpc3: np.ndarray) -> interp1d:
    """Return a log-log interpolator with safe zero fill outside range."""
    if np.any(r_kpc <= 0.0):
        raise ValueError("All radii must be positive for log-log interpolation.")
    if np.any(rho_msun_kpc3 <= 0.0):
        raise ValueError("All densities must be positive for log-log interpolation.")

    return interp1d(
        np.log(r_kpc),
        np.log(rho_msun_kpc3),
        kind="linear",
        bounds_error=False,
        fill_value=(-np.inf, -np.inf),
        assume_sorted=True,
    )


def sigma_of_R(
    r_projected_kpc: np.ndarray,
    r_kpc: np.ndarray,
    rho_msun_kpc3: np.ndarray,
    n_z: int = 800,
    r_max_kpc: float | None = None,
) -> np.ndarray:
    """Compute Sigma(R) in Msun/kpc^2 via stable line-of-sight integration.

    Uses substitution r = sqrt(R^2 + z^2), so
    Sigma(R) = 2 * integral_0^zmax rho(sqrt(R^2 + z^2)) dz.
    """
    r_projected_kpc = np.asarray(r_projected_kpc)
    r_kpc = np.asarray(r_kpc)
    rho_msun_kpc3 = np.asarray(rho_msun_kpc3)

    if not np.all(np.diff(r_projected_kpc) > 0.0):
        raise ValueError("Projected radii must be strictly increasing.")
    if not np.all(np.diff(r_kpc) > 0.0):
        raise ValueError("3D radii must be strictly increasing.")

    interpolation = _loglog_interpolator(r_kpc, rho_msun_kpc3)
    effective_r_max = float(r_max_kpc if r_max_kpc is not None else r_kpc[-1])

    sigma_msun_kpc2 = np.zeros_like(r_projected_kpc)
    for index, radius_projected in enumerate(r_projected_kpc):
        if radius_projected >= effective_r_max:
            sigma_msun_kpc2[index] = 0.0
            continue

        z_max = np.sqrt(max(effective_r_max**2 - radius_projected**2, 0.0))
        z_grid = np.linspace(0.0, z_max, n_z)
        r_grid = np.sqrt(radius_projected**2 + z_grid**2)

        log_rho = interpolation(np.log(r_grid))
        rho_grid = np.exp(log_rho)
        rho_grid[~np.isfinite(rho_grid)] = 0.0

        sigma_msun_kpc2[index] = 2.0 * np.trapz(rho_grid, z_grid)

    return sigma_msun_kpc2


def sigma_bar_of_R(r_projected_kpc: np.ndarray, sigma_msun_kpc2: np.ndarray) -> np.ndarray:
    """Compute mean enclosed Sigma_bar(<R) in Msun/kpc^2."""
    r_projected_kpc = np.asarray(r_projected_kpc)
    sigma_msun_kpc2 = np.asarray(sigma_msun_kpc2)

    if not np.all(np.diff(r_projected_kpc) > 0.0):
        raise ValueError("Projected radii must be strictly increasing.")

    r_extended = np.concatenate(([0.0], r_projected_kpc))
    sigma_extended = np.concatenate(([sigma_msun_kpc2[0]], sigma_msun_kpc2))

    integrand = sigma_extended * r_extended
    cumulative = np.zeros_like(r_extended)
    cumulative[1:] = np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(r_extended))

    sigma_bar_msun_kpc2 = 2.0 * cumulative[1:] / (r_projected_kpc**2)
    return sigma_bar_msun_kpc2


def delta_sigma_of_R(
    r_projected_kpc: np.ndarray,
    r_kpc: np.ndarray,
    rho_msun_kpc3: np.ndarray,
    n_z: int = 800,
    r_max_kpc: float | None = None,
    n_inner_extension: int = 64,
) -> dict[str, np.ndarray]:
    """Compute Sigma(R), Sigma_bar(<R), and DeltaSigma(R) in Msun/kpc^2.

    For stable Sigma_bar at the smallest requested radius, this function evaluates
    Sigma on an internally-extended R-grid below min(R) and interpolates back.
    """
    r_projected_kpc = np.asarray(r_projected_kpc)
    if not np.all(np.diff(r_projected_kpc) > 0.0):
        raise ValueError("Projected radii must be strictly increasing.")

    r_min_requested = r_projected_kpc[0]
    r_min_physical = max(1.0e-6, np.min(r_kpc) * 0.2)

    if r_min_physical < r_min_requested:
        inner_grid = np.geomspace(
            r_min_physical, r_min_requested, n_inner_extension, endpoint=False
        )
        r_work_kpc = np.concatenate((inner_grid, r_projected_kpc))
    else:
        r_work_kpc = r_projected_kpc

    sigma_work = sigma_of_R(
        r_projected_kpc=r_work_kpc,
        r_kpc=r_kpc,
        rho_msun_kpc3=rho_msun_kpc3,
        n_z=n_z,
        r_max_kpc=r_max_kpc,
    )
    sigma_bar_work = sigma_bar_of_R(r_work_kpc, sigma_work)
    delta_sigma_work = sigma_bar_work - sigma_work

    sigma_msun_kpc2 = np.interp(r_projected_kpc, r_work_kpc, sigma_work)
    sigma_bar_msun_kpc2 = np.interp(r_projected_kpc, r_work_kpc, sigma_bar_work)
    delta_sigma_msun_kpc2 = np.interp(r_projected_kpc, r_work_kpc, delta_sigma_work)

    return {
        "r_projected_kpc": r_projected_kpc,
        "sigma_msun_kpc2": sigma_msun_kpc2,
        "sigma_bar_msun_kpc2": sigma_bar_msun_kpc2,
        "delta_sigma_msun_kpc2": delta_sigma_msun_kpc2,
    }


def analytic_nfw_delta_sigma(
    r_projected_kpc: np.ndarray,
    rho_s_msun_kpc3: float,
    rs_kpc: float,
) -> dict[str, np.ndarray]:
    """Analytic NFW Sigma and DeltaSigma (Wright & Brainerd 2000 forms)."""
    x = np.asarray(r_projected_kpc) / rs_kpc

    sigma_factor = 2.0 * rho_s_msun_kpc3 * rs_kpc
    sigma = np.zeros_like(x)
    sigma_bar = np.zeros_like(x)

    less_than_one = x < 1.0
    equal_one = np.isclose(x, 1.0)
    greater_than_one = x > 1.0

    if np.any(less_than_one):
        x_l = x[less_than_one]
        sqrt_term = np.sqrt((1.0 - x_l) / (1.0 + x_l))
        atanh_term = np.arctanh(sqrt_term)
        denominator = np.sqrt(1.0 - x_l**2)
        sigma[less_than_one] = (
            sigma_factor * (1.0 - (2.0 / denominator) * atanh_term) / (x_l**2 - 1.0)
        )
        sigma_bar[less_than_one] = (
            4.0
            * rho_s_msun_kpc3
            * rs_kpc
            * (np.log(x_l / 2.0) + (2.0 / denominator) * atanh_term)
            / x_l**2
        )

    if np.any(equal_one):
        sigma[equal_one] = sigma_factor / 3.0
        sigma_bar[equal_one] = 4.0 * rho_s_msun_kpc3 * rs_kpc * (1.0 + np.log(0.5))

    if np.any(greater_than_one):
        x_g = x[greater_than_one]
        sqrt_term = np.sqrt((x_g - 1.0) / (1.0 + x_g))
        atan_term = np.arctan(sqrt_term)
        denominator = np.sqrt(x_g**2 - 1.0)
        sigma[greater_than_one] = (
            sigma_factor * (1.0 - (2.0 / denominator) * atan_term) / (x_g**2 - 1.0)
        )
        sigma_bar[greater_than_one] = (
            4.0
            * rho_s_msun_kpc3
            * rs_kpc
            * (np.log(x_g / 2.0) + (2.0 / denominator) * atan_term)
            / x_g**2
        )

    return {
        "sigma_msun_kpc2": sigma,
        "sigma_bar_msun_kpc2": sigma_bar,
        "delta_sigma_msun_kpc2": sigma_bar - sigma,
    }
