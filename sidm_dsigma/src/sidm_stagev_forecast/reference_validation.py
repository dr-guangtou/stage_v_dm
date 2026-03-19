"""Optional reference cross-check utilities for NFW lensing projections."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from sidm_stagev_forecast.config import DEFAULT_COSMOLOGY
from sidm_stagev_forecast.profiles import nfw_parameters_from_m_c, nfw_profile_from_m_c
from sidm_stagev_forecast.projection import analytic_nfw_delta_sigma, delta_sigma_of_R


@dataclass(frozen=True)
class BackendStatus:
    """Availability and execution status for an optional reference backend."""

    backend: str
    available: bool
    used: bool
    message: str


def detect_reference_backends() -> dict[str, bool]:
    """Report availability of optional reference backends."""
    availability = {"colossus": False, "pyccl": False}

    try:
        import colossus  # noqa: F401

        availability["colossus"] = True
    except Exception:
        availability["colossus"] = False

    try:
        import pyccl  # noqa: F401

        availability["pyccl"] = True
    except Exception:
        availability["pyccl"] = False

    return availability


def _compute_local_projection(
    m200_msun: float,
    c200: float,
    z: float,
    r_3d_kpc: np.ndarray,
    r_projected_kpc: np.ndarray,
) -> dict[str, np.ndarray]:
    rho_msun_kpc3 = nfw_profile_from_m_c(r_3d_kpc, m200_msun, c200, z, DEFAULT_COSMOLOGY)
    return delta_sigma_of_R(r_projected_kpc, r_3d_kpc, rho_msun_kpc3, n_z=1200)


def _compute_analytic_reference(
    m200_msun: float,
    c200: float,
    z: float,
    r_projected_kpc: np.ndarray,
) -> dict[str, np.ndarray]:
    parameters = nfw_parameters_from_m_c(m200_msun, c200, z, DEFAULT_COSMOLOGY)
    return analytic_nfw_delta_sigma(r_projected_kpc, parameters.rho_s_msun_kpc3, parameters.rs_kpc)


def _compute_colossus_reference(
    m200_msun: float,
    c200: float,
    z: float,
    r_projected_kpc: np.ndarray,
) -> dict[str, np.ndarray]:
    """Best-effort Colossus reference evaluation.

    Notes
    -----
    This adapter is optional and may fail if Colossus conventions/API differ from
    this local environment setup. Failures are handled by the caller.
    """
    from colossus.cosmology import cosmology as colossus_cosmology
    from colossus.halo import profile_nfw

    if "sidm_stagev_crosscheck" not in colossus_cosmology.cosmologies:
        colossus_cosmology.addCosmology(
            "sidm_stagev_crosscheck",
            {
                "flat": True,
                "H0": 100.0 * DEFAULT_COSMOLOGY.h,
                "Om0": DEFAULT_COSMOLOGY.omega_m,
                "Ob0": 0.049,
                "sigma8": 0.81,
                "ns": 0.965,
            },
        )
    colossus_cosmology.setCosmology("sidm_stagev_crosscheck")

    profile = profile_nfw.NFWProfile(M=m200_msun, c=c200, z=z, mdef="200c")

    sigma_msun_kpc2 = np.asarray(profile.surfaceDensity(r_projected_kpc))
    delta_sigma_msun_kpc2 = np.asarray(profile.deltaSigma(r_projected_kpc))
    sigma_bar_msun_kpc2 = sigma_msun_kpc2 + delta_sigma_msun_kpc2

    return {
        "sigma_msun_kpc2": sigma_msun_kpc2,
        "sigma_bar_msun_kpc2": sigma_bar_msun_kpc2,
        "delta_sigma_msun_kpc2": delta_sigma_msun_kpc2,
    }


def run_nfw_reference_crosscheck(
    m200_msun: float = 1.0e14,
    c200: float = 4.0,
    z: float = 0.3,
    r_3d_kpc: np.ndarray | None = None,
    r_projected_kpc: np.ndarray | None = None,
) -> tuple[pd.DataFrame, list[BackendStatus]]:
    """Build NFW cross-check table against analytic and optional external references."""
    if r_3d_kpc is None:
        r_3d_kpc = np.geomspace(0.5, 1.0e4, 600)
    if r_projected_kpc is None:
        r_projected_kpc = np.geomspace(30.0, 2.0e3, 45)

    local_projection = _compute_local_projection(m200_msun, c200, z, r_3d_kpc, r_projected_kpc)
    analytic_projection = _compute_analytic_reference(m200_msun, c200, z, r_projected_kpc)

    table = pd.DataFrame(
        {
            "r_projected_kpc": r_projected_kpc,
            "delta_sigma_local_msun_kpc2": local_projection["delta_sigma_msun_kpc2"],
            "delta_sigma_analytic_msun_kpc2": analytic_projection["delta_sigma_msun_kpc2"],
        }
    )
    table["frac_diff_local_vs_analytic"] = (
        table["delta_sigma_local_msun_kpc2"] / table["delta_sigma_analytic_msun_kpc2"] - 1.0
    )

    backend_status: list[BackendStatus] = []
    availability = detect_reference_backends()

    if availability["colossus"]:
        try:
            colossus_projection = _compute_colossus_reference(m200_msun, c200, z, r_projected_kpc)
            table["delta_sigma_colossus_msun_kpc2"] = colossus_projection["delta_sigma_msun_kpc2"]
            table["frac_diff_local_vs_colossus"] = (
                table["delta_sigma_local_msun_kpc2"] / table["delta_sigma_colossus_msun_kpc2"] - 1.0
            )
            backend_status.append(
                BackendStatus(
                    backend="colossus",
                    available=True,
                    used=True,
                    message="Colossus reference computed successfully.",
                )
            )
        except Exception as error:  # pragma: no cover - external dependency surface
            backend_status.append(
                BackendStatus(
                    backend="colossus",
                    available=True,
                    used=False,
                    message=f"Colossus available but cross-check failed: {error}",
                )
            )
    else:
        backend_status.append(
            BackendStatus(
                backend="colossus",
                available=False,
                used=False,
                message="Colossus is not installed in this environment.",
            )
        )

    if availability["pyccl"]:
        backend_status.append(
            BackendStatus(
                backend="pyccl",
                available=True,
                used=False,
                message="pyccl detected but adapter is not implemented in this phase.",
            )
        )
    else:
        backend_status.append(
            BackendStatus(
                backend="pyccl",
                available=False,
                used=False,
                message="pyccl is not installed in this environment.",
            )
        )

    return table, backend_status
