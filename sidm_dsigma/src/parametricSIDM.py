"""Local compatibility shim for Daneng Yang's parametricSIDM repository.

This module exposes a stable callable expected by the local forecast wrapper:
`density_profile_from_m200_c200`.

It loads `parametricC4.py` from a local parametricSIDM checkout and maps
`(M200c, c200, z, sigma_over_m)` to the C4 profile prescription.
"""

from __future__ import annotations

import importlib.util
import os
from functools import lru_cache
from pathlib import Path

import numpy as np

from sidm_stagev_forecast.config import DEFAULT_COSMOLOGY
from sidm_stagev_forecast.cosmology import rdelta


def _default_parametric_sidm_path() -> Path:
    repository_root = Path(__file__).resolve().parents[1]
    return repository_root / "third_party" / "parametricSIDM"


@lru_cache(maxsize=1)
def _load_parametric_c4_module():
    repo_path = Path(os.environ.get("PARAMETRIC_SIDM_REPO", str(_default_parametric_sidm_path())))
    module_path = repo_path / "parametricC4.py"

    if not module_path.exists():
        raise FileNotFoundError(
            "Could not find parametricSIDM module file parametricC4.py. "
            "Set PARAMETRIC_SIDM_REPO to your local parametricSIDM checkout."
        )

    spec = importlib.util.spec_from_file_location("parametric_c4_module", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Failed to create import spec for parametricC4.py")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _nfw_rho_s_rs_from_m200c(
    m200_msun: float,
    c200: float,
    z: float,
) -> tuple[float, float]:
    r200_kpc = rdelta(m200_msun, z, DEFAULT_COSMOLOGY, definition="200c")
    rs_kpc = r200_kpc / c200
    normalization = np.log(1.0 + c200) - c200 / (1.0 + c200)
    rho_s_msun_kpc3 = m200_msun / (4.0 * np.pi * rs_kpc**3 * normalization)
    return rho_s_msun_kpc3, rs_kpc


def density_profile_from_m200_c200(
    r_kpc: np.ndarray,
    m200_msun: float,
    c200: float,
    z: float,
    sigma_over_m: float,
    elapsed_time_gyr: float | None = None,
) -> np.ndarray:
    """Return SIDM density profile in Msun/kpc^3 using parametricSIDM C4 model."""
    r_kpc = np.asarray(r_kpc)
    rho_s_msun_kpc3, rs_kpc = _nfw_rho_s_rs_from_m200c(m200_msun, c200, z)

    if sigma_over_m <= 0.0:
        x = r_kpc / rs_kpc
        return rho_s_msun_kpc3 / (x * (1.0 + x) ** 2)

    c4_module = _load_parametric_c4_module()

    if elapsed_time_gyr is None:
        elapsed_time_gyr = float(c4_module.tlb(z))

    tc_gyr = float(c4_module.tc(float(sigma_over_m), float(rho_s_msun_kpc3), float(rs_kpc)))
    if tc_gyr <= 0.0:
        tr = 0.0
    else:
        tr = float(np.clip(elapsed_time_gyr / tc_gyr, 0.0, 1.1))

    rho_t = float(c4_module.rhost(tr, float(rho_s_msun_kpc3), float(rs_kpc)))
    rs_t = float(c4_module.rst(tr, float(rho_s_msun_kpc3), float(rs_kpc)))
    rc_t = float(c4_module.rct(tr, float(rho_s_msun_kpc3), float(rs_kpc)))

    return np.asarray(c4_module.frho(r_kpc, rho_t, rs_t, rc_t))
