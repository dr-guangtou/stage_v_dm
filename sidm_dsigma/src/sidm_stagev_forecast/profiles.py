"""Profile generators for CDM (NFW) and SIDM (parametricSIDM wrapper)."""

from __future__ import annotations

import importlib
from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.integrate import cumulative_trapezoid

from sidm_stagev_forecast.config import CosmologyConfig
from sidm_stagev_forecast.cosmology import GRAVITATIONAL_CONSTANT_KPC_KMS, rdelta


@dataclass(frozen=True)
class NfwParameters:
    """NFW structural parameters in physical units."""

    r200_kpc: float
    rs_kpc: float
    rho_s_msun_kpc3: float


def nfw_parameters_from_m_c(
    m200_msun: float,
    c200: float,
    z: float,
    cosmo: CosmologyConfig,
) -> NfwParameters:
    """Compute physical NFW parameters for an M200c halo."""
    r200_kpc = rdelta(m200_msun, z, cosmo, definition="200c")
    rs_kpc = r200_kpc / c200
    normalization = np.log(1.0 + c200) - c200 / (1.0 + c200)
    rho_s_msun_kpc3 = m200_msun / (4.0 * np.pi * rs_kpc**3 * normalization)
    return NfwParameters(r200_kpc=r200_kpc, rs_kpc=rs_kpc, rho_s_msun_kpc3=rho_s_msun_kpc3)


def nfw_profile_from_m_c(
    r_kpc: np.ndarray,
    m200_msun: float,
    c200: float,
    z: float,
    cosmo: CosmologyConfig,
) -> np.ndarray:
    """Return NFW rho(r) in Msun/kpc^3."""
    params = nfw_parameters_from_m_c(m200_msun, c200, z, cosmo)
    x = np.asarray(r_kpc) / params.rs_kpc
    return params.rho_s_msun_kpc3 / (x * (1.0 + x) ** 2)


def nfw_m_enclosed(r_kpc: np.ndarray, params: NfwParameters) -> np.ndarray:
    """Return enclosed mass M(<r) in Msun for NFW."""
    x = np.asarray(r_kpc) / params.rs_kpc
    return (
        4.0 * np.pi * params.rho_s_msun_kpc3 * params.rs_kpc**3 * (np.log(1.0 + x) - x / (1.0 + x))
    )


def circular_velocity_km_s(r_kpc: np.ndarray, m_enclosed_msun: np.ndarray) -> np.ndarray:
    """Return circular velocity in km/s."""
    return np.sqrt(GRAVITATIONAL_CONSTANT_KPC_KMS * m_enclosed_msun / np.asarray(r_kpc))


def nfw_profile_bundle(
    r_kpc: np.ndarray,
    m200_msun: float,
    c200: float,
    z: float,
    cosmo: CosmologyConfig,
) -> dict[str, Any]:
    """Return standard profile dictionary for CDM NFW benchmark."""
    params = nfw_parameters_from_m_c(m200_msun, c200, z, cosmo)
    rho_msun_kpc3 = nfw_profile_from_m_c(r_kpc, m200_msun, c200, z, cosmo)
    m_enclosed_msun = nfw_m_enclosed(r_kpc, params)
    vcirc_km_s = circular_velocity_km_s(r_kpc, m_enclosed_msun)
    return {
        "r_kpc": np.asarray(r_kpc),
        "rho_msun_kpc3": rho_msun_kpc3,
        "m_enclosed_msun": m_enclosed_msun,
        "vcirc_km_s": vcirc_km_s,
        "metadata": {
            "profile": "nfw",
            "mass_definition": "200c",
            "m200_msun": float(m200_msun),
            "c200": float(c200),
            "z": float(z),
        },
    }


def _infer_parametric_sidm_function(module: Any) -> Any:
    """Resolve a callable from parametricSIDM for rho(r) generation.

    The upstream repository API may evolve. This resolver supports two paths:
    1) user-provided callable via model_options['sidm_callable']
    2) module-level function named `density_profile_from_m200_c200`
    """
    if hasattr(module, "density_profile_from_m200_c200"):
        return module.density_profile_from_m200_c200
    raise AttributeError(
        "Could not find density_profile_from_m200_c200 in parametricSIDM. "
        "Pass model_options['sidm_callable'] to sidm_profile_from_parametric_model."
    )


def sidm_profile_from_parametric_model(
    r_kpc: np.ndarray,
    m200_msun: float,
    c200: float,
    z: float,
    sigma_over_m: float,
    model_options: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Return SIDM profile bundle via a stable local wrapper.

    Notes
    -----
    - If `model_options['sidm_callable']` is provided, it is used directly.
    - Otherwise this function tries to import `parametricSIDM` and resolve
      `density_profile_from_m200_c200`.
    - Expected callable signature:
      sidm_callable(r_kpc, m200_msun, c200, z, sigma_over_m, **kwargs) -> rho_msun_kpc3
    """
    model_options = model_options or {}

    sidm_callable = model_options.get("sidm_callable")
    if sidm_callable is None:
        try:
            sidm_module = importlib.import_module("parametricSIDM")
        except ModuleNotFoundError as exc:
            raise RuntimeError(
                "parametricSIDM is not installed or not importable. "
                "Install Daneng Yang's parametricSIDM or pass model_options['sidm_callable']."
            ) from exc
        sidm_callable = _infer_parametric_sidm_function(sidm_module)

    sidm_kwargs = dict(model_options.get("sidm_kwargs", {}))
    rho_msun_kpc3 = np.asarray(
        sidm_callable(
            np.asarray(r_kpc),
            float(m200_msun),
            float(c200),
            float(z),
            float(sigma_over_m),
            **sidm_kwargs,
        )
    )

    if rho_msun_kpc3.shape != np.asarray(r_kpc).shape:
        raise RuntimeError("SIDM wrapper received a density array with unexpected shape.")

    m_integrand = 4.0 * np.pi * np.asarray(r_kpc) ** 2 * rho_msun_kpc3
    m_enclosed_msun = cumulative_trapezoid(m_integrand, np.asarray(r_kpc), initial=0.0)
    vcirc_km_s = circular_velocity_km_s(np.asarray(r_kpc), m_enclosed_msun)

    return {
        "r_kpc": np.asarray(r_kpc),
        "rho_msun_kpc3": rho_msun_kpc3,
        "m_enclosed_msun": m_enclosed_msun,
        "vcirc_km_s": vcirc_km_s,
        "metadata": {
            "profile": "sidm",
            "mass_definition": "200c",
            "m200_msun": float(m200_msun),
            "c200": float(c200),
            "z": float(z),
            "sigma_over_m_cm2_g": float(sigma_over_m),
        },
    }


def build_hybrid_sidm_profile(
    r_kpc: np.ndarray,
    m200_msun: float,
    c200: float,
    z: float,
    sigma_over_m: float,
    cosmo: CosmologyConfig,
    model_options: dict[str, Any] | None = None,
    dk14_params: dict[str, float] | None = None,
    stitch_params: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build Tier-2 hybrid profile components and stitched outputs."""
    from sidm_stagev_forecast.outer_profiles import build_dk14_outer_profile
    from sidm_stagev_forecast.stitch import resolve_match_radius_kpc, stitch_inner_outer_profile

    cdm_inner = nfw_profile_bundle(
        r_kpc=r_kpc,
        m200_msun=m200_msun,
        c200=c200,
        z=z,
        cosmo=cosmo,
    )
    sidm_inner = sidm_profile_from_parametric_model(
        r_kpc=r_kpc,
        m200_msun=m200_msun,
        c200=c200,
        z=z,
        sigma_over_m=sigma_over_m,
        model_options=model_options,
    )

    outer_bundle = build_dk14_outer_profile(
        r_kpc=r_kpc,
        mass_msun=m200_msun,
        z=z,
        concentration=c200,
        mass_def="200c",
        cosmology=cosmo,
        dk14_params=dk14_params,
    )

    stitching_options = {} if stitch_params is None else dict(stitch_params)
    r_match_kpc = resolve_match_radius_kpc(
        mass_msun=m200_msun,
        z=z,
        cosmo=cosmo,
        stitch_config=stitching_options,
    )

    rho_cdm_hybrid = stitch_inner_outer_profile(
        r_kpc=r_kpc,
        rho_inner_sidm_msun_kpc3=cdm_inner["rho_msun_kpc3"],
        rho_outer_reference_msun_kpc3=outer_bundle["rho_total_msun_kpc3"],
        method=str(stitching_options.get("stitch_method", "logistic_logrho_blend")),
        r_match_kpc=r_match_kpc,
        smooth_width_dex=float(stitching_options.get("smooth_width_dex", 0.15)),
        continuity=str(stitching_options.get("continuity", "density")),
    )
    rho_sidm_hybrid = stitch_inner_outer_profile(
        r_kpc=r_kpc,
        rho_inner_sidm_msun_kpc3=sidm_inner["rho_msun_kpc3"],
        rho_outer_reference_msun_kpc3=outer_bundle["rho_total_msun_kpc3"],
        method=str(stitching_options.get("stitch_method", "logistic_logrho_blend")),
        r_match_kpc=r_match_kpc,
        smooth_width_dex=float(stitching_options.get("smooth_width_dex", 0.15)),
        continuity=str(stitching_options.get("continuity", "density")),
    )

    return {
        "r_kpc": np.asarray(r_kpc),
        "cdm_inner": cdm_inner,
        "sidm_inner": sidm_inner,
        "dk14_outer_reference": outer_bundle,
        "rho_cdm_hybrid_msun_kpc3": rho_cdm_hybrid,
        "rho_sidm_hybrid_msun_kpc3": rho_sidm_hybrid,
        "metadata": {
            "profile": "sidm_hybrid_tier2",
            "sigma_over_m_cm2_g": float(sigma_over_m),
            "r_match_kpc": float(r_match_kpc),
            "stitch_method": str(stitching_options.get("stitch_method", "logistic_logrho_blend")),
            "mass_definition": "200c",
        },
    }
