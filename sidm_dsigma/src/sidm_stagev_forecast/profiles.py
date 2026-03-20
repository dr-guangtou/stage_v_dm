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

    metadata: dict[str, Any] = {
        "profile": "sidm",
        "mass_definition": "200c",
        "m200_msun": float(m200_msun),
        "c200": float(c200),
        "z": float(z),
        "sigma_over_m_cm2_g": float(sigma_over_m),
    }
    if "cross_section_parameterization" in sidm_kwargs:
        metadata["cross_section_parameterization"] = str(
            sidm_kwargs["cross_section_parameterization"]
        )
    if "w_km_s" in sidm_kwargs:
        metadata["w_km_s"] = float(sidm_kwargs["w_km_s"])
    if "time_model" in sidm_kwargs:
        metadata["time_model"] = str(sidm_kwargs["time_model"])

    return {
        "r_kpc": np.asarray(r_kpc),
        "rho_msun_kpc3": rho_msun_kpc3,
        "m_enclosed_msun": m_enclosed_msun,
        "vcirc_km_s": vcirc_km_s,
        "metadata": metadata,
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
    from sidm_stagev_forecast.outer_profiles import build_outer_profile
    from sidm_stagev_forecast.stitch import resolve_match_radius_kpc, stitch_inner_outer_profile

    cdm_profile_source = str((model_options or {}).get("cdm_profile_source", "nfw")).lower()
    if cdm_profile_source == "parametric":
        cdm_sigma_over_m = float((model_options or {}).get("cdm_sigma_over_m", 0.0))
        cdm_sidm_kwargs = dict((model_options or {}).get("cdm_sidm_kwargs", {}))
        cdm_inner = sidm_profile_from_parametric_model(
            r_kpc=r_kpc,
            m200_msun=m200_msun,
            c200=c200,
            z=z,
            sigma_over_m=cdm_sigma_over_m,
            model_options={"sidm_kwargs": cdm_sidm_kwargs},
        )
        cdm_inner["metadata"]["profile"] = "cdm_parametric"
    else:
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

    stitching_options = {} if stitch_params is None else dict(stitch_params)
    outer_profile_model = str(stitching_options.get("outer_profile_model", "dk14_like")).lower()

    outer_bundle = build_outer_profile(
        r_kpc=r_kpc,
        mass_msun=m200_msun,
        z=z,
        concentration=c200,
        mass_def="200c",
        cosmology=cosmo,
        outer_profile_model=outer_profile_model,
        outer_params=dk14_params,
    )

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
        "outer_reference": outer_bundle,
        "dk14_outer_reference": outer_bundle,
        "rho_cdm_hybrid_msun_kpc3": rho_cdm_hybrid,
        "rho_sidm_hybrid_msun_kpc3": rho_sidm_hybrid,
        "metadata": {
            "profile": "sidm_hybrid_tier2",
            "sigma_over_m_cm2_g": float(sigma_over_m),
            "r_match_kpc": float(r_match_kpc),
            "stitch_method": str(stitching_options.get("stitch_method", "logistic_logrho_blend")),
            "outer_profile_model": outer_profile_model,
            "mass_definition": "200c",
        },
    }


def build_tier3_sidm_profile(
    r_kpc: np.ndarray,
    m200_msun: float,
    c200: float,
    z: float,
    sigma_over_m: float,
    cosmo: CosmologyConfig,
    model_options: dict[str, Any] | None = None,
    dk14_params: dict[str, float] | None = None,
    stitch_params: dict[str, Any] | None = None,
    tier3_config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build Tier-3 profile: Tier-2 hybrid plus empirical outer correction."""
    from sidm_stagev_forecast.calibration import resolve_tier3_parameters
    from sidm_stagev_forecast.outer_corrections import (
        apply_sidm_outer_correction,
        modify_dk14_parameters_for_sidm,
    )
    from sidm_stagev_forecast.outer_profiles import build_outer_profile
    from sidm_stagev_forecast.stitch import resolve_match_radius_kpc, stitch_inner_outer_profile

    tier2_bundle = build_hybrid_sidm_profile(
        r_kpc=r_kpc,
        m200_msun=m200_msun,
        c200=c200,
        z=z,
        sigma_over_m=sigma_over_m,
        cosmo=cosmo,
        model_options=model_options,
        dk14_params=dk14_params,
        stitch_params=stitch_params,
    )

    config = {} if tier3_config is None else dict(tier3_config)
    if not bool(config.get("enabled", False)):
        return {
            **tier2_bundle,
            "rho_cdm_tier3_msun_kpc3": tier2_bundle["rho_cdm_hybrid_msun_kpc3"],
            "rho_sidm_tier3_msun_kpc3": tier2_bundle["rho_sidm_hybrid_msun_kpc3"],
            "metadata": {**tier2_bundle["metadata"], "profile": "sidm_tier3_disabled"},
        }

    from sidm_stagev_forecast.outer_profiles import default_dk14_like_parameters

    stitch_options = {} if stitch_params is None else dict(stitch_params)
    outer_profile_model = str(stitch_options.get("outer_profile_model", "dk14_like")).lower()
    regime = "cluster"
    if dk14_params is not None and "regime" in dk14_params:
        regime = str(dk14_params["regime"])
    base_dk14 = default_dk14_like_parameters(regime=regime)
    if dk14_params is not None:
        for key, value in dict(dk14_params).items():
            if key == "regime":
                continue
            base_dk14[key] = float(value)
    correction_parameters = resolve_tier3_parameters(config)
    correction_model = str(correction_parameters["correction_model"])

    halo_context = {
        "m200_msun": float(m200_msun),
        "z": float(z),
        "c200": float(c200),
    }
    sidm_context = {"sigma_over_m": float(sigma_over_m)}
    cdm_context = {"sigma_over_m": 0.0}

    dk14_sidm = modify_dk14_parameters_for_sidm(
        halo=halo_context,
        sidm_params=sidm_context,
        base_dk14_params=base_dk14,
        correction_model=correction_model,
        correction_params=correction_parameters,
        cosmology=cosmo,
    )
    dk14_cdm = modify_dk14_parameters_for_sidm(
        halo=halo_context,
        sidm_params=cdm_context,
        base_dk14_params=base_dk14,
        correction_model=correction_model,
        correction_params=correction_parameters,
        cosmology=cosmo,
    )

    outer_sidm = build_outer_profile(
        r_kpc=r_kpc,
        mass_msun=m200_msun,
        z=z,
        concentration=c200,
        mass_def="200c",
        cosmology=cosmo,
        outer_profile_model=outer_profile_model,
        outer_params=dk14_sidm,
    )
    outer_cdm = build_outer_profile(
        r_kpc=r_kpc,
        mass_msun=m200_msun,
        z=z,
        concentration=c200,
        mass_def="200c",
        cosmology=cosmo,
        outer_profile_model=outer_profile_model,
        outer_params=dk14_cdm,
    )

    r_match_kpc = resolve_match_radius_kpc(
        mass_msun=m200_msun,
        z=z,
        cosmo=cosmo,
        stitch_config=stitch_options,
    )
    rho_cdm_tier3 = stitch_inner_outer_profile(
        r_kpc=r_kpc,
        rho_inner_sidm_msun_kpc3=tier2_bundle["cdm_inner"]["rho_msun_kpc3"],
        rho_outer_reference_msun_kpc3=outer_cdm["rho_total_msun_kpc3"],
        method=str(stitch_options.get("stitch_method", "logistic_logrho_blend")),
        r_match_kpc=r_match_kpc,
        smooth_width_dex=float(stitch_options.get("smooth_width_dex", 0.15)),
        continuity=str(stitch_options.get("continuity", "density")),
    )
    rho_sidm_tier3 = stitch_inner_outer_profile(
        r_kpc=r_kpc,
        rho_inner_sidm_msun_kpc3=tier2_bundle["sidm_inner"]["rho_msun_kpc3"],
        rho_outer_reference_msun_kpc3=outer_sidm["rho_total_msun_kpc3"],
        method=str(stitch_options.get("stitch_method", "logistic_logrho_blend")),
        r_match_kpc=r_match_kpc,
        smooth_width_dex=float(stitch_options.get("smooth_width_dex", 0.15)),
        continuity=str(stitch_options.get("continuity", "density")),
    )

    rho_cdm_tier3 = apply_sidm_outer_correction(
        r_kpc=r_kpc,
        rho_tier2_msun_kpc3=rho_cdm_tier3,
        halo=halo_context,
        sidm_params=cdm_context,
        correction_model=correction_model,
        correction_params=correction_parameters,
        cosmology=cosmo,
    )
    rho_sidm_tier3 = apply_sidm_outer_correction(
        r_kpc=r_kpc,
        rho_tier2_msun_kpc3=rho_sidm_tier3,
        halo=halo_context,
        sidm_params=sidm_context,
        correction_model=correction_model,
        correction_params=correction_parameters,
        cosmology=cosmo,
    )

    return {
        **tier2_bundle,
        "rho_cdm_tier3_msun_kpc3": rho_cdm_tier3,
        "rho_sidm_tier3_msun_kpc3": rho_sidm_tier3,
        "dk14_outer_sidm_tier3": outer_sidm,
        "dk14_outer_cdm_tier3": outer_cdm,
        "metadata": {
            **tier2_bundle["metadata"],
            "profile": "sidm_tier3_empirical",
            "outer_profile_model": outer_profile_model,
            "tier3_correction_model": correction_model,
            "tier3_parameters": correction_parameters,
        },
    }
