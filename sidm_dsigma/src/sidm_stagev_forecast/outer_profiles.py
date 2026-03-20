"""DK14-like outer profile components for Tier-2 hybrid forecasts."""

from __future__ import annotations

import os
import warnings
from pathlib import Path
from typing import Any

import numpy as np

from sidm_stagev_forecast.config import CosmologyConfig
from sidm_stagev_forecast.cosmology import rdelta, rho_mean_matter_z
from sidm_stagev_forecast.profiles import nfw_profile_from_m_c


def default_dk14_like_parameters(regime: str | None = None) -> dict[str, float]:
    """Return lightweight DK14-like parameter defaults for a regime."""
    regime_normalized = "cluster" if regime is None else regime.lower()
    if regime_normalized == "cluster":
        return {
            "r_t_fraction_r200m": 1.2,
            "beta": 4.0,
            "gamma": 6.0,
            "outer_norm": 2.2,
            "outer_slope": 1.6,
            "outer_pivot_fraction_r200m": 1.0,
        }
    if regime_normalized == "dwarf":
        return {
            "r_t_fraction_r200m": 1.0,
            "beta": 4.0,
            "gamma": 4.5,
            "outer_norm": 1.4,
            "outer_slope": 1.8,
            "outer_pivot_fraction_r200m": 1.0,
        }
    raise ValueError("regime must be 'cluster', 'dwarf', or None.")


def build_dk14_outer_profile(
    r_kpc: np.ndarray,
    mass_msun: float,
    z: float,
    concentration: float,
    mass_def: str,
    cosmology: CosmologyConfig,
    dk14_params: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Construct a DK14-inspired profile bundle on an input radial grid.

    The profile form follows a practical approximation:
    rho_total(r) = rho_inner_baseline(r) * f_trans(r) + rho_outer_term(r)
    """
    radius_kpc = np.asarray(r_kpc, dtype=float)
    if np.any(radius_kpc <= 0.0):
        raise ValueError("All radii must be positive.")

    mass_definition = mass_def.lower()
    if mass_definition != "200c":
        raise ValueError("Tier-2 currently assumes mass_def='200c' halo inputs.")

    regime = "cluster"
    if dk14_params is not None and "regime" in dk14_params:
        regime = str(dk14_params["regime"])
    parameter_values = default_dk14_like_parameters(regime=regime)
    if dk14_params is not None:
        for key, value in dk14_params.items():
            if key == "regime":
                continue
            parameter_values[key] = float(value)

    r200c_kpc = rdelta(mass_msun, z, cosmology, definition="200c")
    r200m_kpc = rdelta(mass_msun, z, cosmology, definition="200m")

    rho_inner_baseline = nfw_profile_from_m_c(
        r_kpc=radius_kpc,
        m200_msun=mass_msun,
        c200=concentration,
        z=z,
        cosmo=cosmology,
    )

    r_t_kpc = parameter_values["r_t_fraction_r200m"] * r200m_kpc
    beta = max(parameter_values["beta"], 0.1)
    gamma = max(parameter_values["gamma"], 0.1)
    f_trans = (1.0 + (radius_kpc / r_t_kpc) ** beta) ** (-gamma / beta)

    rho_matter_msun_kpc3 = rho_mean_matter_z(z, cosmology)
    pivot_radius_kpc = max(
        parameter_values["outer_pivot_fraction_r200m"] * r200m_kpc,
        1.0e-6,
    )
    outer_slope = max(parameter_values["outer_slope"], 0.1)
    rho_outer_term = (
        parameter_values["outer_norm"]
        * rho_matter_msun_kpc3
        * (radius_kpc / pivot_radius_kpc) ** (-outer_slope)
    )

    rho_total = rho_inner_baseline * f_trans + rho_outer_term

    return {
        "r_kpc": radius_kpc,
        "rho_total_msun_kpc3": rho_total,
        "rho_inner_baseline_msun_kpc3": rho_inner_baseline,
        "rho_outer_term_msun_kpc3": rho_outer_term,
        "f_trans": f_trans,
        "metadata": {
            "model": "dk14_like",
            "mass_definition_input": mass_definition,
            "m200_msun": float(mass_msun),
            "z": float(z),
            "c200": float(concentration),
            "r200c_kpc": float(r200c_kpc),
            "r200m_kpc": float(r200m_kpc),
            "parameters": parameter_values,
        },
    }


def _configure_colossus_runtime(cosmology: CosmologyConfig) -> float:
    """Configure colossus cache/cosmology and return h."""
    from colossus import settings
    from colossus.cosmology import cosmology as colossus_cosmology

    cache_base = Path(os.environ.get("COLOSSUS_CACHE_DIR", ".colossus_runtime"))
    cache_base.mkdir(parents=True, exist_ok=True)
    settings.BASE_DIR = str(cache_base.resolve())
    settings.PERSISTENCE = "rw"

    cosmology_name = "sidm_stagev_planck18"
    if cosmology_name not in colossus_cosmology.cosmologies:
        colossus_cosmology.addCosmology(
            cosmology_name,
            {
                "flat": True,
                "H0": cosmology.h * 100.0,
                "Om0": cosmology.omega_m,
                "Ob0": 0.049,
                "sigma8": 0.811,
                "ns": 0.965,
            },
        )
    colossus_cosmology.setCosmology(cosmology_name)
    return float(cosmology.h)


def build_colossus_diemer23_outer_profile(
    r_kpc: np.ndarray,
    mass_msun: float,
    z: float,
    concentration: float,
    mass_def: str,
    cosmology: CosmologyConfig,
    diemer23_params: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Construct an optional colossus Diemer23+power-law outer profile.

    This backend remains optional. If colossus import or runtime evaluation fails,
    the function falls back to the local DK14-like approximation.
    """
    radius_kpc = np.asarray(r_kpc, dtype=float)
    if np.any(radius_kpc <= 0.0):
        raise ValueError("All radii must be positive.")

    mass_definition = mass_def.lower()
    if mass_definition != "200c":
        raise ValueError("Tier-2 currently assumes mass_def='200c' halo inputs.")

    regime = "cluster"
    if diemer23_params is not None and "regime" in diemer23_params:
        regime = str(diemer23_params["regime"])
    default_parameters = default_dk14_like_parameters(regime=regime)
    parameter_values = dict(default_parameters)
    if diemer23_params is not None:
        for key, value in diemer23_params.items():
            if key == "regime":
                continue
            parameter_values[key] = float(value)

    # Map DK14-like controls to the colossus outer power-law term.
    outer_norm = float(parameter_values.get("outer_norm", default_parameters["outer_norm"]))
    pivot_fraction_r200m = float(
        parameter_values.get(
            "r_t_fraction_r200m",
            parameter_values.get("outer_pivot_fraction_r200m", default_parameters["outer_pivot_fraction_r200m"]),
        )
    )
    slope_base = float(parameter_values.get("outer_slope", default_parameters["outer_slope"]))
    gamma_ratio = float(parameter_values.get("gamma", default_parameters["gamma"])) / float(
        default_parameters["gamma"]
    )
    beta_ratio = float(parameter_values.get("beta", default_parameters["beta"])) / float(
        default_parameters["beta"]
    )
    mapped_outer_slope = max(slope_base * np.sqrt(max(gamma_ratio, 1.0e-3) * max(beta_ratio, 1.0e-3)), 0.1)

    try:
        from colossus.halo import profile_composite

        h = _configure_colossus_runtime(cosmology)
        profile = profile_composite.compositeProfile(
            inner_name="diemer23",
            outer_names=["pl"],
            M=mass_msun * h,
            c=concentration,
            z=z,
            mdef="200c",
            norm=outer_norm,
            slope=mapped_outer_slope,
            pivot="R200m",
            pivot_factor=pivot_fraction_r200m,
        )

        radius_kpc_h = radius_kpc * h
        rho_total_msun_kpc3 = np.asarray(profile.density(radius_kpc_h), dtype=float) * (h**2)
        rho_inner_msun_kpc3 = np.asarray(profile.densityInner(radius_kpc_h), dtype=float) * (h**2)
        rho_outer_msun_kpc3 = np.asarray(profile.densityOuter(radius_kpc_h), dtype=float) * (h**2)
        f_trans = np.clip((rho_total_msun_kpc3 - rho_outer_msun_kpc3) / np.maximum(rho_inner_msun_kpc3, 1.0e-40), 0.0, np.inf)
        r200c_kpc = rdelta(mass_msun, z, cosmology, definition="200c")
        r200m_kpc = rdelta(mass_msun, z, cosmology, definition="200m")

        return {
            "r_kpc": radius_kpc,
            "rho_total_msun_kpc3": rho_total_msun_kpc3,
            "rho_inner_baseline_msun_kpc3": rho_inner_msun_kpc3,
            "rho_outer_term_msun_kpc3": rho_outer_msun_kpc3,
            "f_trans": f_trans,
            "metadata": {
                "model": "colossus_diemer23",
                "mass_definition_input": mass_definition,
                "m200_msun": float(mass_msun),
                "z": float(z),
                "c200": float(concentration),
                "r200c_kpc": float(r200c_kpc),
                "r200m_kpc": float(r200m_kpc),
                "parameters": parameter_values,
                "mapped_outer_norm": outer_norm,
                "mapped_outer_slope": mapped_outer_slope,
                "mapped_pivot_fraction_r200m": pivot_fraction_r200m,
            },
        }
    except Exception as error:
        warnings.warn(
            "colossus_diemer23 outer profile evaluation failed; falling back to local dk14_like. "
            f"Reason: {error!r}",
            RuntimeWarning,
            stacklevel=2,
        )
        return build_dk14_outer_profile(
            r_kpc=radius_kpc,
            mass_msun=mass_msun,
            z=z,
            concentration=concentration,
            mass_def=mass_def,
            cosmology=cosmology,
            dk14_params=diemer23_params,
        )


def build_outer_profile(
    r_kpc: np.ndarray,
    mass_msun: float,
    z: float,
    concentration: float,
    mass_def: str,
    cosmology: CosmologyConfig,
    outer_profile_model: str,
    outer_params: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Build the selected outer-profile model with a stable interface."""
    model_name = str(outer_profile_model).lower()
    if model_name == "dk14_like":
        return build_dk14_outer_profile(
            r_kpc=r_kpc,
            mass_msun=mass_msun,
            z=z,
            concentration=concentration,
            mass_def=mass_def,
            cosmology=cosmology,
            dk14_params=outer_params,
        )
    if model_name in {"colossus_diemer23", "diemer23"}:
        return build_colossus_diemer23_outer_profile(
            r_kpc=r_kpc,
            mass_msun=mass_msun,
            z=z,
            concentration=concentration,
            mass_def=mass_def,
            cosmology=cosmology,
            diemer23_params=outer_params,
        )
    raise ValueError(
        "Unsupported outer_profile_model. Supported models are 'dk14_like' and "
        f"'colossus_diemer23'; received {outer_profile_model!r}."
    )
