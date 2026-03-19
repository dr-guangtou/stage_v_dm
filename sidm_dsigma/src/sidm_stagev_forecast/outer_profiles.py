"""DK14-like outer profile components for Tier-2 hybrid forecasts."""

from __future__ import annotations

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
