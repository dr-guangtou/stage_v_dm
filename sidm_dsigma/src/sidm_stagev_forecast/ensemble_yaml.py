"""YAML parsing and normalization for configurable halo-ensemble runs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class EnsembleYamlConfig:
    """Normalized runtime configuration parsed from a YAML file."""

    label: str
    mode: str
    ensemble_config: dict[str, Any]
    projection_config: dict[str, Any]
    sidm_sigma_grid: tuple[float, ...]
    sidm_config: dict[str, Any]
    tier2_config: dict[str, Any]
    tier3_config: dict[str, Any]


def _build_redshift_configuration(redshift_block: dict[str, Any]) -> dict[str, Any]:
    model = str(redshift_block.get("model", "fixed")).lower()
    if model == "gaussian":
        return {
            "redshift_distribution": {
                "type": "gaussian",
                "mean": float(redshift_block["mean"]),
                "sigma": float(redshift_block["sigma"]),
                "z_min": float(redshift_block["z_min"]),
                "z_max": float(redshift_block["z_max"]),
            }
        }
    if model == "fixed":
        return {"redshift": float(redshift_block["value"])}
    raise ValueError(f"Unsupported redshift model: {model}")


def _normalize_hmf_mode(raw: dict[str, Any], base_config: dict[str, Any]) -> dict[str, Any]:
    mass_function = dict(raw["mass_function"])
    selection = dict(raw["selection"])
    concentration = dict(raw["concentration"])

    ensemble_config = {
        **base_config,
        "mass_min_msun": float(mass_function["M_min"]),
        "mass_max_msun": float(mass_function["M_max"]),
        "n_mass_grid": int(mass_function.get("n_mass_grid", 512)),
        "concentration_scatter_dex": float(concentration["scatter_logc"]),
        "weight_mode": str(raw["stacking"]["weight_scheme"]).lower(),
    }

    hmf_model_name = str(mass_function.get("model", "power_law")).lower()
    if hmf_model_name == "tinker08":
        ensemble_config["hmf_model"] = {"type": "tinker08"}
    else:
        ensemble_config["hmf_model"] = {"type": "power_law", "alpha": 1.9}

    selection_type = str(selection.get("type", "none")).lower()
    if selection_type == "logistic":
        ensemble_config["selection_model"] = {
            "type": "logistic",
            "log10_m_cut": float(selection["logM_cut"]),
            "sigma_log10_m": float(selection["sigma_sel_dex"]),
        }
    elif selection_type == "threshold":
        ensemble_config["selection_model"] = {
            "type": "threshold",
            "log10_m_cut": float(selection["logM_cut"]),
        }
    else:
        ensemble_config["selection_model"] = {"type": "none"}

    concentration_name = str(concentration.get("model", "maccio")).lower()
    if concentration_name == "diemerjoyce19":
        ensemble_config["concentration_model"] = {"type": "diemerjoyce19"}
    else:
        ensemble_config["concentration_model"] = {"type": "maccio"}

    return ensemble_config


def _normalize_shmr_mode(raw: dict[str, Any], base_config: dict[str, Any]) -> dict[str, Any]:
    stellar_mass_distribution = dict(raw["stellar_mass_distribution"])
    shmr_block = dict(raw["shmr"])
    concentration = dict(raw["concentration"])

    ensemble_config = {
        **base_config,
        "stellar_mass_distribution": {
            "type": "lognormal",
            "mean_log10_mstar": float(stellar_mass_distribution["logMstar_mean"]),
            "sigma_log10_mstar": float(stellar_mass_distribution["logMstar_sigma"]),
            "log10_mstar_min": float(stellar_mass_distribution["logMstar_min"]),
            "log10_mstar_max": float(stellar_mass_distribution["logMstar_max"]),
        },
        "halo_scatter_dex": float(shmr_block["scatter_logMhalo"]),
        "concentration_scatter_dex": float(concentration["scatter_logc"]),
        "weight_mode": str(raw["stacking"]["weight_scheme"]).lower(),
        "mhalo_min_msun": 1.0e8,
        "mhalo_max_msun": 1.0e13,
    }

    shmr_model_name = str(shmr_block.get("model", "power_law")).lower()
    if shmr_model_name == "behroozi13":
        ensemble_config["shmr_model"] = {"type": "behroozi13", "pivot_log10_mstar": 7.5}
    else:
        ensemble_config["shmr_model"] = {
            "type": "power_law",
            "a": 10.5,
            "b": 1.4,
            "pivot_log10_mstar": 7.5,
        }

    concentration_name = str(concentration.get("model", "maccio")).lower()
    if concentration_name == "diemerjoyce19":
        ensemble_config["concentration_model"] = {"type": "diemerjoyce19"}
    else:
        ensemble_config["concentration_model"] = {"type": "maccio"}

    return ensemble_config


def load_ensemble_yaml_config(path: Path) -> EnsembleYamlConfig:
    """Parse and normalize one YAML ensemble configuration file."""
    with path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)

    mode = str(raw["mode"]).upper()
    label = path.stem.replace("_ensemble_config", "")

    base_config: dict[str, Any] = {
        "n_halos": int(raw["ensemble"]["N_halos"]),
        "seed": int(raw["ensemble"]["random_seed"]),
        **_build_redshift_configuration(dict(raw["redshift"])),
    }

    if mode == "HMF":
        ensemble_config = _normalize_hmf_mode(raw, base_config)
    elif mode == "SHMR":
        ensemble_config = _normalize_shmr_mode(raw, base_config)
    else:
        raise ValueError(f"Unsupported mode in YAML: {mode}")

    projection = dict(raw["projection"])
    projection_config = {
        "r_min_kpc": float(projection["R_min_kpc"]),
        "r_max_kpc": float(projection["R_max_kpc"]),
        "n_r_bins": int(projection["N_R_bins"]),
    }

    sidm_block = dict(raw["sidm"])
    sidm_parameterization = str(sidm_block.get("parameterization", "effective")).lower()
    if sidm_parameterization == "effective":
        sidm_sigma_grid = tuple(float(value) for value in sidm_block["sigma_over_m_grid"])
    elif sidm_parameterization == "velocity_dependent":
        sidm_sigma_grid = tuple(float(value) for value in sidm_block["sigma0_over_m_grid"])
    else:
        raise ValueError(
            "sidm.parameterization must be 'effective' or 'velocity_dependent', "
            f"received {sidm_parameterization!r}."
        )

    cdm_block = dict(sidm_block.get("cdm_reference", {}))
    sidm_config = {
        "parameterization": sidm_parameterization,
        "w_km_s": float(sidm_block.get("w_km_s", 0.0)),
        "time_model": str(sidm_block.get("time_model", "lookback_to_z")),
        "mass_definition": str(sidm_block.get("mass_definition", "200c")).lower(),
        "cdm_profile_source": str(cdm_block.get("profile_source", "nfw")).lower(),
        "cdm_sigma_over_m": float(cdm_block.get("sigma0_over_m", 0.0)),
        "cdm_w_km_s": float(cdm_block.get("w_km_s", 0.0)),
        "cdm_time_model": str(cdm_block.get("time_model", sidm_block.get("time_model", "lookback_to_z"))),
    }

    tier2_block = dict(raw.get("tier2", {}))
    tier2_config = {
        "enabled": bool(tier2_block.get("enabled", False)),
        "outer_profile_model": str(tier2_block.get("outer_profile_model", "dk14_like")),
        "stitch_method": str(tier2_block.get("stitch_method", "logistic_logrho_blend")),
        "r_match_mode": str(tier2_block.get("r_match_mode", "fraction_r200m")),
        "r_match_value": float(tier2_block.get("r_match_value", 0.8)),
        "smooth_width_dex": float(tier2_block.get("smooth_width_dex", 0.15)),
        "continuity": str(tier2_block.get("continuity", "density")),
        "regime": str(tier2_block.get("regime", "cluster" if mode == "HMF" else "dwarf")),
    }

    regime_overrides = dict(tier2_block.get("regime_overrides", {}))
    active_regime_override = regime_overrides.get(tier2_config["regime"], None)
    if isinstance(active_regime_override, dict):
        for key, value in active_regime_override.items():
            tier2_config[key] = value

    tier3_block = dict(raw.get("tier3", {}))
    tier3_config = {
        "enabled": bool(tier3_block.get("enabled", False)),
        "correction_model": str(tier3_block.get("correction_model", "rt_gamma_shift")),
        "sigma_pivot": float(tier3_block.get("sigma_pivot", 1.0)),
        "apply_to_regimes": tuple(
            str(value)
            for value in tier3_block.get("apply_to_regimes", ["cluster"])
        ),
        "calibration_mode": str(tier3_block.get("calibration_mode", "manual_preset")),
        "preset": str(tier3_block.get("preset", "none")),
        "regime": str(tier3_block.get("regime", "cluster" if mode == "HMF" else "dwarf")),
        "rt_shift": dict(tier3_block.get("rt_shift", {})),
        "gamma_shift": dict(tier3_block.get("gamma_shift", {})),
        "beta_shift": dict(tier3_block.get("beta_shift", {})),
        "outer_window": dict(tier3_block.get("outer_window", {})),
    }

    return EnsembleYamlConfig(
        label=label,
        mode=mode,
        ensemble_config=ensemble_config,
        projection_config=projection_config,
        sidm_sigma_grid=sidm_sigma_grid,
        sidm_config=sidm_config,
        tier2_config=tier2_config,
        tier3_config=tier3_config,
    )
