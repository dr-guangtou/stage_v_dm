"""
YAML-based run configuration loader.

Parses a YAML config file into a RunConfig dataclass that drives the
entire forecast pipeline. This separates I/O and validation from the
pure dataclass definitions in config.py.

Example YAML structure::

    shmr_params:
      use_mass_dependent_scatter: true

    forecast:
      systematic_floor_fraction: 0.05
      include_nuisance_params: true

    lensing:
      n_source_per_arcmin2: 25.0

    nuisance:
      sigma_m: 0.02

    surveys:
      stage4_low_z:
        name: "Stage-IV Low-z"
        area_deg2: 14000
        z_min: 0.05
        z_max: 0.4
        n_gal_total: 10_000_000
        log_Mstar_min: 9.0

      stage5_wide:
        name: "Stage-V Wide"
        area_deg2: 10000
        z_min: 0.05
        z_max: 1.0
        n_gal_total: 50_000_000
        log_Mstar_min: 8.5
        log_Mstar_min_z_dep:
          base: 8.5
          slope: 1.0
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, fields
from pathlib import Path

import yaml

from .config import (
    ForecastConfig,
    LensingConfig,
    NuisanceConfig,
    SHMRParams,
    SurveyConfig,
)


@dataclass
class RunConfig:
    """
    Top-level container for a complete forecast run.

    Created by loading a YAML config file via ``load_run_config()``.

    Attributes
    ----------
    run_name : str
        Name of this run, derived from the YAML filename.
    shmr_params : SHMRParams
        SHMR model parameters.
    forecast_config : ForecastConfig
        Fisher forecast settings.
    lensing_config : LensingConfig
        Fixed lensing survey assumptions.
    nuisance_config : NuisanceConfig or None
        Nuisance parameter priors. None if not included.
    surveys : dict of str to SurveyConfig
        Survey configurations keyed by short name.
    output_dir : Path
        Output directory: ``outputs/{run_name}/``.
    """

    run_name: str
    shmr_params: SHMRParams
    forecast_config: ForecastConfig
    lensing_config: LensingConfig
    nuisance_config: NuisanceConfig | None
    surveys: dict[str, SurveyConfig]
    output_dir: Path


# ---------------------------------------------------------------------------
# Required fields for each survey entry in the YAML
# ---------------------------------------------------------------------------
_REQUIRED_SURVEY_FIELDS = {
    "name", "area_deg2", "z_min", "z_max", "n_gal_total", "log_Mstar_min",
}


def load_run_config(yaml_path: str | Path) -> RunConfig:
    """
    Load a YAML config file and return a fully validated RunConfig.

    Parameters
    ----------
    yaml_path : str or Path
        Path to the YAML configuration file.

    Returns
    -------
    run_config : RunConfig

    Raises
    ------
    FileNotFoundError
        If the YAML file does not exist.
    ValueError
        If the config is malformed (missing surveys, wrong number of
        surveys, missing required fields, etc.).
    """
    yaml_path = Path(yaml_path)
    if not yaml_path.exists():
        raise FileNotFoundError(f"Config file not found: {yaml_path}")

    with open(yaml_path) as f:
        raw = yaml.safe_load(f)

    if not isinstance(raw, dict):
        raise ValueError(f"Config file must be a YAML mapping, got {type(raw)}")

    # Run name: explicit override or filename stem
    run_name = raw.get("run_name", yaml_path.stem)

    # Build component configs
    shmr_params = _build_shmr_params(raw.get("shmr_params", {}))
    forecast_config = _build_forecast_config(raw.get("forecast", {}))
    lensing_config = _build_lensing_config(raw.get("lensing", {}))

    # Nuisance config: only build if forecast says to include it
    nuisance_config = None
    if forecast_config.include_nuisance_params:
        nuisance_config = _build_nuisance_config(raw.get("nuisance", {}))
        if not raw.get("nuisance"):
            warnings.warn(
                "include_nuisance_params is True but no 'nuisance' section "
                "in config — using NuisanceConfig defaults.",
                stacklevel=2,
            )

    # Surveys
    surveys_raw = raw.get("surveys", {})
    if not surveys_raw:
        raise ValueError("Config must define at least one survey in 'surveys' section.")
    if len(surveys_raw) > 4:
        raise ValueError(
            f"Config defines {len(surveys_raw)} surveys; maximum is 4."
        )

    surveys = {}
    for key, sdict in surveys_raw.items():
        surveys[key] = _build_survey_config(key, sdict)

    output_dir = Path("outputs") / run_name

    return RunConfig(
        run_name=run_name,
        shmr_params=shmr_params,
        forecast_config=forecast_config,
        lensing_config=lensing_config,
        nuisance_config=nuisance_config,
        surveys=surveys,
        output_dir=output_dir,
    )


# ---------------------------------------------------------------------------
# Internal builders
# ---------------------------------------------------------------------------

def _build_shmr_params(d: dict | None) -> SHMRParams:
    """Construct SHMRParams from a dict, filtering to valid fields."""
    if not d:
        return SHMRParams()
    valid_keys = {f.name for f in fields(SHMRParams)}
    filtered = {k: v for k, v in d.items() if k in valid_keys}
    return SHMRParams(**filtered)


def _build_forecast_config(d: dict | None) -> ForecastConfig:
    """Construct ForecastConfig from a dict, filtering to valid fields."""
    if not d:
        return ForecastConfig()
    valid_keys = {f.name for f in fields(ForecastConfig)}
    filtered = {k: v for k, v in d.items() if k in valid_keys}
    return ForecastConfig(**filtered)


def _build_lensing_config(d: dict | None) -> LensingConfig:
    """Construct LensingConfig from a dict, filtering to valid fields."""
    if not d:
        return LensingConfig()
    valid_keys = {f.name for f in fields(LensingConfig)}
    filtered = {k: v for k, v in d.items() if k in valid_keys}
    return LensingConfig(**filtered)


def _build_nuisance_config(d: dict | None) -> NuisanceConfig:
    """Construct NuisanceConfig from a dict, filtering to valid fields."""
    if not d:
        return NuisanceConfig()
    valid_keys = {f.name for f in fields(NuisanceConfig)}
    filtered = {k: v for k, v in d.items() if k in valid_keys}
    return NuisanceConfig(**filtered)


def _build_survey_config(key: str, d: dict) -> SurveyConfig:
    """
    Construct a SurveyConfig from a YAML survey dict.

    Handles the special ``log_Mstar_min_z_dep`` field by constructing
    a lambda from its ``base`` and ``slope`` values.

    Parameters
    ----------
    key : str
        Survey key name (for error messages).
    d : dict
        Raw YAML dict for this survey.

    Returns
    -------
    survey : SurveyConfig

    Raises
    ------
    ValueError
        If required fields are missing.
    """
    missing = _REQUIRED_SURVEY_FIELDS - set(d.keys())
    if missing:
        raise ValueError(
            f"Survey '{key}' is missing required fields: {sorted(missing)}"
        )

    # Extract z-dependent completeness if present
    z_dep = d.pop("log_Mstar_min_z_dep", None)
    log_Mstar_min_func = None
    if z_dep is not None:
        base = z_dep["base"]
        slope = z_dep["slope"]
        # Use default args to capture values in closure
        log_Mstar_min_func = lambda z, _b=base, _s=slope: _b + _s * z

    # Filter to valid SurveyConfig fields
    valid_keys = {f.name for f in fields(SurveyConfig)}
    filtered = {k: v for k, v in d.items() if k in valid_keys}

    if log_Mstar_min_func is not None:
        filtered["log_Mstar_min_func"] = log_Mstar_min_func

    return SurveyConfig(**filtered)
