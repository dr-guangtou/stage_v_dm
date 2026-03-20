"""
Main forecast driver.

Runs the Fisher forecast for survey configurations defined in a YAML
config file (preferred) or from the hard-coded survey_configs module
(legacy mode). Saves results as JSON and NPZ, and generates summary
figures.

Usage:
    uv run python scripts/run_forecast.py --config configs/default.yaml
    uv run python scripts/run_forecast.py --config configs/stage4_vs_stage5.yaml
    uv run python scripts/run_forecast.py --surveys stage4_low_z stage5_wide
    uv run python scripts/run_forecast.py --systematics
"""

import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shmr_fisher.config import (
    ForecastConfig,
    LensingConfig,
    NuisanceConfig,
    SHMRParams,
    SurveyConfig,
)
from shmr_fisher.config_io import RunConfig, load_run_config
from shmr_fisher.fisher import (
    compute_fisher_matrix,
    extract_shmr_constraints,
    marginalized_errors,
)
from shmr_fisher.survey_configs import surveys as all_surveys


def run_forecast(
    survey_names: list[str] | None = None,
    survey_dict: dict[str, SurveyConfig] | None = None,
    forecast_config: ForecastConfig | None = None,
    lensing_config: LensingConfig | None = None,
    nuisance_config: NuisanceConfig | None = None,
    shmr_params: SHMRParams | None = None,
    output_dir: str = "outputs",
) -> dict:
    """
    Run the Fisher forecast for specified surveys.

    Parameters
    ----------
    survey_names : list of str or None
        Survey keys from survey_configs.surveys. Ignored if survey_dict
        is provided.
    survey_dict : dict of str to SurveyConfig or None
        Explicit survey configurations (from YAML). Takes priority over
        survey_names.
    forecast_config : ForecastConfig or None
        Uses defaults if None.
    lensing_config : LensingConfig or None
    nuisance_config : NuisanceConfig or None
    shmr_params : SHMRParams or None
    output_dir : str
        Directory for output files.

    Returns
    -------
    results : dict
        Keyed by survey name, containing fisher, errors, metadata.
    """
    if shmr_params is None:
        shmr_params = SHMRParams()
    if forecast_config is None:
        forecast_config = ForecastConfig()
    if lensing_config is None:
        lensing_config = LensingConfig()

    # Resolve active surveys: explicit dict > named list > all defaults
    if survey_dict is not None:
        active_surveys = survey_dict
    elif survey_names is not None:
        active_surveys = {}
        for name in survey_names:
            if name not in all_surveys:
                print(f"WARNING: survey '{name}' not found, skipping.")
                continue
            active_surveys[name] = all_surveys[name]
    else:
        active_surveys = dict(all_surveys)

    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    results = {}

    for name, survey in active_surveys.items():
        print(f"\n{'='*60}")
        print(f"Running forecast: {survey.name}")
        print(f"  Area: {survey.area_deg2} deg^2, z=[{survey.z_min}, {survey.z_max}]")
        print(f"  N_gal: {survey.n_gal_total:.0e}, log M*_min: {survey.log_Mstar_min}")
        print(f"  Systematic floor: {forecast_config.systematic_floor_fraction}")
        print(f"  Nuisance params: {forecast_config.include_nuisance_params}")
        print(f"{'='*60}")

        nuis = nuisance_config if forecast_config.include_nuisance_params else None

        fisher, param_names, metadata = compute_fisher_matrix(
            shmr_params, survey, lensing_config, forecast_config, nuis,
        )

        eigvals = np.linalg.eigvalsh(fisher)
        pos_def = bool(np.all(eigvals > 0))

        if pos_def:
            errors = marginalized_errors(fisher)
            # Also extract SHMR-only errors if nuisance params present
            shmr_errors, shmr_names = extract_shmr_constraints(fisher, param_names)
        else:
            errors = np.full(len(param_names), np.inf)
            shmr_errors = errors
            shmr_names = param_names
            print("  WARNING: Fisher matrix is not positive definite!")

        results[name] = {
            "survey_name": survey.name,
            "param_names": param_names,
            "shmr_param_names": shmr_names,
            "errors": dict(zip(param_names, errors.tolist())),
            "shmr_errors": dict(zip(shmr_names, shmr_errors.tolist())),
            "condition_number": metadata.get("condition_number", np.inf),
            "positive_definite": pos_def,
            "n_z_bins": len(metadata.get("z_bins", [])),
            "n_Mstar_bins_total": sum(
                len(v) for v in metadata.get("Mstar_bins_per_z", {}).values()
            ),
            "fisher": fisher,
            "metadata": metadata,
        }

        # Print summary
        print(f"\n  Condition number: {results[name]['condition_number']:.2e}")
        print(f"  Positive definite: {pos_def}")
        print(f"  z-bins: {results[name]['n_z_bins']}, "
              f"M*-bins total: {results[name]['n_Mstar_bins_total']}")
        print(f"\n  SHMR parameter constraints:")
        for pname, err in zip(shmr_names, shmr_errors):
            fid = abs(getattr(shmr_params, pname))
            frac = err / fid * 100 if fid > 1e-10 else 0
            print(f"    {pname:15s}: sigma = {err:.6f} ({frac:.3f}%)")

    # Save results
    _save_results(results, shmr_params, forecast_config, outdir)

    return results


def run_forecast_from_config(run_config: RunConfig) -> dict:
    """
    Run the Fisher forecast using a RunConfig loaded from YAML.

    This is the preferred entry point for YAML-driven runs. It also
    copies the config file into the output directory for reproducibility.

    Parameters
    ----------
    run_config : RunConfig
        Loaded configuration.

    Returns
    -------
    results : dict
        Keyed by survey name.
    """
    return run_forecast(
        survey_dict=run_config.surveys,
        forecast_config=run_config.forecast_config,
        lensing_config=run_config.lensing_config,
        nuisance_config=run_config.nuisance_config,
        shmr_params=run_config.shmr_params,
        output_dir=str(run_config.output_dir),
    )


def _save_results(
    results: dict,
    shmr_params: SHMRParams,
    forecast_config: ForecastConfig,
    outdir: Path,
) -> None:
    """Save forecast results as JSON and NPZ."""

    # JSON: human-readable summary (no numpy arrays)
    json_data = {
        "fiducial_params": shmr_params.to_dict(),
        "forecast_settings": {
            "systematic_floor_fraction": forecast_config.systematic_floor_fraction,
            "include_nuisance_params": forecast_config.include_nuisance_params,
            "n_R_bins": forecast_config.n_R_bins,
            "dz": forecast_config.dz,
            "frac_step": forecast_config.frac_step,
        },
        "surveys": {},
    }

    npz_data = {}

    for name, r in results.items():
        json_data["surveys"][name] = {
            "survey_name": r["survey_name"],
            "param_names": r["param_names"],
            "shmr_param_names": r["shmr_param_names"],
            "errors": r["errors"],
            "shmr_errors": r["shmr_errors"],
            "condition_number": float(r["condition_number"]),
            "positive_definite": r["positive_definite"],
            "n_z_bins": r["n_z_bins"],
            "n_Mstar_bins_total": r["n_Mstar_bins_total"],
        }

        npz_data[f"{name}_fisher"] = r["fisher"]
        npz_data[f"{name}_param_names"] = np.array(r["param_names"])

    json_path = outdir / "forecast_results.json"
    with open(json_path, "w") as f:
        json.dump(json_data, f, indent=2)
    print(f"\nSaved: {json_path}")

    npz_path = outdir / "forecast_results.npz"
    np.savez(npz_path, **npz_data)
    print(f"Saved: {npz_path}")


def main():
    parser = argparse.ArgumentParser(description="SHMR Fisher Forecast")
    parser.add_argument(
        "--config", type=str, default=None,
        help="Path to YAML config file (preferred over --surveys)",
    )
    parser.add_argument(
        "--surveys", nargs="+", default=None,
        help="Survey names to run (legacy mode, default: all)",
    )
    parser.add_argument(
        "--systematics", action="store_true",
        help="Enable both systematic floor (5%%) and nuisance params",
    )
    parser.add_argument(
        "--floor", type=float, default=0.0,
        help="Systematic floor fraction (default: 0.0)",
    )
    parser.add_argument(
        "--nuisance", action="store_true",
        help="Enable nuisance parameter marginalization",
    )
    parser.add_argument(
        "--sigma-m", type=float, default=0.02,
        help="Prior on shear calibration bias (default: 0.02)",
    )
    parser.add_argument(
        "--sigma-dz", type=float, default=0.03,
        help="Prior on source photo-z bias (default: 0.03)",
    )
    parser.add_argument(
        "--output-dir", default=None,
        help="Output directory override (default: from config or 'outputs')",
    )
    args = parser.parse_args()

    # YAML config mode (preferred)
    if args.config is not None:
        run_config = load_run_config(args.config)

        # Allow output-dir override from CLI
        if args.output_dir is not None:
            run_config.output_dir = Path(args.output_dir)

        # Copy config into output dir for reproducibility
        run_config.output_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(args.config, run_config.output_dir / "config.yaml")

        results = run_forecast_from_config(run_config)

    # Legacy CLI mode
    else:
        if args.systematics:
            floor = 0.05
            nuisance = True
        else:
            floor = args.floor
            nuisance = args.nuisance

        fc = ForecastConfig(
            n_R_bins=10,
            n_Mh_bins=200,
            systematic_floor_fraction=floor,
            include_nuisance_params=nuisance,
        )

        nc = NuisanceConfig(
            sigma_m=args.sigma_m,
            sigma_dz_source=args.sigma_dz,
        ) if nuisance else None

        output_dir = args.output_dir or "outputs"

        results = run_forecast(
            survey_names=args.surveys,
            forecast_config=fc,
            nuisance_config=nc,
            output_dir=output_dir,
        )

    print(f"\n{'='*60}")
    print(f"Forecast complete: {len(results)} surveys processed.")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
