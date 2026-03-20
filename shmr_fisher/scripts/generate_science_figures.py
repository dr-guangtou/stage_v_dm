"""
Generate all science figures for a forecast run.

Reads survey configuration from a YAML config file and produces:
    - DeltaSigma with error bars
    - Two-regime summary (z=0 shape + evolution)
    - Improvement factor (first vs last survey)
    - Parameter sweeps (area, log_Mstar_min, z_max)

Usage:
    uv run python scripts/generate_science_figures.py --config configs/default.yaml
    uv run python scripts/generate_science_figures.py  # uses configs/default.yaml
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np

from shmr_fisher.config import ForecastConfig, NuisanceConfig, SHMRParams
from shmr_fisher.config_io import load_run_config
from shmr_fisher.plot_results import (
    plot_delta_sigma_with_errors,
    plot_two_regime_summary,
    plot_improvement_factor,
)
from scripts.run_forecast import run_forecast, run_forecast_from_config
from scripts.run_sweep import parameter_sweep, plot_sweep


def main():
    parser = argparse.ArgumentParser(
        description="Generate science figures for an SHMR forecast run"
    )
    parser.add_argument(
        "--config", type=str, default="configs/default.yaml",
        help="Path to YAML config file (default: configs/default.yaml)",
    )
    parser.add_argument(
        "--skip-sweeps", action="store_true",
        help="Skip parameter sweep figures (faster)",
    )
    args = parser.parse_args()

    # Load configuration
    run_config = load_run_config(args.config)
    outdir = run_config.output_dir
    outdir.mkdir(parents=True, exist_ok=True)

    survey_keys = list(run_config.surveys.keys())
    n_surveys = len(survey_keys)

    print("=" * 60)
    print(f"Run: {run_config.run_name}")
    print(f"Surveys: {survey_keys}")
    print(f"Output: {outdir}")
    print("=" * 60)

    # -------------------------------------------------------------------
    # Run the forecast
    # -------------------------------------------------------------------
    print("\nRunning forecast...")
    results = run_forecast_from_config(run_config)

    # -------------------------------------------------------------------
    # Figure: DeltaSigma with error bars
    # -------------------------------------------------------------------
    print("\n--- Generating DeltaSigma with error bars ---")
    plot_delta_sigma_with_errors(
        results, run_config.shmr_params,
        z=0.3, mstar_bin=(10.5, 11.0),
        save_path=outdir / "delta_sigma_with_errors.png",
    )

    # -------------------------------------------------------------------
    # Figure: Two-regime summary
    # -------------------------------------------------------------------
    print("\n--- Generating two-regime summary ---")
    plot_two_regime_summary(
        results, run_config.shmr_params,
        save_path=outdir / "two_regime_summary.png",
    )

    # -------------------------------------------------------------------
    # Figure: Improvement factor (first vs last survey)
    # -------------------------------------------------------------------
    if n_surveys >= 2:
        print("\n--- Generating improvement factor ---")
        baseline_key = survey_keys[0]
        target_key = survey_keys[-1]

        # Use 5-param Fisher for fair comparison
        from dataclasses import replace
        fc_5p = replace(
            run_config.forecast_config,
            n_R_bins=8, n_Mh_bins=150,
            vary_z_evolution=False,
        )
        results_5p = run_forecast(
            survey_dict={
                baseline_key: run_config.surveys[baseline_key],
                target_key: run_config.surveys[target_key],
            },
            forecast_config=fc_5p,
            nuisance_config=run_config.nuisance_config,
            shmr_params=run_config.shmr_params,
            output_dir=str(outdir),
        )
        plot_improvement_factor(
            results_5p[baseline_key], results_5p[target_key],
            baseline_label=results[baseline_key]["survey_name"],
            target_label=results[target_key]["survey_name"],
            save_path=outdir / "improvement_factor.png",
        )

    # -------------------------------------------------------------------
    # Parameter sweeps
    # -------------------------------------------------------------------
    if not args.skip_sweeps and n_surveys >= 1:
        fc_sweep = ForecastConfig(n_R_bins=8, n_Mh_bins=150)

        # Pick a representative survey for sweeps (first one)
        sweep_key = survey_keys[0]
        sweep_survey = run_config.surveys[sweep_key]

        print(f"\n--- Sweep: area_deg2 (base: {sweep_key}) ---")
        r_area = parameter_sweep(
            sweep_key, "area_deg2",
            [1000, 3000, 7000, 10000, 14000],
            forecast_config=fc_sweep,
            shmr_params=run_config.shmr_params,
            base_survey=sweep_survey,
        )
        plot_sweep(r_area, output_dir=str(outdir))

        print(f"\n--- Sweep: log_Mstar_min (base: {sweep_key}) ---")
        r_mmin = parameter_sweep(
            sweep_key, "log_Mstar_min",
            [9.5, 9.0, 8.5, 8.0, 7.5],
            forecast_config=fc_sweep,
            shmr_params=run_config.shmr_params,
            base_survey=sweep_survey,
        )
        plot_sweep(r_mmin, output_dir=str(outdir))

        print(f"\n--- Sweep: z_max (base: {sweep_key}) ---")
        r_zmax = parameter_sweep(
            sweep_key, "z_max",
            [0.2, 0.3, 0.4, 0.6, 0.8, 1.0],
            forecast_config=fc_sweep,
            shmr_params=run_config.shmr_params,
            base_survey=sweep_survey,
        )
        plot_sweep(r_zmax, output_dir=str(outdir))

    print("\n" + "=" * 60)
    print(f"All science figures saved to {outdir}")
    print("=" * 60)


if __name__ == "__main__":
    main()
