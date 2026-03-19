"""Run benchmark SIDM-vs-CDM forecasts and save figures/tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from sidm_stagev_forecast.config import (
    BENCHMARKS,
    DEFAULT_COSMOLOGY,
    DEFAULT_FORECAST_CONFIG,
    log_radius_grid,
)
from sidm_stagev_forecast.forecast import (
    delta_chi2,
    fractional_error_model,
    required_uniform_fractional_precision,
)
from sidm_stagev_forecast.io import ensure_output_directories, save_table
from sidm_stagev_forecast.plotting import plot_delta_sigma_profiles, plot_rho_profiles
from sidm_stagev_forecast.profiles import (
    nfw_parameters_from_m_c,
    nfw_profile_bundle,
    sidm_profile_from_parametric_model,
)
from sidm_stagev_forecast.projection import delta_sigma_of_R


def surrogate_sidm_callable(
    r_kpc: np.ndarray,
    m200_msun: float,
    c200: float,
    z: float,
    sigma_over_m: float,
) -> np.ndarray:
    """Simple cored-NFW surrogate used only when explicit fallback is requested."""
    cdm = nfw_profile_bundle(r_kpc, m200_msun, c200, z, DEFAULT_COSMOLOGY)
    rho_cdm = cdm["rho_msun_kpc3"]
    nfw_parameters = nfw_parameters_from_m_c(m200_msun, c200, z, DEFAULT_COSMOLOGY)

    response = sigma_over_m / (1.0 + sigma_over_m)
    r_core_kpc = 0.5 * nfw_parameters.rs_kpc * response
    core_suppression = (r_kpc / (r_kpc + r_core_kpc)) ** 2
    return rho_cdm * core_suppression


def run_benchmarks(output_root: Path, sidm_backend: str) -> None:
    """Execute full dwarf+cluster benchmark suite."""
    output_paths = ensure_output_directories(output_root)

    summary_rows: list[dict[str, float | str]] = []

    for benchmark in BENCHMARKS:
        r_3d_kpc = log_radius_grid(
            benchmark.r_3d_min_kpc,
            benchmark.r_3d_max_kpc,
            DEFAULT_FORECAST_CONFIG.n_r_3d,
        )
        r_lensing_kpc = log_radius_grid(
            benchmark.r_lensing_min_kpc,
            benchmark.r_lensing_max_kpc,
            DEFAULT_FORECAST_CONFIG.n_r_lensing,
        )

        cdm_profile = nfw_profile_bundle(
            r_kpc=r_3d_kpc,
            m200_msun=benchmark.m200_msun,
            c200=benchmark.c200,
            z=benchmark.redshift,
            cosmo=DEFAULT_COSMOLOGY,
        )

        cdm_projection = delta_sigma_of_R(
            r_projected_kpc=r_lensing_kpc,
            r_kpc=r_3d_kpc,
            rho_msun_kpc3=cdm_profile["rho_msun_kpc3"],
            n_z=900,
        )

        sidm_rho_profiles: dict[float, np.ndarray] = {}
        sidm_delta_sigma_profiles: dict[float, np.ndarray] = {}

        fractional_errors = fractional_error_model(
            r_lensing_kpc, benchmark.label, scenario="baseline"
        )

        for sigma_over_m in benchmark.sigma_over_m_grid_cm2_g:
            if np.isclose(sigma_over_m, 0.0):
                sidm_profile = cdm_profile
            else:
                model_options = None
                if sidm_backend == "surrogate":
                    model_options = {"sidm_callable": surrogate_sidm_callable}

                sidm_profile = sidm_profile_from_parametric_model(
                    r_kpc=r_3d_kpc,
                    m200_msun=benchmark.m200_msun,
                    c200=benchmark.c200,
                    z=benchmark.redshift,
                    sigma_over_m=sigma_over_m,
                    model_options=model_options,
                )

            sidm_projection = delta_sigma_of_R(
                r_projected_kpc=r_lensing_kpc,
                r_kpc=r_3d_kpc,
                rho_msun_kpc3=sidm_profile["rho_msun_kpc3"],
                n_z=900,
            )

            sidm_rho_profiles[sigma_over_m] = sidm_profile["rho_msun_kpc3"]
            sidm_delta_sigma_profiles[sigma_over_m] = sidm_projection["delta_sigma_msun_kpc2"]

            sigma_abs = fractional_errors * np.abs(cdm_projection["delta_sigma_msun_kpc2"])
            sigma_abs = np.maximum(sigma_abs, 1.0e-12)
            metric_delta_chi2 = delta_chi2(
                sidm_projection["delta_sigma_msun_kpc2"],
                cdm_projection["delta_sigma_msun_kpc2"],
                sigma_abs,
            )

            required_fraction = required_uniform_fractional_precision(
                sidm_projection["delta_sigma_msun_kpc2"],
                cdm_projection["delta_sigma_msun_kpc2"],
                target_sigma_significance=3.0,
            )

            max_ratio_shift_percent = 100.0 * np.max(
                np.abs(
                    sidm_projection["delta_sigma_msun_kpc2"]
                    / cdm_projection["delta_sigma_msun_kpc2"]
                    - 1.0
                )
            )

            summary_rows.append(
                {
                    "benchmark": benchmark.label,
                    "m200_msun": benchmark.m200_msun,
                    "c200": benchmark.c200,
                    "z": benchmark.redshift,
                    "sigma_over_m_cm2_g": sigma_over_m,
                    "delta_chi2_baseline": metric_delta_chi2,
                    "snr_baseline": np.sqrt(metric_delta_chi2),
                    "required_uniform_fraction_for_3sigma": required_fraction,
                    "max_abs_delta_sigma_ratio_shift_percent": max_ratio_shift_percent,
                }
            )

            intermediate_profile_table = pd.DataFrame(
                {
                    "r_kpc": r_3d_kpc,
                    "rho_cdm_msun_kpc3": cdm_profile["rho_msun_kpc3"],
                    "rho_sidm_msun_kpc3": sidm_profile["rho_msun_kpc3"],
                }
            )
            intermediate_profile_name = (
                f"{benchmark.label}_profile_sigma_{sigma_over_m:.3f}".replace(".", "p") + ".csv"
            )
            save_table(
                intermediate_profile_table, output_paths["intermediate"] / intermediate_profile_name
            )

            intermediate_delta_sigma_table = pd.DataFrame(
                {
                    "r_projected_kpc": r_lensing_kpc,
                    "delta_sigma_cdm_msun_kpc2": cdm_projection["delta_sigma_msun_kpc2"],
                    "delta_sigma_sidm_msun_kpc2": sidm_projection["delta_sigma_msun_kpc2"],
                }
            )
            intermediate_delta_sigma_name = (
                f"{benchmark.label}_delta_sigma_sigma_{sigma_over_m:.3f}".replace(".", "p") + ".csv"
            )
            save_table(
                intermediate_delta_sigma_table,
                output_paths["intermediate"] / intermediate_delta_sigma_name,
            )

        rho_figure_path = output_paths["figures"] / f"{benchmark.label}_rho_profiles.png"
        plot_rho_profiles(
            r_kpc=r_3d_kpc,
            rho_cdm_msun_kpc3=cdm_profile["rho_msun_kpc3"],
            sidm_profiles=sidm_rho_profiles,
            output_path=rho_figure_path,
            title=benchmark.label.capitalize(),
        )

        delta_sigma_figure_path = (
            output_paths["figures"] / f"{benchmark.label}_delta_sigma_profiles.png"
        )
        plot_delta_sigma_profiles(
            r_projected_kpc=r_lensing_kpc,
            delta_sigma_cdm_msun_kpc2=cdm_projection["delta_sigma_msun_kpc2"],
            sidm_delta_sigma_profiles=sidm_delta_sigma_profiles,
            fractional_error=fractional_errors,
            output_path=delta_sigma_figure_path,
            title=benchmark.label.capitalize(),
        )

    summary_table = pd.DataFrame(summary_rows)
    save_table(summary_table, output_paths["tables"] / "benchmark_summary.csv")


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("outputs"),
        help="Output directory root (default: outputs)",
    )
    parser.add_argument(
        "--sidm-backend",
        choices=["parametric", "surrogate"],
        default="parametric",
        help="Use parametricSIDM backend or a local surrogate fallback.",
    )
    return parser.parse_args()


def main() -> None:
    """Program entrypoint."""
    arguments = parse_arguments()
    run_benchmarks(arguments.output_root, arguments.sidm_backend)


if __name__ == "__main__":
    main()
