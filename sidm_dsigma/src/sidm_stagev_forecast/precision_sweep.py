"""Precision sweep utilities for target-significance SIDM separation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from sidm_stagev_forecast.config import (
    BENCHMARKS,
    DEFAULT_COSMOLOGY,
    DEFAULT_FORECAST_CONFIG,
    BenchmarkConfig,
    log_radius_grid,
)
from sidm_stagev_forecast.forecast import required_uniform_fractional_precision
from sidm_stagev_forecast.profiles import nfw_profile_bundle, sidm_profile_from_parametric_model
from sidm_stagev_forecast.projection import delta_sigma_of_R


@dataclass(frozen=True)
class RadialWindow:
    """Named radial window in projected radius space."""

    name: str
    r_min_kpc: float
    r_max_kpc: float


def default_radial_windows(benchmark: BenchmarkConfig) -> tuple[RadialWindow, ...]:
    """Define standard full/inner/outer windows for a benchmark."""
    r_min = benchmark.r_lensing_min_kpc
    r_max = benchmark.r_lensing_max_kpc
    split = np.sqrt(r_min * r_max)
    return (
        RadialWindow("full", r_min, r_max),
        RadialWindow("inner", r_min, split),
        RadialWindow("outer", split, r_max),
    )


def _window_mask(r_projected_kpc: np.ndarray, window: RadialWindow) -> np.ndarray:
    mask = (r_projected_kpc >= window.r_min_kpc) & (r_projected_kpc <= window.r_max_kpc)
    if not np.any(mask):
        raise ValueError(f"No points in window '{window.name}'")
    return mask


def _resolve_sidm_model_options(sidm_backend: str):
    if sidm_backend == "parametric":
        return None
    if sidm_backend == "surrogate":

        def surrogate_sidm_callable(
            r_kpc: np.ndarray,
            m200_msun: float,
            c200: float,
            z: float,
            sigma_over_m: float,
        ) -> np.ndarray:
            cdm = nfw_profile_bundle(r_kpc, m200_msun, c200, z, DEFAULT_COSMOLOGY)
            response = sigma_over_m / (1.0 + sigma_over_m)
            r200_guess = np.max(r_kpc)
            r_core_kpc = 0.02 * r200_guess * response
            return cdm["rho_msun_kpc3"] * (r_kpc / (r_kpc + r_core_kpc)) ** 2

        return {"sidm_callable": surrogate_sidm_callable}
    raise ValueError("sidm_backend must be 'parametric' or 'surrogate'.")


def run_precision_sweep(
    target_significance_levels: tuple[float, ...] = (2.0, 3.0, 5.0),
    sidm_backend: str = "parametric",
) -> pd.DataFrame:
    """Compute required precision for multiple significance targets and windows."""
    summary_rows: list[dict[str, float | str | int]] = []

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

        windows = default_radial_windows(benchmark)

        for sigma_over_m in benchmark.sigma_over_m_grid_cm2_g:
            if np.isclose(sigma_over_m, 0.0):
                continue

            sidm_profile = sidm_profile_from_parametric_model(
                r_kpc=r_3d_kpc,
                m200_msun=benchmark.m200_msun,
                c200=benchmark.c200,
                z=benchmark.redshift,
                sigma_over_m=sigma_over_m,
                model_options=_resolve_sidm_model_options(sidm_backend),
            )
            sidm_projection = delta_sigma_of_R(
                r_projected_kpc=r_lensing_kpc,
                r_kpc=r_3d_kpc,
                rho_msun_kpc3=sidm_profile["rho_msun_kpc3"],
                n_z=900,
            )

            for window in windows:
                mask = _window_mask(r_lensing_kpc, window)
                reference_slice = cdm_projection["delta_sigma_msun_kpc2"][mask]
                model_slice = sidm_projection["delta_sigma_msun_kpc2"][mask]

                max_shift_percent = float(
                    100.0 * np.max(np.abs(model_slice / np.maximum(reference_slice, 1.0e-12) - 1.0))
                )

                for target_sigma in target_significance_levels:
                    required_fraction = required_uniform_fractional_precision(
                        model=model_slice,
                        reference=reference_slice,
                        target_sigma_significance=target_sigma,
                    )
                    summary_rows.append(
                        {
                            "benchmark": benchmark.label,
                            "m200_msun": benchmark.m200_msun,
                            "sigma_over_m_cm2_g": sigma_over_m,
                            "window": window.name,
                            "window_r_min_kpc": window.r_min_kpc,
                            "window_r_max_kpc": window.r_max_kpc,
                            "n_bins": int(np.sum(mask)),
                            "target_sigma": target_sigma,
                            "required_uniform_fraction": required_fraction,
                            "required_uniform_percent": 100.0 * required_fraction,
                            "max_abs_ratio_shift_percent": max_shift_percent,
                        }
                    )

    return pd.DataFrame(summary_rows)
