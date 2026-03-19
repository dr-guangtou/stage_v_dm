"""Run configurable Tier-1 halo-ensemble stacked DeltaSigma forecasts."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from sidm_stagev_forecast.config import (
    CLUSTER_HMF_ENSEMBLE_CONFIG_EXAMPLE,
    DEFAULT_COSMOLOGY,
    DEFAULT_ENSEMBLE_BENCHMARK,
    DEFAULT_TIER2_CONFIG,
    DWARF_SHMR_ENSEMBLE_CONFIG_EXAMPLE,
    log_radius_grid,
)
from sidm_stagev_forecast.cosmology import rdelta
from sidm_stagev_forecast.ensemble import generate_ensemble, summarize_ensemble
from sidm_stagev_forecast.ensemble_yaml import load_ensemble_yaml_config
from sidm_stagev_forecast.forecast import (
    evaluate_stacked_distinguishability,
    required_uniform_fractional_precision,
)
from sidm_stagev_forecast.io import ensure_output_directories, save_table
from sidm_stagev_forecast.io import append_figure_caption_entries, append_inventory_entries
from sidm_stagev_forecast.plotting import (
    plot_delta_chi2_tier_comparison,
    plot_single_halo_tier2_diagnostics,
    plot_stacked_delta_sigma,
    plot_tier1_summary,
    plot_tier1_tier2_stacked_comparison,
)
from sidm_stagev_forecast.profiles import (
    build_hybrid_sidm_profile,
    nfw_profile_bundle,
    sidm_profile_from_parametric_model,
)
from sidm_stagev_forecast.projection import delta_sigma_of_R
from sidm_stagev_forecast.stacking import (
    interpolate_profile_to_common_grid,
    stack_delta_sigma_profiles,
    weighted_stack_profiles,
)


def surrogate_sidm_callable(
    r_kpc: np.ndarray,
    m200_msun: float,
    c200: float,
    z: float,
    sigma_over_m: float,
) -> np.ndarray:
    """Simple cored-NFW fallback model for deterministic local runs."""
    cdm = nfw_profile_bundle(r_kpc, m200_msun, c200, z, DEFAULT_COSMOLOGY)
    response = sigma_over_m / (1.0 + sigma_over_m)
    r200_kpc = rdelta(m200_msun, z, DEFAULT_COSMOLOGY, definition="200c")
    r_core_kpc = 0.05 * r200_kpc * response
    return cdm["rho_msun_kpc3"] * (r_kpc / (r_kpc + r_core_kpc)) ** 2


def halo_specific_r_grids(
    m200_msun: float,
    z: float,
    n_r_3d: int,
    n_r_lensing_per_halo: int,
    r_common_min_kpc: float,
    r_common_max_kpc: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Build per-halo 3D and projected radius grids."""
    r200_kpc = rdelta(m200_msun, z, DEFAULT_COSMOLOGY, definition="200c")

    r_3d_min_kpc = max(0.01 * r200_kpc, 0.1)
    r_3d_max_kpc = max(5.0 * r200_kpc, r_common_max_kpc)
    r_3d_kpc = np.geomspace(r_3d_min_kpc, r_3d_max_kpc, n_r_3d)

    r_lensing_min_kpc = max(0.03 * r200_kpc, 0.6 * r_common_min_kpc)
    r_lensing_max_kpc = min(2.5 * r200_kpc, 1.25 * r_common_max_kpc)
    if r_lensing_max_kpc <= r_lensing_min_kpc:
        r_lensing_min_kpc = r_common_min_kpc
        r_lensing_max_kpc = r_common_max_kpc
    r_lensing_kpc = np.geomspace(r_lensing_min_kpc, r_lensing_max_kpc, n_r_lensing_per_halo)

    return r_3d_kpc, r_lensing_kpc


def _resolve_default_sampling_configuration(
    ensemble_mode: str,
    n_halos: int | None,
    seed: int | None,
) -> dict[str, Any]:
    if ensemble_mode == "HMF":
        configuration = dict(CLUSTER_HMF_ENSEMBLE_CONFIG_EXAMPLE)
    elif ensemble_mode == "SHMR":
        configuration = dict(DWARF_SHMR_ENSEMBLE_CONFIG_EXAMPLE)
    else:
        raise ValueError("ensemble_mode must be 'HMF' or 'SHMR'.")
    configuration["n_halos"] = int(
        DEFAULT_ENSEMBLE_BENCHMARK.n_halos if n_halos is None else n_halos
    )
    configuration["seed"] = int(DEFAULT_ENSEMBLE_BENCHMARK.seed if seed is None else seed)
    return configuration


def _resolve_runtime_configuration(
    ensemble_mode: str,
    n_halos: int | None,
    seed: int | None,
    config_path: Path | None,
) -> dict[str, Any]:
    if config_path is None:
        return {
            "label": ensemble_mode.lower(),
            "mode": ensemble_mode,
            "sampling_configuration": _resolve_default_sampling_configuration(
                ensemble_mode=ensemble_mode,
                n_halos=n_halos,
                seed=seed,
            ),
            "projection_configuration": None,
            "sigma_grid": tuple(DEFAULT_ENSEMBLE_BENCHMARK.sigma_over_m_grid_cm2_g),
            "tier2_configuration": {
                **DEFAULT_TIER2_CONFIG.__dict__,
                "regime": "cluster" if ensemble_mode == "HMF" else "dwarf",
            },
        }

    parsed = load_ensemble_yaml_config(config_path)
    sampling_configuration = dict(parsed.ensemble_config)
    if n_halos is not None:
        sampling_configuration["n_halos"] = int(n_halos)
    if seed is not None:
        sampling_configuration["seed"] = int(seed)
    return {
        "label": parsed.label,
        "mode": parsed.mode,
        "sampling_configuration": sampling_configuration,
        "projection_configuration": dict(parsed.projection_config),
        "sigma_grid": tuple(parsed.sidm_sigma_grid),
        "tier2_configuration": dict(parsed.tier2_config),
    }


def _resolve_common_radii(
    ensemble_mode: str,
    benchmark_n_common: int,
    projection_override: dict[str, float] | None,
) -> tuple[np.ndarray, np.ndarray]:
    if projection_override is None:
        if ensemble_mode == "HMF":
            r_common_lensing_kpc = log_radius_grid(30.0, 3000.0, benchmark_n_common)
        elif ensemble_mode == "SHMR":
            r_common_lensing_kpc = log_radius_grid(3.0, 300.0, benchmark_n_common)
        else:
            raise ValueError("ensemble_mode must be 'HMF' or 'SHMR'.")
    else:
        r_common_lensing_kpc = log_radius_grid(
            projection_override["r_min_kpc"],
            projection_override["r_max_kpc"],
            int(projection_override["n_r_bins"]),
        )

    if ensemble_mode == "HMF":
        r_common_3d_kpc = np.geomspace(2.0, 8.0e3, 240)
    elif ensemble_mode == "SHMR":
        r_common_3d_kpc = np.geomspace(0.1, 1.0e3, 240)
    else:
        raise ValueError("ensemble_mode must be 'HMF' or 'SHMR'.")
    return r_common_lensing_kpc, r_common_3d_kpc


def _validation_table(
    label: str,
    ensemble_mode: str,
    halos: list[dict[str, float]],
    r_common_lensing_kpc: np.ndarray,
    expected_projection_bins: int,
    sampling_configuration: dict[str, Any],
) -> pd.DataFrame:
    masses = np.asarray([halo["m200_msun"] for halo in halos], dtype=float)
    median_mass = float(np.median(masses))

    if ensemble_mode == "HMF":
        expected_mass = 3.0e14
        low, high = 1.5e14, 6.0e14
    else:
        expected_mass = 1.0e10
        low, high = 5.0e9, 2.0e10

    median_mass_check_pass = low <= median_mass <= high
    if not median_mass_check_pass:
        raise RuntimeError(
            f"{label}: median mass {median_mass:.3e} outside expected range "
            f"[{low:.3e}, {high:.3e}] for mode {ensemble_mode}."
        )

    radial_bin_check_pass = int(len(r_common_lensing_kpc)) == int(expected_projection_bins)
    if not radial_bin_check_pass:
        raise RuntimeError(
            f"{label}: expected {expected_projection_bins} radial bins, "
            f"got {len(r_common_lensing_kpc)}."
        )

    reproduce_first = generate_ensemble(mode=ensemble_mode, config_dict=sampling_configuration)
    reproduce_second = generate_ensemble(mode=ensemble_mode, config_dict=sampling_configuration)
    reproducibility_check_pass = reproduce_first == reproduce_second
    if not reproducibility_check_pass:
        raise RuntimeError(f"{label}: reproducibility check failed.")

    single_halo_config = dict(sampling_configuration)
    single_halo_config["n_halos"] = 1
    single_halo = generate_ensemble(mode=ensemble_mode, config_dict=single_halo_config)
    single_halo_limit_check_pass = len(single_halo) == 1 and np.isclose(
        single_halo[0]["weight"],
        1.0,
    )
    if not single_halo_limit_check_pass:
        raise RuntimeError(f"{label}: single-halo limit check failed.")

    return pd.DataFrame(
        [
            {
                "label": label,
                "mode": ensemble_mode,
                "median_mass_msun": median_mass,
                "expected_median_mass_msun": expected_mass,
                "median_mass_check_pass": bool(median_mass_check_pass),
                "expected_projection_bins": int(expected_projection_bins),
                "actual_projection_bins": int(len(r_common_lensing_kpc)),
                "radial_bin_check_pass": bool(radial_bin_check_pass),
                "reproducibility_check_pass": bool(reproducibility_check_pass),
                "single_halo_limit_check_pass": bool(single_halo_limit_check_pass),
            }
        ]
    )


def run_ensemble_forecast(
    output_root: Path,
    sidm_backend: str,
    n_halos: int | None,
    seed: int | None,
    ensemble_mode: str,
    config_path: Path | None = None,
) -> None:
    """Execute configurable ensemble stacking forecast and save outputs."""
    benchmark = DEFAULT_ENSEMBLE_BENCHMARK
    output_paths = ensure_output_directories(output_root)

    runtime_configuration = _resolve_runtime_configuration(
        ensemble_mode=ensemble_mode,
        n_halos=n_halos,
        seed=seed,
        config_path=config_path,
    )

    mode_tag = str(runtime_configuration["label"])
    active_mode = str(runtime_configuration["mode"])
    sampling_configuration = dict(runtime_configuration["sampling_configuration"])
    projection_configuration = runtime_configuration["projection_configuration"]
    sigma_grid = tuple(float(value) for value in runtime_configuration["sigma_grid"])
    tier2_configuration = dict(runtime_configuration["tier2_configuration"])
    tier2_enabled = bool(tier2_configuration.get("enabled", False))

    halos = generate_ensemble(mode=active_mode, config_dict=sampling_configuration)
    halo_summary = summarize_ensemble(halos)

    halo_catalog_table = pd.DataFrame(halos)
    save_table(
        halo_catalog_table,
        output_paths["intermediate"] / f"{mode_tag}_ensemble_halo_catalog.csv",
    )

    r_common_lensing_kpc, r_common_3d_kpc = _resolve_common_radii(
        ensemble_mode=active_mode,
        benchmark_n_common=benchmark.n_r_common,
        projection_override=projection_configuration,
    )
    r_common_min_kpc = float(r_common_lensing_kpc[0])
    r_common_max_kpc = float(r_common_lensing_kpc[-1])

    weights = np.asarray([halo["weight"] for halo in halos])
    cdm_projected_profiles: list[dict[str, np.ndarray]] = []
    sidm_projected_profiles: dict[float, list[dict[str, np.ndarray]]] = {
        sigma: [] for sigma in sigma_grid
    }
    cdm_rho_rows: list[np.ndarray] = []
    sidm_rho_rows: dict[float, list[np.ndarray]] = {sigma: [] for sigma in sigma_grid}
    cdm_projected_profiles_tier2: list[dict[str, np.ndarray]] = []
    sidm_projected_profiles_tier2: dict[float, list[dict[str, np.ndarray]]] = {
        sigma: [] for sigma in sigma_grid
    }
    cdm_rho_rows_tier2: list[np.ndarray] = []
    sidm_rho_rows_tier2: dict[float, list[np.ndarray]] = {sigma: [] for sigma in sigma_grid}

    model_options = None
    if sidm_backend == "surrogate":
        model_options = {"sidm_callable": surrogate_sidm_callable}

    single_halo_records: dict[str, dict[str, np.ndarray]] = {}

    for halo_index, halo in enumerate(halos):
        r_3d_kpc, r_lensing_kpc = halo_specific_r_grids(
            m200_msun=halo["m200_msun"],
            z=halo["z"],
            n_r_3d=benchmark.n_r_3d,
            n_r_lensing_per_halo=benchmark.n_r_lensing_per_halo,
            r_common_min_kpc=r_common_min_kpc,
            r_common_max_kpc=r_common_max_kpc,
        )

        cdm_profile = nfw_profile_bundle(
            r_kpc=r_3d_kpc,
            m200_msun=halo["m200_msun"],
            c200=halo["c200"],
            z=halo["z"],
            cosmo=DEFAULT_COSMOLOGY,
        )
        cdm_projected = delta_sigma_of_R(
            r_projected_kpc=r_lensing_kpc,
            r_kpc=r_3d_kpc,
            rho_msun_kpc3=cdm_profile["rho_msun_kpc3"],
            n_z=benchmark.projection_n_z,
        )
        cdm_projected_profiles.append(cdm_projected)
        cdm_rho_rows.append(
            interpolate_profile_to_common_grid(
                r_input_kpc=r_3d_kpc,
                profile_input=cdm_profile["rho_msun_kpc3"],
                r_common_kpc=r_common_3d_kpc,
            )
        )

        for sigma_over_m in sigma_grid:
            sidm_profile = sidm_profile_from_parametric_model(
                r_kpc=r_3d_kpc,
                m200_msun=halo["m200_msun"],
                c200=halo["c200"],
                z=halo["z"],
                sigma_over_m=sigma_over_m,
                model_options=model_options,
            )
            sidm_projected = delta_sigma_of_R(
                r_projected_kpc=r_lensing_kpc,
                r_kpc=r_3d_kpc,
                rho_msun_kpc3=sidm_profile["rho_msun_kpc3"],
                n_z=benchmark.projection_n_z,
            )
            sidm_projected_profiles[sigma_over_m].append(sidm_projected)
            sidm_rho_rows[sigma_over_m].append(
                interpolate_profile_to_common_grid(
                    r_input_kpc=r_3d_kpc,
                    profile_input=sidm_profile["rho_msun_kpc3"],
                    r_common_kpc=r_common_3d_kpc,
                )
            )

            if tier2_enabled:
                hybrid = build_hybrid_sidm_profile(
                    r_kpc=r_3d_kpc,
                    m200_msun=halo["m200_msun"],
                    c200=halo["c200"],
                    z=halo["z"],
                    sigma_over_m=sigma_over_m,
                    cosmo=DEFAULT_COSMOLOGY,
                    model_options=model_options,
                    dk14_params={"regime": tier2_configuration.get("regime", "cluster")},
                    stitch_params=tier2_configuration,
                )
                if sigma_over_m == sigma_grid[0]:
                    cdm_projected_tier2 = delta_sigma_of_R(
                        r_projected_kpc=r_lensing_kpc,
                        r_kpc=r_3d_kpc,
                        rho_msun_kpc3=hybrid["rho_cdm_hybrid_msun_kpc3"],
                        n_z=benchmark.projection_n_z,
                    )
                    cdm_projected_profiles_tier2.append(cdm_projected_tier2)
                    cdm_rho_rows_tier2.append(
                        interpolate_profile_to_common_grid(
                            r_input_kpc=r_3d_kpc,
                            profile_input=hybrid["rho_cdm_hybrid_msun_kpc3"],
                            r_common_kpc=r_common_3d_kpc,
                        )
                    )

                sidm_projected_tier2 = delta_sigma_of_R(
                    r_projected_kpc=r_lensing_kpc,
                    r_kpc=r_3d_kpc,
                    rho_msun_kpc3=hybrid["rho_sidm_hybrid_msun_kpc3"],
                    n_z=benchmark.projection_n_z,
                )
                sidm_projected_profiles_tier2[sigma_over_m].append(sidm_projected_tier2)
                sidm_rho_rows_tier2[sigma_over_m].append(
                    interpolate_profile_to_common_grid(
                        r_input_kpc=r_3d_kpc,
                        profile_input=hybrid["rho_sidm_hybrid_msun_kpc3"],
                        r_common_kpc=r_common_3d_kpc,
                    )
                )

                if halo_index == 0 and sigma_over_m == sigma_grid[-1]:
                    single_halo_records["tier2"] = {
                        "r_kpc": r_3d_kpc,
                        "rho_cdm_inner": cdm_profile["rho_msun_kpc3"],
                        "rho_sidm_inner": sidm_profile["rho_msun_kpc3"],
                        "rho_dk14_reference": hybrid["dk14_outer_reference"][
                            "rho_total_msun_kpc3"
                        ],
                        "rho_hybrid": hybrid["rho_sidm_hybrid_msun_kpc3"],
                        "r_projected_kpc": r_lensing_kpc,
                        "delta_sigma_tier1": sidm_projected["delta_sigma_msun_kpc2"],
                        "delta_sigma_tier2": sidm_projected_tier2["delta_sigma_msun_kpc2"],
                    }

    stacked_cdm = stack_delta_sigma_profiles(
        r_common_kpc=r_common_lensing_kpc,
        projected_profiles=cdm_projected_profiles,
        weights=weights,
    )
    stacked_cdm_rho = weighted_stack_profiles(np.vstack(cdm_rho_rows), weights)

    stacked_sidm_by_sigma: dict[float, np.ndarray] = {}
    stacked_sidm_rho_by_sigma: dict[float, np.ndarray] = {}
    for sigma_over_m in sigma_grid:
        stacked_sidm_by_sigma[sigma_over_m] = stack_delta_sigma_profiles(
            r_common_kpc=r_common_lensing_kpc,
            projected_profiles=sidm_projected_profiles[sigma_over_m],
            weights=weights,
        )["delta_sigma_msun_kpc2"]
        stacked_sidm_rho_by_sigma[sigma_over_m] = weighted_stack_profiles(
            np.vstack(sidm_rho_rows[sigma_over_m]),
            weights,
        )

    stacked_cdm_tier2: dict[str, np.ndarray] | None = None
    stacked_cdm_rho_tier2: np.ndarray | None = None
    stacked_sidm_by_sigma_tier2: dict[float, np.ndarray] = {}
    stacked_sidm_rho_by_sigma_tier2: dict[float, np.ndarray] = {}
    if tier2_enabled:
        stacked_cdm_tier2 = stack_delta_sigma_profiles(
            r_common_kpc=r_common_lensing_kpc,
            projected_profiles=cdm_projected_profiles_tier2,
            weights=weights,
        )
        stacked_cdm_rho_tier2 = weighted_stack_profiles(np.vstack(cdm_rho_rows_tier2), weights)
        for sigma_over_m in sigma_grid:
            stacked_sidm_by_sigma_tier2[sigma_over_m] = stack_delta_sigma_profiles(
                r_common_kpc=r_common_lensing_kpc,
                projected_profiles=sidm_projected_profiles_tier2[sigma_over_m],
                weights=weights,
            )["delta_sigma_msun_kpc2"]
            stacked_sidm_rho_by_sigma_tier2[sigma_over_m] = weighted_stack_profiles(
                np.vstack(sidm_rho_rows_tier2[sigma_over_m]),
                weights,
            )

    scenario_fractions = {"optimistic_5pct": 0.05, "conservative_10pct": 0.10}
    summary_rows: list[dict[str, float | int | str]] = []
    delta_chi2_by_sigma_tier1: dict[float, float] = {}
    delta_chi2_by_sigma_tier2: dict[float, float] = {}
    for sigma_over_m in sigma_grid:
        metrics_tier1 = evaluate_stacked_distinguishability(
            model=stacked_sidm_by_sigma[sigma_over_m],
            reference=stacked_cdm["delta_sigma_msun_kpc2"],
            fractional_error_scenarios=scenario_fractions,
        )
        required_for_3sigma_tier1 = required_uniform_fractional_precision(
            model=stacked_sidm_by_sigma[sigma_over_m],
            reference=stacked_cdm["delta_sigma_msun_kpc2"],
            target_sigma_significance=3.0,
        )
        max_ratio_shift_percent_tier1 = 100.0 * np.max(
            np.abs(stacked_sidm_by_sigma[sigma_over_m] / stacked_cdm["delta_sigma_msun_kpc2"] - 1.0)
        )

        summary_rows.append(
            {
                "label": mode_tag,
                "n_halos": int(sampling_configuration["n_halos"]),
                "seed": int(sampling_configuration["seed"]),
                "ensemble_mode": active_mode,
                "tier": "tier1",
                "sigma_over_m_cm2_g": sigma_over_m,
                "delta_chi2_optimistic_5pct": metrics_tier1["delta_chi2_optimistic_5pct"],
                "sigma_separation_optimistic_5pct": metrics_tier1["sigma_separation_optimistic_5pct"],
                "delta_chi2_conservative_10pct": metrics_tier1["delta_chi2_conservative_10pct"],
                "sigma_separation_conservative_10pct": metrics_tier1[
                    "sigma_separation_conservative_10pct"
                ],
                "required_uniform_fraction_for_3sigma": required_for_3sigma_tier1,
                "required_uniform_percent_for_3sigma": 100.0 * required_for_3sigma_tier1,
                "max_abs_delta_sigma_ratio_shift_percent": max_ratio_shift_percent_tier1,
                "ensemble_mean_mass_msun": halo_summary["mean_mass_msun"],
                "ensemble_median_mass_msun": float(
                    np.median([halo["m200_msun"] for halo in halos])
                ),
                "ensemble_std_log10_mass_dex": halo_summary["std_log10_mass_dex"],
            }
        )
        delta_chi2_by_sigma_tier1[sigma_over_m] = metrics_tier1["delta_chi2_optimistic_5pct"]

        if tier2_enabled and stacked_cdm_tier2 is not None:
            metrics_tier2 = evaluate_stacked_distinguishability(
                model=stacked_sidm_by_sigma_tier2[sigma_over_m],
                reference=stacked_cdm_tier2["delta_sigma_msun_kpc2"],
                fractional_error_scenarios=scenario_fractions,
            )
            required_for_3sigma_tier2 = required_uniform_fractional_precision(
                model=stacked_sidm_by_sigma_tier2[sigma_over_m],
                reference=stacked_cdm_tier2["delta_sigma_msun_kpc2"],
                target_sigma_significance=3.0,
            )
            max_ratio_shift_percent_tier2 = 100.0 * np.max(
                np.abs(
                    stacked_sidm_by_sigma_tier2[sigma_over_m]
                    / stacked_cdm_tier2["delta_sigma_msun_kpc2"]
                    - 1.0
                )
            )
            summary_rows.append(
                {
                    "label": mode_tag,
                    "n_halos": int(sampling_configuration["n_halos"]),
                    "seed": int(sampling_configuration["seed"]),
                    "ensemble_mode": active_mode,
                    "tier": "tier2",
                    "sigma_over_m_cm2_g": sigma_over_m,
                    "delta_chi2_optimistic_5pct": metrics_tier2["delta_chi2_optimistic_5pct"],
                    "sigma_separation_optimistic_5pct": metrics_tier2[
                        "sigma_separation_optimistic_5pct"
                    ],
                    "delta_chi2_conservative_10pct": metrics_tier2[
                        "delta_chi2_conservative_10pct"
                    ],
                    "sigma_separation_conservative_10pct": metrics_tier2[
                        "sigma_separation_conservative_10pct"
                    ],
                    "required_uniform_fraction_for_3sigma": required_for_3sigma_tier2,
                    "required_uniform_percent_for_3sigma": 100.0 * required_for_3sigma_tier2,
                    "max_abs_delta_sigma_ratio_shift_percent": max_ratio_shift_percent_tier2,
                    "ensemble_mean_mass_msun": halo_summary["mean_mass_msun"],
                    "ensemble_median_mass_msun": float(
                        np.median([halo["m200_msun"] for halo in halos])
                    ),
                    "ensemble_std_log10_mass_dex": halo_summary["std_log10_mass_dex"],
                }
            )
            delta_chi2_by_sigma_tier2[sigma_over_m] = metrics_tier2["delta_chi2_optimistic_5pct"]

    stacked_delta_sigma_table = pd.DataFrame(
        {
            "r_projected_kpc": r_common_lensing_kpc,
            "delta_sigma_cdm_tier1_msun_kpc2": stacked_cdm["delta_sigma_msun_kpc2"],
            **{
                f"delta_sigma_sidm_tier1_sigma_{sigma_value:.3f}".replace(".", "p"): values
                for sigma_value, values in stacked_sidm_by_sigma.items()
            },
        }
    )
    if tier2_enabled and stacked_cdm_tier2 is not None:
        stacked_delta_sigma_table["delta_sigma_cdm_tier2_msun_kpc2"] = stacked_cdm_tier2[
            "delta_sigma_msun_kpc2"
        ]
        for sigma_value, values in stacked_sidm_by_sigma_tier2.items():
            stacked_delta_sigma_table[
                f"delta_sigma_sidm_tier2_sigma_{sigma_value:.3f}".replace(".", "p")
            ] = values
    save_table(
        stacked_delta_sigma_table,
        output_paths["intermediate"] / f"{mode_tag}_ensemble_stacked_delta_sigma_profiles.csv",
    )

    stacked_rho_table = pd.DataFrame(
        {
            "r_kpc": r_common_3d_kpc,
            "rho_cdm_tier1_msun_kpc3": stacked_cdm_rho,
            **{
                f"rho_sidm_tier1_sigma_{sigma_value:.3f}".replace(".", "p"): values
                for sigma_value, values in stacked_sidm_rho_by_sigma.items()
            },
        }
    )
    if tier2_enabled and stacked_cdm_rho_tier2 is not None:
        stacked_rho_table["rho_cdm_tier2_msun_kpc3"] = stacked_cdm_rho_tier2
        for sigma_value, values in stacked_sidm_rho_by_sigma_tier2.items():
            stacked_rho_table[f"rho_sidm_tier2_sigma_{sigma_value:.3f}".replace(".", "p")] = values
    save_table(
        stacked_rho_table,
        output_paths["intermediate"] / f"{mode_tag}_ensemble_stacked_rho_profiles.csv",
    )

    save_table(
        pd.DataFrame(summary_rows),
        output_paths["tables"] / f"{mode_tag}_ensemble_delta_chi2_summary.csv",
    )

    expected_bins = (
        int(projection_configuration["n_r_bins"])
        if projection_configuration is not None
        else int(len(r_common_lensing_kpc))
    )
    validation_table = _validation_table(
        label=mode_tag,
        ensemble_mode=active_mode,
        halos=halos,
        r_common_lensing_kpc=r_common_lensing_kpc,
        expected_projection_bins=expected_bins,
        sampling_configuration=sampling_configuration,
    )
    save_table(
        validation_table,
        output_paths["tables"] / f"{mode_tag}_ensemble_validation_checks.csv",
    )

    plot_stacked_delta_sigma(
        r_projected_kpc=r_common_lensing_kpc,
        delta_sigma_cdm_msun_kpc2=stacked_cdm["delta_sigma_msun_kpc2"],
        sidm_delta_sigma_profiles=stacked_sidm_by_sigma,
        output_path=output_paths["figures"] / f"{mode_tag}_ensemble_stacked_delta_sigma.png",
    )
    plot_tier1_summary(
        r_3d_kpc=r_common_3d_kpc,
        rho_cdm_msun_kpc3=stacked_cdm_rho,
        rho_sidm_profiles=stacked_sidm_rho_by_sigma,
        r_projected_kpc=r_common_lensing_kpc,
        delta_sigma_cdm_msun_kpc2=stacked_cdm["delta_sigma_msun_kpc2"],
        delta_sigma_sidm_profiles=stacked_sidm_by_sigma,
        sigma_over_m_grid_cm2_g=sigma_grid,
        delta_chi2_by_sigma=delta_chi2_by_sigma_tier1,
        output_path=output_paths["figures"] / f"{mode_tag}_ensemble_tier1_summary.png",
    )

    figure_paths = [
        output_paths["figures"] / f"{mode_tag}_ensemble_stacked_delta_sigma.png",
        output_paths["figures"] / f"{mode_tag}_ensemble_tier1_summary.png",
    ]
    if tier2_enabled and stacked_cdm_tier2 is not None:
        plot_tier1_tier2_stacked_comparison(
            r_projected_kpc=r_common_lensing_kpc,
            cdm_tier1=stacked_cdm["delta_sigma_msun_kpc2"],
            sidm_tier1_by_sigma=stacked_sidm_by_sigma,
            cdm_tier2=stacked_cdm_tier2["delta_sigma_msun_kpc2"],
            sidm_tier2_by_sigma=stacked_sidm_by_sigma_tier2,
            output_path=output_paths["figures"] / f"{mode_tag}_ensemble_tier1_vs_tier2_stacked.png",
            title=f"{mode_tag} ensemble",
        )
        plot_delta_chi2_tier_comparison(
            sigma_over_m_grid_cm2_g=sigma_grid,
            delta_chi2_tier1=delta_chi2_by_sigma_tier1,
            delta_chi2_tier2=delta_chi2_by_sigma_tier2,
            output_path=output_paths["figures"] / f"{mode_tag}_ensemble_tier1_vs_tier2_delta_chi2.png",
            title=f"{mode_tag} ensemble",
        )
        figure_paths.extend(
            [
                output_paths["figures"] / f"{mode_tag}_ensemble_tier1_vs_tier2_stacked.png",
                output_paths["figures"] / f"{mode_tag}_ensemble_tier1_vs_tier2_delta_chi2.png",
            ]
        )
        if "tier2" in single_halo_records:
            single_halo_output = output_paths["figures"] / f"{mode_tag}_single_halo_tier2_diagnostic.png"
            plot_single_halo_tier2_diagnostics(
                r_kpc=single_halo_records["tier2"]["r_kpc"],
                rho_cdm_inner=single_halo_records["tier2"]["rho_cdm_inner"],
                rho_sidm_inner=single_halo_records["tier2"]["rho_sidm_inner"],
                rho_dk14_reference=single_halo_records["tier2"]["rho_dk14_reference"],
                rho_hybrid=single_halo_records["tier2"]["rho_hybrid"],
                r_projected_kpc=single_halo_records["tier2"]["r_projected_kpc"],
                delta_sigma_tier1=single_halo_records["tier2"]["delta_sigma_tier1"],
                delta_sigma_tier2=single_halo_records["tier2"]["delta_sigma_tier2"],
                output_path=single_halo_output,
                title=f"{mode_tag} representative halo",
            )
            figure_paths.append(single_halo_output)

    caption_entries: list[dict[str, str]] = []
    for figure_path in figure_paths:
        timestamp = pd.Timestamp(figure_path.stat().st_mtime, unit="s").tz_localize("UTC").tz_convert(
            "America/Phoenix"
        )
        caption_entries.append(
            {
                "filename": figure_path.name,
                "created_at": timestamp.isoformat(),
                "caption": (
                    "Tier forecast figure. Axes are in physical kpc and Msun-based density units; "
                    "Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts "
                    "attachment with smooth log-density stitching."
                ),
            }
        )
    append_figure_caption_entries(
        caption_path=output_paths["figures"] / "CAPTION.md",
        entries=caption_entries,
    )

    inventory_entries = [
        {
            "path": str((output_paths["intermediate"] / f"{mode_tag}_ensemble_halo_catalog.csv").relative_to(output_root)),
            "description": "Sampled halo catalog with M200, c200, z, and stack weights.",
        },
        {
            "path": str((output_paths["intermediate"] / f"{mode_tag}_ensemble_stacked_delta_sigma_profiles.csv").relative_to(output_root)),
            "description": "Stacked DeltaSigma profiles for CDM/SIDM and Tier-1/Tier-2 when enabled.",
        },
        {
            "path": str((output_paths["intermediate"] / f"{mode_tag}_ensemble_stacked_rho_profiles.csv").relative_to(output_root)),
            "description": "Stacked rho(r) profiles for CDM/SIDM and Tier-1/Tier-2 when enabled.",
        },
        {
            "path": str((output_paths["tables"] / f"{mode_tag}_ensemble_delta_chi2_summary.csv").relative_to(output_root)),
            "description": "Distinguishability summary with DeltaChi2 metrics by sigma/m and tier.",
        },
        {
            "path": str((output_paths["tables"] / f"{mode_tag}_ensemble_validation_checks.csv").relative_to(output_root)),
            "description": "Validation checks for median mass, bins, reproducibility, and single-halo limit.",
        },
    ]
    append_inventory_entries(
        inventory_path=output_root / "INVENTORY.md",
        entries=inventory_entries,
    )


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("outputs"),
        help="Output directory root (default: outputs).",
    )
    parser.add_argument(
        "--sidm-backend",
        choices=["parametric", "surrogate"],
        default="parametric",
        help="Use parametricSIDM backend or deterministic surrogate fallback.",
    )
    parser.add_argument(
        "--n-halos",
        type=int,
        default=None,
        help="Optional halo count override. Uses YAML/default value when omitted.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Optional random-seed override. Uses YAML/default value when omitted.",
    )
    parser.add_argument(
        "--ensemble-mode",
        choices=["HMF", "SHMR"],
        default="HMF",
        help="Choose HMF-based or SHMR-based halo-ensemble sampling.",
    )
    parser.add_argument(
        "--config-path",
        type=Path,
        default=None,
        help="Optional YAML config path. If set, mode and sampler settings come from YAML.",
    )
    return parser.parse_args()


def main() -> None:
    """Program entrypoint."""
    arguments = parse_arguments()
    run_ensemble_forecast(
        output_root=arguments.output_root,
        sidm_backend=arguments.sidm_backend,
        n_halos=arguments.n_halos,
        seed=arguments.seed,
        ensemble_mode=arguments.ensemble_mode,
        config_path=arguments.config_path,
    )


if __name__ == "__main__":
    main()
