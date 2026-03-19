"""Build summary-only Tier-3 overlay figures for split ensembles across redshift."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from sidm_stagev_forecast.config import DEFAULT_COSMOLOGY
from sidm_stagev_forecast.io import append_figure_caption_entries
from sidm_stagev_forecast.plotting import plot_tier3_redshift_overlay_summary


def _load_tier3_dataset(output_root: Path, label: str, redshift_label: str) -> dict[str, object]:
    stacked_rho = pd.read_csv(output_root / "intermediate" / f"{label}_ensemble_stacked_rho_profiles.csv")
    stacked_ds = pd.read_csv(output_root / "intermediate" / f"{label}_ensemble_stacked_delta_sigma_profiles.csv")
    summary = pd.read_csv(output_root / "tables" / f"{label}_ensemble_delta_chi2_summary.csv")
    halo_catalog = pd.read_csv(output_root / "intermediate" / f"{label}_ensemble_halo_catalog.csv")

    tier3_summary = summary[summary["tier"] == "tier3"].copy()
    sigma_grid = tuple(sorted(float(value) for value in tier3_summary["sigma_over_m_cm2_g"].unique()))

    rho_sidm_profiles: dict[float, pd.Series] = {}
    ds_sidm_profiles: dict[float, pd.Series] = {}
    for sigma in sigma_grid:
        sigma_tag = f"{sigma:.3f}".replace(".", "p")
        rho_sidm_profiles[sigma] = stacked_rho[f"rho_sidm_tier3_sigma_{sigma_tag}"]
        ds_sidm_profiles[sigma] = stacked_ds[f"delta_sigma_sidm_tier3_sigma_{sigma_tag}"]

    delta_chi2_by_sigma = {
        float(row["sigma_over_m_cm2_g"]): float(row["delta_chi2_optimistic_5pct"])
        for _, row in tier3_summary.iterrows()
    }

    return {
        "label": redshift_label,
        "sigma_grid": sigma_grid,
        "r_3d_kpc": stacked_rho["r_kpc"],
        "rho_cdm_msun_kpc3": stacked_rho["rho_cdm_tier3_msun_kpc3"],
        "rho_sidm_profiles": rho_sidm_profiles,
        "r_projected_kpc": stacked_ds["r_projected_kpc"],
        "delta_sigma_cdm_msun_kpc2": stacked_ds["delta_sigma_cdm_tier3_msun_kpc2"],
        "delta_sigma_sidm_profiles": ds_sidm_profiles,
        "delta_chi2_by_sigma": delta_chi2_by_sigma,
        "m200_msun_samples": halo_catalog["m200_msun"].to_numpy(dtype=float),
        "c200_samples": halo_catalog["c200"].to_numpy(dtype=float),
        "weights": halo_catalog["weight"].to_numpy(dtype=float),
    }


def _build_annotation_text(config_path: Path) -> str:
    import yaml

    with config_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    concentration_scatter = float(config["concentration"]["scatter_logc"])
    mode = str(config["mode"]).upper()
    if mode == "HMF":
        mass_line = (
            "M_halo [Msun]: "
            f"[{float(config['mass_function']['M_min']):.2e}, {float(config['mass_function']['M_max']):.2e}]"
        )
        stellar_line = "log10(M*/Msun): N/A (HMF mode)"
    else:
        stellar_line = (
            "log10(M*/Msun): "
            f"[{float(config['stellar_mass_distribution']['logMstar_min']):.1f}, "
            f"{float(config['stellar_mass_distribution']['logMstar_max']):.1f}]"
        )
        # Keep halo range broad but explicit from SHMR normalization block.
        mass_line = "M_halo [Msun]: SHMR-mapped sample"

    return f"{stellar_line}\n{mass_line}\nc scatter: {concentration_scatter:.2f} dex"


def _weighted_quantile(values: pd.Series, quantile: float, weights: pd.Series) -> float:
    ordered = np.argsort(values.to_numpy(dtype=float))
    sorted_values = values.to_numpy(dtype=float)[ordered]
    sorted_weights = weights.to_numpy(dtype=float)[ordered]
    cumulative = np.cumsum(sorted_weights)
    normalized = cumulative / cumulative[-1]
    return float(np.interp(quantile, normalized, sorted_values))


def _build_distribution_stats_text(datasets: list[dict[str, object]]) -> str:
    lines: list[str] = []
    for dataset in datasets:
        label = str(dataset["label"])
        masses = pd.Series(dataset["m200_msun_samples"], dtype=float)
        concentrations = pd.Series(dataset["c200_samples"], dtype=float)
        weights = pd.Series(dataset["weights"], dtype=float)
        log_mass = np.log10(masses)

        mass_mean = float(np.average(log_mass, weights=weights))
        mass_p16 = _weighted_quantile(log_mass, 0.16, weights)
        mass_p84 = _weighted_quantile(log_mass, 0.84, weights)

        concentration_mean = float(np.average(concentrations, weights=weights))
        concentration_p16 = _weighted_quantile(concentrations, 0.16, weights)
        concentration_p84 = _weighted_quantile(concentrations, 0.84, weights)

        lines.append(
            (
                f"{label}: "
                f"log10M mean={mass_mean:.2f}, 16-84=[{mass_p16:.2f}, {mass_p84:.2f}]; "
                f"c mean={concentration_mean:.2f}, 16-84=[{concentration_p16:.2f}, {concentration_p84:.2f}]"
            )
        )
    return "\n".join(lines)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    output_root = root / "outputs"

    ensemble_groups = {
        "cluster1": [
            ("cluster1_z0p1", "z=0.1", root / "docs" / "cluster1_z0p1_ensemble_config.yaml"),
            ("cluster1_z0p5", "z=0.5", root / "docs" / "cluster1_z0p5_ensemble_config.yaml"),
        ],
        "cluster2": [
            ("cluster2_z0p1", "z=0.1", root / "docs" / "cluster2_z0p1_ensemble_config.yaml"),
            ("cluster2_z0p5", "z=0.5", root / "docs" / "cluster2_z0p5_ensemble_config.yaml"),
        ],
        "dwarf1": [
            ("dwarf1_z0p05", "z=0.05", root / "docs" / "dwarf1_z0p05_ensemble_config.yaml"),
            ("dwarf1_z0p2", "z=0.2", root / "docs" / "dwarf1_z0p2_ensemble_config.yaml"),
        ],
        "dwarf2": [
            ("dwarf2_z0p05", "z=0.05", root / "docs" / "dwarf2_z0p05_ensemble_config.yaml"),
            ("dwarf2_z0p2", "z=0.2", root / "docs" / "dwarf2_z0p2_ensemble_config.yaml"),
        ],
    }

    caption_entries: list[dict[str, str]] = []

    for ensemble_name, members in ensemble_groups.items():
        datasets = []
        sigma_grid_ref: tuple[float, ...] | None = None
        base_annotation_text = _build_annotation_text(members[0][2])

        for label, redshift_label, _ in members:
            dataset = _load_tier3_dataset(output_root=output_root, label=label, redshift_label=redshift_label)
            sigma_grid = dataset["sigma_grid"]
            if sigma_grid_ref is None:
                sigma_grid_ref = sigma_grid
            elif sigma_grid != sigma_grid_ref:
                raise RuntimeError(f"Sigma grids do not match for {ensemble_name}: {sigma_grid_ref} vs {sigma_grid}")
            datasets.append(dataset)

        if sigma_grid_ref is None:
            raise RuntimeError(f"No datasets loaded for {ensemble_name}.")
        distribution_stats_text = _build_distribution_stats_text(datasets)
        annotation_text = f"{base_annotation_text}\n{distribution_stats_text}"

        output_path = output_root / "figures" / f"{ensemble_name}_ensemble_tier3_summary.png"
        plot_tier3_redshift_overlay_summary(
            redshift_datasets=datasets,
            sigma_over_m_grid_cm2_g=sigma_grid_ref,
            output_path=output_path,
            annotation_text=annotation_text,
            h=DEFAULT_COSMOLOGY.h,
        )

        created_at = pd.Timestamp(output_path.stat().st_mtime, unit="s").tz_localize("UTC").tz_convert(
            "America/Phoenix"
        ).isoformat()
        caption_entries.append(
            {
                "filename": output_path.name,
                "created_at": created_at,
                "caption": (
                    f"Tier-3 summary for {ensemble_name} with redshift overlay. Two line styles correspond to the two "
                    "fixed-redshift runs in that ensemble split; color encodes sigma/m. Top-right panel uses "
                    "DeltaSigma units of h Msun / pc^2. Annotation box lists stellar-mass range (or N/A for HMF), "
                    "halo-mass range definition, and concentration scatter used in the ensemble configuration."
                ),
            }
        )

    append_figure_caption_entries(output_root / "figures" / "CAPTION.md", caption_entries)


if __name__ == "__main__":
    main()
