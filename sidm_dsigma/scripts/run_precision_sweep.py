"""Run precision sweep for target significance over selected radial windows."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from sidm_stagev_forecast.io import ensure_output_directories, save_table
from sidm_stagev_forecast.precision_sweep import run_precision_sweep


def _plot_precision_sweep_for_benchmark(
    table: pd.DataFrame,
    benchmark_label: str,
    output_path: Path,
) -> None:
    subset = table[table["benchmark"] == benchmark_label]
    windows = ["full", "inner", "outer"]
    targets = [2.0, 3.0, 5.0]

    figure, axes = plt.subplots(1, 3, figsize=(14, 4.5), constrained_layout=True, sharey=True)

    for axis, window in zip(axes, windows, strict=True):
        data_window = subset[subset["window"] == window]
        for target in targets:
            data_target = data_window[data_window["target_sigma"] == target].sort_values(
                "sigma_over_m_cm2_g"
            )
            axis.plot(
                data_target["sigma_over_m_cm2_g"],
                data_target["required_uniform_percent"],
                marker="o",
                lw=1.7,
                label=f"{int(target)}$\\sigma$",
            )

        axis.set_xlabel(r"$\sigma/m$ [cm$^2$/g]")
        axis.set_title(f"{benchmark_label.capitalize()} - {window}")
        axis.grid(alpha=0.25)

    axes[0].set_ylabel("Required uniform precision [%]")
    axes[-1].legend(fontsize=8)

    figure.savefig(output_path, dpi=220)
    plt.close(figure)


def parse_arguments() -> argparse.Namespace:
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
    arguments = parse_arguments()
    output_paths = ensure_output_directories(arguments.output_root)

    table = run_precision_sweep(sidm_backend=arguments.sidm_backend)
    table_path = output_paths["tables"] / "precision_sweep.csv"
    save_table(table, table_path)

    for benchmark_label in sorted(table["benchmark"].unique()):
        figure_path = output_paths["figures"] / f"{benchmark_label}_precision_sweep.png"
        _plot_precision_sweep_for_benchmark(table, benchmark_label, figure_path)
        print(f"Saved: {figure_path}")

    print(f"Saved: {table_path}")


if __name__ == "__main__":
    main()
