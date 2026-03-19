"""Run optional NFW projection reference cross-check and save artifacts."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sidm_stagev_forecast.io import ensure_output_directories, save_table
from sidm_stagev_forecast.reference_validation import run_nfw_reference_crosscheck


def _plot_crosscheck(table: pd.DataFrame, output_path: Path) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

    axes[0].loglog(
        table["r_projected_kpc"],
        table["delta_sigma_local_msun_kpc2"],
        color="black",
        lw=2.0,
        label="Local numeric",
    )
    axes[0].loglog(
        table["r_projected_kpc"],
        table["delta_sigma_analytic_msun_kpc2"],
        color="tab:blue",
        lw=1.8,
        label="Analytic NFW",
    )
    if "delta_sigma_colossus_msun_kpc2" in table.columns:
        axes[0].loglog(
            table["r_projected_kpc"],
            table["delta_sigma_colossus_msun_kpc2"],
            color="tab:green",
            lw=1.8,
            label="Colossus",
        )
    axes[0].set_xlabel("R [kpc]")
    axes[0].set_ylabel(r"$\Delta\Sigma$ [Msun / kpc$^2$]")
    axes[0].set_title("NFW Reference Cross-Check")
    axes[0].legend(fontsize=8)

    axes[1].axhline(0.0, color="black", lw=1.0, ls="--")
    axes[1].semilogx(
        table["r_projected_kpc"],
        100.0 * table["frac_diff_local_vs_analytic"],
        color="tab:blue",
        lw=1.8,
        label="Local vs Analytic",
    )
    if "frac_diff_local_vs_colossus" in table.columns:
        axes[1].semilogx(
            table["r_projected_kpc"],
            100.0 * table["frac_diff_local_vs_colossus"],
            color="tab:green",
            lw=1.8,
            label="Local vs Colossus",
        )
    axes[1].set_xlabel("R [kpc]")
    axes[1].set_ylabel("Fractional difference [%]")
    axes[1].set_title("Cross-Check Residuals")
    axes[1].legend(fontsize=8)

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
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    output_paths = ensure_output_directories(arguments.output_root)

    table, backend_status = run_nfw_reference_crosscheck()

    table_path = output_paths["tables"] / "nfw_reference_crosscheck.csv"
    save_table(table, table_path)

    figure_path = output_paths["figures"] / "nfw_reference_crosscheck.png"
    _plot_crosscheck(table, figure_path)

    status_table = pd.DataFrame([status.__dict__ for status in backend_status])
    status_path = output_paths["tables"] / "nfw_reference_backend_status.csv"
    save_table(status_table, status_path)

    max_local_vs_analytic_percent = float(
        100.0 * np.max(np.abs(table["frac_diff_local_vs_analytic"]))
    )
    print(f"Saved: {table_path}")
    print(f"Saved: {figure_path}")
    print(f"Saved: {status_path}")
    print(f"Max |local-analytic| fractional difference: {max_local_vs_analytic_percent:.3f}%")

    if "frac_diff_local_vs_colossus" in table.columns:
        max_local_vs_colossus_percent = float(
            100.0 * np.max(np.abs(table["frac_diff_local_vs_colossus"]))
        )
        print(f"Max |local-colossus| fractional difference: {max_local_vs_colossus_percent:.3f}%")


if __name__ == "__main__":
    main()
