"""
Generate key summary figures for the dwarf-regime forecast.

Produces:
    1. Parameter error comparison bar chart (the key figure)
    2. DeltaSigma(R) with error bars in a representative dwarf bin
    3. Improvement factor table (printed)

Usage:
    MPLBACKEND=Agg python scripts/generate_dwarf_figures.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from shmr_fisher.config_io import load_run_config
from shmr_fisher.plot_results import (
    _assign_survey_colors,
    _save_or_show,
    PARAM_LATEX,
)
from scripts.run_forecast import run_forecast_from_config


def plot_dwarf_error_comparison(
    forecast_results: dict,
    shmr_params,
    save_path: str | Path | None = None,
) -> None:
    """
    Bar chart of fractional parameter errors for the dwarf-regime forecast.

    Only shows the varied parameters (N_0, beta_0, scatter_sigma_high,
    scatter_sigma_rise). Uses log scale since errors span orders of magnitude.

    Parameters
    ----------
    forecast_results : dict
        Output of run_forecast(), keyed by survey name.
    shmr_params : SHMRParams
        Fiducial SHMR parameters.
    save_path : str, Path, or None
        If provided, save figure here. Otherwise display.
    """
    # Get varied params from the first result
    first_result = next(iter(forecast_results.values()))
    first_errors = first_result.get("shmr_errors", first_result.get("errors", {}))
    varied_params = list(first_errors.keys())

    survey_keys = list(forecast_results.keys())
    survey_color_map = _assign_survey_colors(survey_keys)
    n_surveys = len(survey_keys)

    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(varied_params))
    width = min(0.8 / max(n_surveys, 1), 0.25)

    for i, sname in enumerate(survey_keys):
        r = forecast_results[sname]
        color = survey_color_map[sname]
        errs = r.get("shmr_errors", r.get("errors", {}))

        frac_errors = []
        for p in varied_params:
            sigma = errs.get(p, 0)
            fid = abs(getattr(shmr_params, p))
            frac = sigma / fid * 100 if fid > 1e-10 else 0
            frac_errors.append(frac)

        bars = ax.bar(
            x + i * width, frac_errors, width,
            label=r["survey_name"], color=color, edgecolor="white", lw=0.5,
        )
        # Add value labels on bars
        for bar, val in zip(bars, frac_errors):
            if val > 0:
                ax.text(
                    bar.get_x() + bar.get_width() / 2, bar.get_height() * 1.1,
                    f"{val:.1f}%", ha="center", va="bottom", fontsize=7,
                    rotation=45,
                )

    ax.set_xlabel("SHMR Parameter", fontsize=13)
    ax.set_ylabel(r"$\sigma / |\theta_\mathrm{fid}|$ [%]", fontsize=13)
    ax.set_title(
        "Dwarf-regime SHMR constraints: DESI → DESI-II → Spec-S5",
        fontsize=13,
    )
    ax.set_xticks(x + width * (n_surveys - 1) / 2)
    ax.set_xticklabels(
        [PARAM_LATEX.get(p, p) for p in varied_params], fontsize=12,
    )
    ax.set_yscale("log")
    ax.set_ylim(bottom=0.01)
    ax.legend(fontsize=10, loc="upper left")
    ax.grid(True, alpha=0.3, axis="y")

    # Add horizontal reference lines
    ax.axhline(1, color="gray", ls="--", lw=0.8, alpha=0.5)
    ax.axhline(10, color="gray", ls="--", lw=0.8, alpha=0.5)
    ax.axhline(100, color="gray", ls="--", lw=0.8, alpha=0.5)

    plt.tight_layout()
    _save_or_show(fig, save_path)


def plot_dwarf_delta_sigma(
    forecast_results: dict,
    shmr_params,
    save_path: str | Path | None = None,
) -> None:
    """
    DeltaSigma(R) with error bars in a representative dwarf stellar mass bin.

    Uses the (8.0, 8.5) bin at z=0.05, which is common to all three surveys.

    Parameters
    ----------
    forecast_results : dict
        Output of run_forecast().
    shmr_params : SHMRParams
        Fiducial SHMR parameters.
    save_path : str, Path, or None
    """
    from shmr_fisher.halo_model import delta_sigma_bin
    from shmr_fisher.covariance import lensing_covariance
    from shmr_fisher.config import LensingConfig

    z = 0.05
    ms_lo, ms_hi = 8.0, 8.5

    R_Mpc = np.logspace(np.log10(0.05), np.log10(10.0), 30)
    ds = delta_sigma_bin(R_Mpc, ms_lo, ms_hi, shmr_params, z)
    ds_pc2 = ds / 1e12  # Msun/Mpc^2 -> Msun/pc^2

    R_edges = np.logspace(np.log10(0.05), np.log10(10.0), 11)
    R_centers = np.sqrt(R_edges[:-1] * R_edges[1:])
    ds_centers = delta_sigma_bin(R_centers, ms_lo, ms_hi, shmr_params, z)
    ds_centers_pc2 = ds_centers / 1e12

    lc = LensingConfig()
    survey_keys = list(forecast_results.keys())
    survey_color_map = _assign_survey_colors(survey_keys)
    n_surveys = len(survey_keys)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(R_Mpc, ds_pc2, "k-", lw=2.5, label="Model (NFW + 2-halo)", zorder=10)

    # Offset surveys in R for visibility
    offsets = np.linspace(-0.03, 0.03, n_surveys)

    for i, sname in enumerate(survey_keys):
        r = forecast_results[sname]
        color = survey_color_map[sname]
        meta = r.get("metadata", {})

        # Get N_lens for this survey (sum over all bins, use average)
        N_lens_dict = meta.get("N_lens_per_bin", {})
        if N_lens_dict:
            N_lens = np.mean(list(N_lens_dict.values()))
        else:
            N_lens = r.get("n_gal_total", 1e5) / 3

        area = meta.get("area_deg2", 10000)

        var_ds = lensing_covariance(R_edges, z, N_lens, area, lc)
        f_sys = meta.get("systematic_floor_fraction", 0.05)
        if f_sys > 0:
            var_ds = var_ds + (f_sys * ds_centers) ** 2

        sigma_pc2 = np.sqrt(var_ds) / 1e12
        R_plot = R_centers * 10 ** offsets[i]

        ax.errorbar(
            R_plot, ds_centers_pc2, yerr=sigma_pc2,
            fmt="o", color=color, ms=4, capsize=2, lw=1.5,
            label=f"{r['survey_name']}",
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$R$ [physical Mpc]", fontsize=13)
    ax.set_ylabel(r"$\Delta\Sigma$ [$M_\odot/\mathrm{pc}^2$]", fontsize=13)
    ax.set_title(
        rf"Dwarf lensing signal: $\log M_* \in [{ms_lo}, {ms_hi}]$ at $z={z}$",
        fontsize=13,
    )
    ax.set_xlim(0.03, 15)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which="both")

    plt.tight_layout()
    _save_or_show(fig, save_path)


def main():
    config_path = "configs/dwarf_regime.yaml"
    run_config = load_run_config(config_path)
    outdir = run_config.output_dir
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print(f"Dwarf-regime forecast: {list(run_config.surveys.keys())}")
    print(f"Output: {outdir}")
    print("=" * 60)

    # Run the forecast
    print("\nRunning forecast...")
    results = run_forecast_from_config(run_config)

    # Figure 1: Parameter error comparison
    print("\n--- Generating parameter error comparison ---")
    plot_dwarf_error_comparison(
        results, run_config.shmr_params,
        save_path=outdir / "dwarf_error_comparison.png",
    )

    # Figure 2: DeltaSigma with error bars
    print("\n--- Generating DeltaSigma with errors ---")
    plot_dwarf_delta_sigma(
        results, run_config.shmr_params,
        save_path=outdir / "dwarf_delta_sigma.png",
    )

    # Print improvement table
    print("\n" + "=" * 60)
    print("Improvement summary (fractional error %)")
    print("=" * 60)
    header = f"{'Parameter':<25}"
    for sname in results:
        header += f"{results[sname]['survey_name']:>15}"
    print(header)
    print("-" * len(header))

    first_errors = next(iter(results.values()))
    params = list(first_errors.get("shmr_errors", first_errors.get("errors", {})).keys())
    for p in params:
        row = f"{p:<25}"
        for sname in results:
            errs = results[sname].get("shmr_errors", results[sname].get("errors", {}))
            fid = abs(getattr(run_config.shmr_params, p))
            sigma = errs.get(p, 0)
            frac = sigma / fid * 100 if fid > 1e-10 else 0
            row += f"{frac:>14.1f}%"
        print(row)

    print("\n" + "=" * 60)
    print(f"All figures saved to {outdir}")
    print("=" * 60)


if __name__ == "__main__":
    main()
