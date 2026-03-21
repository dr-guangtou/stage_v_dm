"""
Generate key summary figures for the dwarf-regime forecast.

Produces:
    1. Parameter error comparison bar chart (the key figure)
    2. DeltaSigma(R) with error bars in a representative dwarf bin
    3. Fisher uncertainty ellipses for key parameter pairs
    4. Relative improvement comparison (DESI as baseline)
    5. Improvement factor table (printed)

Usage:
    MPLBACKEND=Agg python scripts/generate_dwarf_figures.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
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


def _shmr_covariance(result: dict, param_names_filter: list[str]) -> np.ndarray:
    """
    Extract the SHMR parameter sub-block of the full covariance matrix.

    Inverts the full Fisher matrix (including nuisance params) and returns
    the covariance sub-block for the requested SHMR parameters.

    Parameters
    ----------
    result : dict
        Single survey result dict from run_forecast().
    param_names_filter : list of str
        SHMR parameter names to extract (must be a subset of param_names).

    Returns
    -------
    cov_sub : array, shape (len(param_names_filter), len(param_names_filter))
    """
    fisher = result["fisher"]
    all_names = result["param_names"]
    cov = np.linalg.inv(fisher)
    idx = [all_names.index(p) for p in param_names_filter]
    return cov[np.ix_(idx, idx)]


def plot_dwarf_fisher_ellipses(
    forecast_results: dict,
    shmr_params,
    save_path: str | Path | None = None,
) -> None:
    """
    2D Fisher forecast ellipses (1σ and 2σ) for key dwarf parameter pairs.

    Shows three panels: (N_0, β_0), (σ_high, σ_rise), (N_0, σ_high).
    Each panel overlays ellipses from all surveys.

    Parameters
    ----------
    forecast_results : dict
        Output of run_forecast(), keyed by survey name.
    shmr_params : SHMRParams
        Fiducial SHMR parameters.
    save_path : str, Path, or None
    """
    # Parameter pairs to display
    pairs = [
        ("N_0", "beta_0"),
        ("scatter_sigma_high", "scatter_sigma_rise"),
        ("N_0", "scatter_sigma_high"),
    ]

    survey_keys = list(forecast_results.keys())
    survey_color_map = _assign_survey_colors(survey_keys)

    # Get the SHMR param names from the first result
    first_result = next(iter(forecast_results.values()))
    shmr_names = first_result["shmr_param_names"]

    # Filter pairs to only those whose params exist
    pairs = [(a, b) for a, b in pairs if a in shmr_names and b in shmr_names]
    n_panels = len(pairs)
    if n_panels == 0:
        return

    fig, axes = plt.subplots(1, n_panels, figsize=(6 * n_panels, 5.5))
    if n_panels == 1:
        axes = [axes]

    for ax, (px, py) in zip(axes, pairs):
        fid_x = getattr(shmr_params, px)
        fid_y = getattr(shmr_params, py)

        for sname in survey_keys:
            r = forecast_results[sname]
            color = survey_color_map[sname]

            try:
                cov_2x2 = _shmr_covariance(r, [px, py])
            except (ValueError, np.linalg.LinAlgError):
                continue

            # Eigendecomposition for ellipse orientation and size
            eigvals, eigvecs = np.linalg.eigh(cov_2x2)
            # Ensure positive eigenvalues (can be tiny negative from numerics)
            eigvals = np.maximum(eigvals, 0)
            order = eigvals.argsort()[::-1]
            eigvals = eigvals[order]
            eigvecs = eigvecs[:, order]

            # Ellipse angle (degrees) — angle of the major axis
            angle = np.degrees(np.arctan2(eigvecs[1, 0], eigvecs[0, 0]))

            # 1σ (Δχ²=2.30 for 2 params) and 2σ (Δχ²=6.18)
            for n_sigma, chi2_scale, alpha in [(1, 2.30, 0.25), (2, 6.18, 0.10)]:
                w = 2 * np.sqrt(chi2_scale * eigvals[0])
                h = 2 * np.sqrt(chi2_scale * eigvals[1])
                ell = Ellipse(
                    xy=(fid_x, fid_y), width=w, height=h, angle=angle,
                    facecolor=color, alpha=alpha,
                    edgecolor=color, lw=1.5,
                    label=r["survey_name"] if n_sigma == 1 else None,
                )
                ax.add_patch(ell)

        # Fiducial point
        ax.plot(fid_x, fid_y, "k+", ms=12, mew=2, zorder=10)

        # Axis labels
        ax.set_xlabel(PARAM_LATEX.get(px, px), fontsize=14)
        ax.set_ylabel(PARAM_LATEX.get(py, py), fontsize=14)
        ax.grid(True, alpha=0.3)

        # Auto-scale axes: use the largest survey's 2σ extent
        max_dx, max_dy = 0, 0
        for sname in survey_keys:
            r = forecast_results[sname]
            try:
                cov_2x2 = _shmr_covariance(r, [px, py])
                sx = np.sqrt(cov_2x2[0, 0])
                sy = np.sqrt(cov_2x2[1, 1])
                max_dx = max(max_dx, sx)
                max_dy = max(max_dy, sy)
            except (ValueError, np.linalg.LinAlgError):
                continue
        # 3σ margin
        margin = 3.5
        ax.set_xlim(fid_x - margin * max_dx, fid_x + margin * max_dx)
        ax.set_ylim(fid_y - margin * max_dy, fid_y + margin * max_dy)

    # Single legend from the first panel
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="upper center",
        ncol=min(len(survey_keys), 4), fontsize=10,
        bbox_to_anchor=(0.5, 1.02),
    )

    fig.suptitle(
        "Fisher forecast: 1σ & 2σ confidence ellipses",
        fontsize=14, y=1.06,
    )
    plt.tight_layout()
    _save_or_show(fig, save_path)


def plot_dwarf_relative_improvement(
    forecast_results: dict,
    shmr_params,
    baseline_key: str | None = None,
    save_path: str | Path | None = None,
) -> None:
    """
    Relative improvement in parameter constraints compared to a baseline survey.

    Shows σ(baseline) / σ(survey) as a grouped bar chart. Values > 1 mean the
    survey is better than the baseline.

    Parameters
    ----------
    forecast_results : dict
        Output of run_forecast(), keyed by survey name.
    shmr_params : SHMRParams
        Fiducial SHMR parameters.
    baseline_key : str or None
        Key of the baseline survey. If None, uses the first survey.
    save_path : str, Path, or None
    """
    survey_keys = list(forecast_results.keys())
    if baseline_key is None:
        baseline_key = survey_keys[0]
    survey_color_map = _assign_survey_colors(survey_keys)

    baseline = forecast_results[baseline_key]
    base_errs = baseline.get("shmr_errors", baseline.get("errors", {}))
    varied_params = list(base_errs.keys())

    # Exclude the baseline from the comparison bars
    compare_keys = [k for k in survey_keys if k != baseline_key]
    n_compare = len(compare_keys)

    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(varied_params))
    width = min(0.8 / max(n_compare, 1), 0.25)

    for i, sname in enumerate(compare_keys):
        r = forecast_results[sname]
        color = survey_color_map[sname]
        errs = r.get("shmr_errors", r.get("errors", {}))

        ratios = []
        for p in varied_params:
            sigma_base = base_errs.get(p, np.inf)
            sigma_survey = errs.get(p, np.inf)
            if sigma_survey > 0 and sigma_base > 0:
                ratios.append(sigma_base / sigma_survey)
            else:
                ratios.append(0)

        bars = ax.bar(
            x + i * width, ratios, width,
            label=r["survey_name"], color=color,
            edgecolor="white", lw=0.5,
        )
        # Value labels
        for bar, val in zip(bars, ratios):
            if val > 0:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.05,
                    f"{val:.1f}×",
                    ha="center", va="bottom", fontsize=8, fontweight="bold",
                )

    # Reference line at 1
    ax.axhline(1.0, color="gray", ls="--", lw=1.5, alpha=0.7, zorder=0)

    baseline_name = baseline["survey_name"]
    ax.set_xlabel("SHMR Parameter", fontsize=13)
    ax.set_ylabel(
        rf"$\sigma(\mathrm{{{baseline_name}}}) \,/\, \sigma(\mathrm{{survey}})$",
        fontsize=13,
    )
    ax.set_title(
        f"Relative improvement over {baseline_name}",
        fontsize=14,
    )
    ax.set_xticks(x + width * (n_compare - 1) / 2)
    ax.set_xticklabels(
        [PARAM_LATEX.get(p, p) for p in varied_params], fontsize=12,
    )
    ax.legend(fontsize=10, loc="upper left")
    ax.grid(True, alpha=0.3, axis="y")

    # Log scale if range spans > 10x
    all_max = max(
        base_errs.get(p, 0) / max(
            forecast_results[k].get("shmr_errors", forecast_results[k].get("errors", {})).get(p, np.inf),
            1e-30,
        )
        for k in compare_keys
        for p in varied_params
    )
    if all_max > 10:
        ax.set_yscale("log")
        ax.set_ylim(bottom=0.5)

    plt.tight_layout()
    _save_or_show(fig, save_path)


def run_and_plot(config_path: str, prefix: str = "dwarf") -> dict:
    """
    Run a forecast from a YAML config and generate all four figures.

    Parameters
    ----------
    config_path : str
        Path to the YAML config file.
    prefix : str
        Filename prefix for the output figures.

    Returns
    -------
    results : dict
        Forecast results keyed by survey name.
    """
    run_config = load_run_config(config_path)
    outdir = run_config.output_dir
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print(f"Forecast: {list(run_config.surveys.keys())}")
    print(f"Config: {config_path}")
    print(f"Output: {outdir}")
    print("=" * 60)

    # Run the forecast
    print("\nRunning forecast...")
    results = run_forecast_from_config(run_config)

    # Figure 1: Parameter error comparison
    print("\n--- Generating parameter error comparison ---")
    plot_dwarf_error_comparison(
        results, run_config.shmr_params,
        save_path=outdir / f"{prefix}_error_comparison.png",
    )

    # Figure 2: DeltaSigma with error bars
    print("\n--- Generating DeltaSigma with errors ---")
    plot_dwarf_delta_sigma(
        results, run_config.shmr_params,
        save_path=outdir / f"{prefix}_delta_sigma.png",
    )

    # Figure 3: Fisher uncertainty ellipses
    print("\n--- Generating Fisher ellipses ---")
    plot_dwarf_fisher_ellipses(
        results, run_config.shmr_params,
        save_path=outdir / f"{prefix}_fisher_ellipses.png",
    )

    # Figure 4: Relative improvement (DESI as baseline)
    print("\n--- Generating relative improvement ---")
    baseline_key = list(results.keys())[0]  # First survey is baseline
    plot_dwarf_relative_improvement(
        results, run_config.shmr_params,
        baseline_key=baseline_key,
        save_path=outdir / f"{prefix}_relative_improvement.png",
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

    return results


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate dwarf-regime forecast figures.",
    )
    parser.add_argument(
        "--config", default="configs/dwarf_regime.yaml",
        help="Path to YAML config (default: configs/dwarf_regime.yaml)",
    )
    parser.add_argument(
        "--prefix", default="dwarf",
        help="Filename prefix for output figures (default: dwarf)",
    )
    args = parser.parse_args()

    run_and_plot(args.config, prefix=args.prefix)


if __name__ == "__main__":
    main()
