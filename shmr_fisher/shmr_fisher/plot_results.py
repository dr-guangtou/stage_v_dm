"""
Visualization functions for the SHMR forecast pipeline.

All plot functions accept a save_path argument: if provided, the figure
is saved to that path; otherwise it is displayed interactively.
Default format: PNG at 300 dpi.

Functions
---------
plot_shmr_validation : M*/Mh vs Mh at multiple redshifts
plot_delta_sigma_with_errors : DS(R) with error bars per survey tier
plot_fisher_comparison : bar chart of sigma(theta) across surveys
plot_fisher_ellipses : 2D Fisher ellipses for parameter pairs
plot_two_regime_summary : z=0 shape (left) vs z-evolution (right)
plot_improvement_factor : Stage-V / Stage-IV error ratio
plot_scaling : sigma(param) vs a swept survey parameter
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .config import LensingConfig, SHMRParams, SurveyConfig, ForecastConfig
from .shmr_model import mean_log_Mstar


# Map parameter names to LaTeX labels
PARAM_LATEX = {
    "log_M1_0": r"$\log M_1$",
    "N_0": r"$N_0$",
    "beta_0": r"$\beta_0$",
    "gamma_0": r"$\gamma_0$",
    "sigma_logMs": r"$\sigma_{\log M_*}$",
    "nu_M1": r"$\nu_{M_1}$",
    "nu_N": r"$\nu_N$",
    "nu_beta": r"$\nu_\beta$",
    "nu_gamma": r"$\nu_\gamma$",
}


def _save_or_show(fig: plt.Figure, save_path: str | Path | None) -> None:
    """Save figure to path or show interactively."""
    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Saved: {save_path}")
        plt.close(fig)
    else:
        plt.show()


# -------------------------------------------------------------------
# Phase 1: SHMR validation
# -------------------------------------------------------------------

def plot_shmr_validation(
    params: SHMRParams | None = None,
    save_path: str | Path | None = None,
) -> None:
    """
    Plot M*/Mh vs Mh at z=0, 0.5, 1.0, 2.0.

    Parameters
    ----------
    params : SHMRParams or None
        If None, uses fiducial defaults.
    save_path : str, Path, or None
        If provided, save figure to this path.
    """
    if params is None:
        params = SHMRParams()

    log_Mh = np.linspace(10.0, 15.0, 500)
    redshifts = [0.0, 0.5, 1.0, 2.0]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]

    fig, ax = plt.subplots(figsize=(8, 6))

    for z, color in zip(redshifts, colors):
        log_Mstar = mean_log_Mstar(log_Mh, params, z)
        ratio = 10.0 ** (log_Mstar - log_Mh)
        ax.plot(log_Mh, ratio, color=color, lw=2, label=f"z = {z:.1f}")

    ax.set_xlabel(r"$\log_{10}(M_h / M_\odot)$", fontsize=14)
    ax.set_ylabel(r"$M_* / M_h$", fontsize=14)
    ax.set_title("SHMR: Moster+2013 parameterization", fontsize=14)
    ax.set_xlim(10.0, 15.0)
    ax.set_ylim(0, 0.07)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=12)

    log_Mstar_z0 = mean_log_Mstar(log_Mh, params, 0.0)
    ratio_z0 = 10.0 ** (log_Mstar_z0 - log_Mh)
    peak_idx = np.argmax(ratio_z0)
    ax.annotate(
        f"Peak: {ratio_z0[peak_idx]:.3f}\nat log Mh = {log_Mh[peak_idx]:.1f}",
        xy=(log_Mh[peak_idx], ratio_z0[peak_idx]),
        xytext=(13.0, 0.05),
        fontsize=10,
        arrowprops=dict(arrowstyle="->", color="gray"),
        ha="center",
    )

    plt.tight_layout()
    _save_or_show(fig, save_path)


# -------------------------------------------------------------------
# Phase 2: DeltaSigma with error bars
# -------------------------------------------------------------------

def plot_delta_sigma_with_errors(
    forecast_results: dict,
    shmr_params: SHMRParams | None = None,
    z: float = 0.3,
    mstar_bin: tuple[float, float] = (10.5, 11.0),
    save_path: str | Path | None = None,
) -> None:
    """
    Plot DeltaSigma(R) model prediction with error bars for different survey tiers.

    Parameters
    ----------
    forecast_results : dict
        Output of run_forecast(), keyed by survey name.
    shmr_params : SHMRParams or None
    z : float
        Redshift for the signal.
    mstar_bin : tuple of (lo, hi)
        Stellar mass bin [log10(Msun)].
    save_path : str, Path, or None
    """
    from .halo_model import delta_sigma_bin
    from .covariance import lensing_covariance

    if shmr_params is None:
        shmr_params = SHMRParams()

    R_Mpc = np.logspace(np.log10(0.1), np.log10(30.0), 25)
    ms_lo, ms_hi = mstar_bin
    ds = delta_sigma_bin(R_Mpc, ms_lo, ms_hi, shmr_params, z)
    ds_pc2 = ds / 1e12

    fig, ax = plt.subplots(figsize=(9, 6))

    # Model prediction
    ax.plot(R_Mpc, ds_pc2, "k-", lw=2.5, label="Model", zorder=10)

    # Error bars per survey
    survey_colors = {
        "stage3_shallow_wide": ("#1f77b4", "Stage-III"),
        "stage4_low_z": ("#ff7f0e", "Stage-IV"),
        "stage5_wide": ("#2ca02c", "Stage-V"),
    }

    R_edges = np.logspace(np.log10(0.1), np.log10(30.0), 11)
    R_centers = np.sqrt(R_edges[:-1] * R_edges[1:])
    ds_at_centers = delta_sigma_bin(R_centers, ms_lo, ms_hi, shmr_params, z)
    ds_at_centers_pc2 = ds_at_centers / 1e12

    lc = LensingConfig()
    offsets = {"stage3_shallow_wide": -0.03, "stage4_low_z": 0.0, "stage5_wide": 0.03}

    for sname, (color, label) in survey_colors.items():
        if sname not in forecast_results:
            continue
        meta = forecast_results[sname].get("metadata", {})
        # Find N_lens for a bin near our target
        N_lens_total = 0
        for key, obs in meta.get("fiducial_observables", {}).items():
            N_lens_total += obs.get("N_lens", 0)
        # Use a representative N_lens (total / number of bins)
        n_bins = max(1, len(meta.get("fiducial_observables", {})))
        N_lens = N_lens_total / n_bins

        var_ds = lensing_covariance(R_edges, z, N_lens, 10000, lc)
        # Add systematic floor if present
        f_sys = meta.get("systematic_floor_fraction",
                         forecast_results[sname].get("metadata", {}).get(
                             "systematic_floor_fraction", 0))
        if f_sys and f_sys > 0:
            var_ds = var_ds + (f_sys * ds_at_centers) ** 2

        sigma_pc2 = np.sqrt(var_ds) / 1e12

        # Offset slightly in R for visibility
        R_plot = R_centers * 10 ** offsets.get(sname, 0)
        ax.errorbar(R_plot, ds_at_centers_pc2, yerr=sigma_pc2,
                     fmt="o", color=color, ms=4, capsize=2, lw=1.5,
                     label=f"{label} (N$_{{lens}}$={N_lens:.0e})")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$R$ [physical Mpc]", fontsize=13)
    ax.set_ylabel(r"$\Delta\Sigma$ [$M_\odot/\mathrm{pc}^2$]", fontsize=13)
    ax.set_title(
        rf"$\log M_* \in [{ms_lo}, {ms_hi}]$ at $z={z}$", fontsize=13)
    ax.set_xlim(0.08, 40)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which="both")

    plt.tight_layout()
    _save_or_show(fig, save_path)


# -------------------------------------------------------------------
# Phase 4: Science figures
# -------------------------------------------------------------------

def plot_two_regime_summary(
    forecast_results: dict,
    shmr_params: SHMRParams | None = None,
    save_path: str | Path | None = None,
) -> None:
    """
    Two-regime summary: z=0 shape constraints (left) vs z-evolution (right).

    Left panel shows all surveys on the 5 z=0 shape + scatter parameters.
    Right panel shows only surveys that constrain evolution (9-param Fisher).

    Parameters
    ----------
    forecast_results : dict
        Output of run_forecast().
    shmr_params : SHMRParams or None
    save_path : str, Path, or None
    """
    if shmr_params is None:
        shmr_params = SHMRParams()

    z0_params = ["log_M1_0", "N_0", "beta_0", "gamma_0", "sigma_logMs"]
    evo_params = ["nu_M1", "nu_N", "nu_beta", "nu_gamma"]

    survey_colors = {
        "stage3_shallow_wide": "#1f77b4",
        "stage4_low_z": "#ff7f0e",
        "stage4_high_z": "#9467bd",
        "stage5_wide": "#2ca02c",
        "stage5_deep": "#d62728",
    }

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # --- Left: z=0 shape parameters ---
    ax = axes[0]
    x = np.arange(len(z0_params))
    width = 0.15
    i = 0
    for sname, r in forecast_results.items():
        color = survey_colors.get(sname, "gray")
        errs = []
        for p in z0_params:
            e = r.get("shmr_errors", r.get("errors", {})).get(p, 0)
            fid = abs(getattr(shmr_params, p))
            errs.append(e / fid * 100 if fid > 1e-10 else 0)
        ax.bar(x + i * width, errs, width, label=r["survey_name"],
               color=color)
        i += 1

    ax.set_xlabel("SHMR Parameter", fontsize=12)
    ax.set_ylabel(r"$\sigma / |\theta_\mathrm{fid}|$ [%]", fontsize=12)
    ax.set_title("z=0 SHMR shape + scatter (all surveys)", fontsize=13)
    ax.set_xticks(x + width * (i - 1) / 2)
    ax.set_xticklabels([PARAM_LATEX.get(p, p) for p in z0_params], fontsize=11)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, alpha=0.3, axis="y")
    ax.set_yscale("log")

    # --- Right: evolution parameters ---
    ax = axes[1]
    x = np.arange(len(evo_params))
    i = 0
    for sname, r in forecast_results.items():
        # Only plot surveys that have evolution parameters
        has_evo = any(p in r.get("errors", {}) for p in evo_params)
        if not has_evo:
            continue
        color = survey_colors.get(sname, "gray")
        errs = []
        for p in evo_params:
            e = r.get("shmr_errors", r.get("errors", {})).get(p, 0)
            fid = abs(getattr(shmr_params, p))
            errs.append(e / fid * 100 if fid > 1e-10 else 0)
        ax.bar(x + i * width, errs, width, label=r["survey_name"],
               color=color)
        i += 1

    ax.set_xlabel("Evolution Parameter", fontsize=12)
    ax.set_ylabel(r"$\sigma / |\theta_\mathrm{fid}|$ [%]", fontsize=12)
    ax.set_title("Redshift evolution (wide-z surveys only)", fontsize=13)
    if i > 0:
        ax.set_xticks(x + width * (i - 1) / 2)
    else:
        ax.set_xticks(x)
    ax.set_xticklabels([PARAM_LATEX.get(p, p) for p in evo_params], fontsize=11)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, alpha=0.3, axis="y")
    ax.set_yscale("log")

    plt.tight_layout()
    _save_or_show(fig, save_path)


def plot_improvement_factor(
    results_baseline: dict,
    results_target: dict,
    baseline_label: str = "Stage-IV",
    target_label: str = "Stage-V",
    save_path: str | Path | None = None,
) -> None:
    """
    Improvement factor: ratio of baseline to target errors per parameter.

    Parameters
    ----------
    results_baseline : dict
        Single survey result dict (containing 'shmr_errors').
    results_target : dict
        Single survey result dict.
    baseline_label, target_label : str
        Labels for the bar chart.
    save_path : str, Path, or None
    """
    # Find shared parameters
    base_errs = results_baseline.get("shmr_errors", results_baseline.get("errors", {}))
    tgt_errs = results_target.get("shmr_errors", results_target.get("errors", {}))
    shared = [p for p in base_errs if p in tgt_errs]

    ratios = []
    for p in shared:
        if tgt_errs[p] > 0:
            ratios.append(base_errs[p] / tgt_errs[p])
        else:
            ratios.append(0)

    fig, ax = plt.subplots(figsize=(10, 5))

    x = np.arange(len(shared))
    bars = ax.bar(x, ratios, color="#2ca02c", edgecolor="black", lw=0.8)

    # Reference line at 1
    ax.axhline(1.0, color="gray", ls="--", lw=1.5, alpha=0.7)

    ax.set_xlabel("SHMR Parameter", fontsize=13)
    ax.set_ylabel(
        rf"$\sigma({baseline_label}) / \sigma({target_label})$", fontsize=13)
    ax.set_title(
        f"Improvement factor: {target_label} over {baseline_label}",
        fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels([PARAM_LATEX.get(p, p) for p in shared], fontsize=11)
    ax.grid(True, alpha=0.3, axis="y")

    # Annotate bars
    for bar, r in zip(bars, ratios):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.05,
                f"{r:.1f}x", ha="center", va="bottom", fontsize=10,
                fontweight="bold")

    plt.tight_layout()
    _save_or_show(fig, save_path)
