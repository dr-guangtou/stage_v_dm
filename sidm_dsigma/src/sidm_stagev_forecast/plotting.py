"""Plotting helpers for benchmark figures."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def plot_rho_profiles(
    r_kpc: np.ndarray,
    rho_cdm_msun_kpc3: np.ndarray,
    sidm_profiles: dict[float, np.ndarray],
    output_path: Path,
    title: str,
) -> None:
    """Save a slide-ready rho(r) comparison plot."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

    axes[0].loglog(r_kpc, rho_cdm_msun_kpc3, color="black", lw=2.2, label="CDM")
    for sigma_over_m, rho_values in sidm_profiles.items():
        axes[0].loglog(r_kpc, rho_values, lw=1.7, label=f"SIDM {sigma_over_m:g} cm^2/g")
    axes[0].set_xlabel("r [kpc]")
    axes[0].set_ylabel(r"$\rho(r)$ [Msun / kpc$^3$]")
    axes[0].set_title(f"{title}: 3D Density")
    axes[0].legend(fontsize=8)

    axes[1].axhline(1.0, color="black", lw=1.0, ls="--")
    for sigma_over_m, rho_values in sidm_profiles.items():
        ratio = rho_values / rho_cdm_msun_kpc3
        axes[1].semilogx(r_kpc, ratio, lw=1.7, label=f"{sigma_over_m:g}")
    axes[1].set_xlabel("r [kpc]")
    axes[1].set_ylabel(r"$\rho_{\mathrm{SIDM}}/\rho_{\mathrm{CDM}}$")
    axes[1].set_title(f"{title}: Ratio to CDM")

    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def plot_delta_sigma_profiles(
    r_projected_kpc: np.ndarray,
    delta_sigma_cdm_msun_kpc2: np.ndarray,
    sidm_delta_sigma_profiles: dict[float, np.ndarray],
    fractional_error: np.ndarray,
    output_path: Path,
    title: str,
) -> None:
    """Save a slide-ready DeltaSigma and ratio figure."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

    axes[0].loglog(
        r_projected_kpc,
        delta_sigma_cdm_msun_kpc2,
        color="black",
        lw=2.2,
        label="CDM",
    )

    cdm_error = fractional_error * np.abs(delta_sigma_cdm_msun_kpc2)
    axes[0].fill_between(
        r_projected_kpc,
        delta_sigma_cdm_msun_kpc2 - cdm_error,
        delta_sigma_cdm_msun_kpc2 + cdm_error,
        color="gray",
        alpha=0.2,
        label=r"Toy 1$\sigma$",
    )

    for sigma_over_m, delta_sigma_values in sidm_delta_sigma_profiles.items():
        axes[0].loglog(
            r_projected_kpc,
            delta_sigma_values,
            lw=1.7,
            label=f"SIDM {sigma_over_m:g} cm^2/g",
        )

    axes[0].set_xlabel("R [kpc]")
    axes[0].set_ylabel(r"$\Delta\Sigma$ [Msun / kpc$^2$]")
    axes[0].set_title(f"{title}: Lensing Signal")
    axes[0].legend(fontsize=8)

    axes[1].axhline(1.0, color="black", lw=1.0, ls="--")
    for sigma_over_m, delta_sigma_values in sidm_delta_sigma_profiles.items():
        ratio = delta_sigma_values / delta_sigma_cdm_msun_kpc2
        axes[1].semilogx(r_projected_kpc, ratio, lw=1.7, label=f"{sigma_over_m:g}")
    axes[1].set_xlabel("R [kpc]")
    axes[1].set_ylabel(r"$\Delta\Sigma_{\mathrm{SIDM}}/\Delta\Sigma_{\mathrm{CDM}}$")
    axes[1].set_title(f"{title}: Ratio to CDM")

    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def plot_stacked_delta_sigma(
    r_projected_kpc: np.ndarray,
    delta_sigma_cdm_msun_kpc2: np.ndarray,
    sidm_delta_sigma_profiles: dict[float, np.ndarray],
    output_path: Path,
) -> None:
    """Save stacked DeltaSigma and ratio-to-CDM panels."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

    axes[0].loglog(r_projected_kpc, delta_sigma_cdm_msun_kpc2, color="black", lw=2.2, label="CDM")
    for sigma_over_m, values in sidm_delta_sigma_profiles.items():
        axes[0].loglog(r_projected_kpc, values, lw=1.8, label=f"SIDM {sigma_over_m:g} cm^2/g")
    axes[0].set_xlabel("R [kpc]")
    axes[0].set_ylabel(r"Stacked $\Delta\Sigma$ [Msun / kpc$^2$]")
    axes[0].set_title("Tier-1 Ensemble: Stacked Lensing")
    axes[0].legend(fontsize=8)

    axes[1].axhline(1.0, color="black", lw=1.0, ls="--")
    for sigma_over_m, values in sidm_delta_sigma_profiles.items():
        axes[1].semilogx(
            r_projected_kpc,
            values / delta_sigma_cdm_msun_kpc2,
            lw=1.8,
            label=f"{sigma_over_m:g}",
        )
    axes[1].set_xlabel("R [kpc]")
    axes[1].set_ylabel(r"$\Delta\Sigma_{\mathrm{SIDM}} / \Delta\Sigma_{\mathrm{CDM}}$")
    axes[1].set_title("Ratio to CDM")

    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, quantile: float) -> float:
    """Return the weighted quantile for one-dimensional input arrays."""
    if values.ndim != 1 or weights.ndim != 1 or values.size != weights.size:
        raise ValueError("values and weights must be one-dimensional with equal length.")
    if values.size == 0:
        raise ValueError("values cannot be empty.")
    if not 0.0 <= quantile <= 1.0:
        raise ValueError("quantile must be in [0, 1].")

    order = np.argsort(values)
    sorted_values = values[order]
    sorted_weights = np.clip(weights[order], 0.0, np.inf)
    total_weight = float(np.sum(sorted_weights))
    if not np.isfinite(total_weight) or total_weight <= 0.0:
        sorted_weights = np.ones_like(sorted_values, dtype=float)
        total_weight = float(sorted_weights.size)
    cumulative_weight = np.cumsum(sorted_weights)
    target = quantile * total_weight
    return float(np.interp(target, cumulative_weight, sorted_values))


def _weighted_summary(values: np.ndarray, weights: np.ndarray) -> tuple[float, float, float]:
    """Return weighted mean, p16, and p84 for one-dimensional arrays."""
    clipped_weights = np.clip(np.asarray(weights, dtype=float), 0.0, np.inf)
    if not np.isfinite(clipped_weights).all() or float(np.sum(clipped_weights)) <= 0.0:
        clipped_weights = np.ones_like(values, dtype=float)
    weighted_mean = float(np.average(values, weights=clipped_weights))
    p16 = _weighted_quantile(values, clipped_weights, 0.16)
    p84 = _weighted_quantile(values, clipped_weights, 0.84)
    return weighted_mean, p16, p84


def plot_tier1_summary(
    r_3d_kpc: np.ndarray,
    rho_cdm_msun_kpc3: np.ndarray,
    rho_sidm_profiles: dict[float, np.ndarray],
    r_projected_kpc: np.ndarray,
    delta_sigma_cdm_msun_kpc2: np.ndarray,
    delta_sigma_sidm_profiles: dict[float, np.ndarray],
    sigma_over_m_grid_cm2_g: tuple[float, ...],
    delta_chi2_by_sigma: dict[float, float],
    output_path: Path,
    h: float,
    m200_msun_samples: np.ndarray,
    c200_samples: np.ndarray,
    weights: np.ndarray,
    annotation_text: str | None = None,
) -> None:
    """Save a 4-panel ensemble summary figure."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)

    axis_rho = axes[0, 0]
    axis_delta_sigma = axes[0, 1]
    axis_ratio = axes[1, 0]
    axis_chi2 = axes[1, 1]
    mass_inset = inset_axes(axis_rho, width="50%", height="33%", loc="lower left", borderpad=0.8)
    concentration_inset = inset_axes(
        axis_delta_sigma,
        width="50%",
        height="33%",
        loc="lower center",
        borderpad=0.8,
    )

    rho_scale = (h**2) / 1.0e9
    delta_sigma_scale = h / 1.0e6
    rho_cdm_h2_msun_pc3 = rho_cdm_msun_kpc3 * rho_scale
    delta_sigma_cdm_h_msun_pc2 = delta_sigma_cdm_msun_kpc2 * delta_sigma_scale

    axis_rho.loglog(r_3d_kpc, rho_cdm_h2_msun_pc3, color="black", lw=2.2, label="CDM")
    for sigma_over_m in sigma_over_m_grid_cm2_g:
        axis_rho.loglog(
            r_3d_kpc,
            rho_sidm_profiles[sigma_over_m] * rho_scale,
            lw=1.6,
            label=f"SIDM {sigma_over_m:g}",
        )
    axis_rho.set_xlabel("r [kpc]")
    axis_rho.set_ylabel(r"Stacked $\rho$ [$h^2$ Msun / pc$^3$]")
    axis_rho.set_title("Ensemble 3D Density")
    axis_rho.legend(fontsize=7)

    axis_delta_sigma.loglog(
        r_projected_kpc,
        delta_sigma_cdm_h_msun_pc2,
        color="black",
        lw=2.2,
        label="CDM",
    )
    for sigma_over_m in sigma_over_m_grid_cm2_g:
        axis_delta_sigma.loglog(
            r_projected_kpc,
            delta_sigma_sidm_profiles[sigma_over_m] * delta_sigma_scale,
            lw=1.6,
            label=f"{sigma_over_m:g}",
        )
    axis_delta_sigma.set_xlabel("R [kpc]")
    axis_delta_sigma.set_ylabel(r"Stacked $\Delta\Sigma$ [$h$ Msun / pc$^2$]")
    axis_delta_sigma.set_title("Ensemble Lensing Signal")

    axis_ratio.axhline(1.0, color="black", lw=1.0, ls="--")
    for sigma_over_m in sigma_over_m_grid_cm2_g:
        axis_ratio.semilogx(
            r_projected_kpc,
            delta_sigma_sidm_profiles[sigma_over_m] / delta_sigma_cdm_msun_kpc2,
            lw=1.6,
            label=f"{sigma_over_m:g}",
        )
    axis_ratio.set_xlabel("R [kpc]")
    axis_ratio.set_ylabel(r"$\Delta\Sigma_{\mathrm{SIDM}} / \Delta\Sigma_{\mathrm{CDM}}$")
    axis_ratio.set_title("Stacked Ratio")

    chi2_x = np.asarray(sigma_over_m_grid_cm2_g, dtype=float)
    chi2_y = np.asarray([delta_chi2_by_sigma[value] for value in sigma_over_m_grid_cm2_g])
    axis_chi2.plot(chi2_x, chi2_y, marker="o", lw=2.0, color="tab:red")
    axis_chi2.set_xlabel(r"$\sigma/m$ [cm$^2$/g]")
    axis_chi2.set_ylabel(r"$\Delta\chi^2$ (5% fractional error)")
    axis_chi2.set_title("Distinguishability vs SIDM Strength")
    axis_chi2.grid(alpha=0.25)

    mass_values = np.log10(np.clip(np.asarray(m200_msun_samples, dtype=float), 1.0, np.inf))
    concentration_values = np.asarray(c200_samples, dtype=float)
    histogram_weights = np.asarray(weights, dtype=float)

    mass_mean, mass_p16, mass_p84 = _weighted_summary(mass_values, histogram_weights)
    concentration_mean, concentration_p16, concentration_p84 = _weighted_summary(
        concentration_values, histogram_weights
    )

    mass_inset.hist(
        mass_values,
        bins=20,
        weights=histogram_weights,
        color="tab:blue",
        alpha=0.18,
        edgecolor="tab:blue",
        lw=0.8,
    )
    mass_inset.axvspan(mass_p16, mass_p84, color="tab:orange", alpha=0.24)
    mass_inset.axvline(mass_mean, color="tab:red", lw=1.2, ls="--")
    mass_inset.set_title(r"$\log_{10}(M_{200c}/M_\odot)$", fontsize=7)
    mass_inset.tick_params(labelsize=6)
    mass_inset.grid(alpha=0.2, ls=":")
    mass_inset.text(
        0.03,
        0.95,
        f"mean={mass_mean:.2f}\n16-84=[{mass_p16:.2f}, {mass_p84:.2f}]",
        transform=mass_inset.transAxes,
        va="top",
        fontsize=6,
    )

    concentration_inset.hist(
        concentration_values,
        bins=20,
        weights=histogram_weights,
        color="tab:green",
        alpha=0.18,
        edgecolor="tab:green",
        lw=0.8,
    )
    concentration_inset.axvspan(concentration_p16, concentration_p84, color="tab:orange", alpha=0.24)
    concentration_inset.axvline(concentration_mean, color="tab:red", lw=1.2, ls="--")
    concentration_inset.set_title(r"$c_{200c}$", fontsize=7)
    concentration_inset.tick_params(labelsize=6)
    concentration_inset.grid(alpha=0.2, ls=":")
    concentration_inset.text(
        0.03,
        0.95,
        (
            f"mean={concentration_mean:.2f}\n"
            f"16-84=[{concentration_p16:.2f}, {concentration_p84:.2f}]"
        ),
        transform=concentration_inset.transAxes,
        va="top",
        fontsize=6,
    )

    if annotation_text:
        fig.text(
            0.02,
            0.02,
            annotation_text,
            ha="left",
            va="bottom",
            fontsize=8,
            bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.85, "edgecolor": "gray"},
        )

    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def log_slope(r_kpc: np.ndarray, rho_msun_kpc3: np.ndarray) -> np.ndarray:
    """Return dln(rho)/dln(r) on the input radial grid."""
    return np.gradient(np.log(rho_msun_kpc3), np.log(r_kpc))


def plot_single_halo_tier2_diagnostics(
    r_kpc: np.ndarray,
    rho_cdm_inner: np.ndarray,
    rho_sidm_inner: np.ndarray,
    rho_dk14_reference: np.ndarray,
    rho_hybrid: np.ndarray,
    r_projected_kpc: np.ndarray,
    delta_sigma_tier1: np.ndarray,
    delta_sigma_tier2: np.ndarray,
    output_path: Path,
    title: str,
) -> None:
    """Save Tier-2 single-halo diagnostics for density, slope, and DeltaSigma."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)

    axes[0].loglog(r_kpc, rho_cdm_inner, color="black", lw=2.0, label="CDM inner (NFW)")
    axes[0].loglog(r_kpc, rho_sidm_inner, color="tab:blue", lw=1.8, label="SIDM inner")
    axes[0].loglog(r_kpc, rho_dk14_reference, color="tab:orange", lw=1.8, label="DK14-like ref")
    axes[0].loglog(r_kpc, rho_hybrid, color="tab:red", lw=2.0, ls="--", label="SIDM+DK14 hybrid")
    axes[0].set_xlabel("r [kpc]")
    axes[0].set_ylabel(r"$\rho$ [Msun / kpc$^3$]")
    axes[0].set_title(f"{title}: Density Components")
    axes[0].legend(fontsize=8)

    axes[1].semilogx(r_kpc, log_slope(r_kpc, rho_cdm_inner), color="black", lw=2.0, label="CDM inner")
    axes[1].semilogx(r_kpc, log_slope(r_kpc, rho_sidm_inner), color="tab:blue", lw=1.8, label="SIDM inner")
    axes[1].semilogx(
        r_kpc,
        log_slope(r_kpc, rho_dk14_reference),
        color="tab:orange",
        lw=1.8,
        label="DK14-like ref",
    )
    axes[1].semilogx(r_kpc, log_slope(r_kpc, rho_hybrid), color="tab:red", lw=2.0, ls="--", label="Hybrid")
    axes[1].set_xlabel("r [kpc]")
    axes[1].set_ylabel(r"$d\ln \rho / d\ln r$")
    axes[1].set_title("Log-Slope Diagnostic")
    axes[1].axhline(-3.0, color="gray", lw=0.9, ls=":")

    axes[2].loglog(r_projected_kpc, delta_sigma_tier1, color="tab:blue", lw=2.0, label="Tier-1 inner-only")
    axes[2].loglog(r_projected_kpc, delta_sigma_tier2, color="tab:red", lw=2.0, ls="--", label="Tier-2 hybrid")
    axes[2].set_xlabel("R [kpc]")
    axes[2].set_ylabel(r"$\Delta\Sigma$ [Msun / kpc$^2$]")
    axes[2].set_title("Projection Stability")
    axes[2].legend(fontsize=8)

    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def plot_tier1_tier2_stacked_comparison(
    r_projected_kpc: np.ndarray,
    cdm_tier1: np.ndarray,
    sidm_tier1_by_sigma: dict[float, np.ndarray],
    cdm_tier2: np.ndarray,
    sidm_tier2_by_sigma: dict[float, np.ndarray],
    output_path: Path,
    title: str,
) -> None:
    """Save stacked Tier-1 vs Tier-2 DeltaSigma comparison and ratio panels."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)

    axes[0].loglog(r_projected_kpc, cdm_tier1, color="black", lw=2.2, label="CDM Tier-1")
    axes[0].loglog(r_projected_kpc, cdm_tier2, color="black", lw=1.8, ls="--", label="CDM Tier-2")
    for sigma_over_m, values in sidm_tier1_by_sigma.items():
        axes[0].loglog(
            r_projected_kpc,
            values,
            lw=1.4,
            alpha=0.7,
            label=f"SIDM {sigma_over_m:g} Tier-1",
        )
        axes[0].loglog(
            r_projected_kpc,
            sidm_tier2_by_sigma[sigma_over_m],
            lw=1.8,
            ls="--",
            label=f"SIDM {sigma_over_m:g} Tier-2",
        )
    axes[0].set_xlabel("R [kpc]")
    axes[0].set_ylabel(r"Stacked $\Delta\Sigma$ [Msun / kpc$^2$]")
    axes[0].set_title(f"{title}: Tier-1 vs Tier-2")
    axes[0].legend(fontsize=7, ncol=2)

    axes[1].axhline(1.0, color="black", lw=1.0, ls="--")
    axes[1].semilogx(r_projected_kpc, cdm_tier2 / cdm_tier1, color="black", lw=2.0, label="CDM T2/T1")
    for sigma_over_m, values in sidm_tier1_by_sigma.items():
        axes[1].semilogx(
            r_projected_kpc,
            sidm_tier2_by_sigma[sigma_over_m] / values,
            lw=1.8,
            label=f"SIDM {sigma_over_m:g} T2/T1",
        )
    axes[1].set_xlabel("R [kpc]")
    axes[1].set_ylabel(r"$\Delta\Sigma_{\mathrm{Tier2}} / \Delta\Sigma_{\mathrm{Tier1}}$")
    axes[1].set_title("Outskirts Attachment Impact")
    axes[1].legend(fontsize=8)

    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def plot_delta_chi2_tier_comparison(
    sigma_over_m_grid_cm2_g: tuple[float, ...],
    delta_chi2_tier1: dict[float, float],
    delta_chi2_tier2: dict[float, float],
    output_path: Path,
    title: str,
    label_tier1: str = "Tier-1 inner-only",
    label_tier2: str = "Tier-2 hybrid",
) -> None:
    """Save DeltaChi2 versus sigma/m comparison between Tier-1 and Tier-2."""
    x_values = np.asarray(sigma_over_m_grid_cm2_g, dtype=float)
    y_tier1 = np.asarray([delta_chi2_tier1[value] for value in sigma_over_m_grid_cm2_g])
    y_tier2 = np.asarray([delta_chi2_tier2[value] for value in sigma_over_m_grid_cm2_g])

    fig, axis = plt.subplots(1, 1, figsize=(6.5, 4.5), constrained_layout=True)
    axis.plot(x_values, y_tier1, marker="o", lw=2.0, label=label_tier1)
    axis.plot(x_values, y_tier2, marker="s", lw=2.0, ls="--", label=label_tier2)
    axis.set_xlabel(r"$\sigma/m$ [cm$^2$/g]")
    axis.set_ylabel(r"$\Delta\chi^2$ (5% fractional error)")
    axis.set_title(f"{title}: Distinguishability")
    axis.grid(alpha=0.25)
    axis.legend()
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def plot_single_halo_tier3_diagnostics(
    r_kpc: np.ndarray,
    rho_tier2: np.ndarray,
    rho_tier3: np.ndarray,
    r_projected_kpc: np.ndarray,
    sigma_tier2: np.ndarray,
    sigma_tier3: np.ndarray,
    delta_sigma_tier2: np.ndarray,
    delta_sigma_tier3: np.ndarray,
    output_path: Path,
    title: str,
) -> None:
    """Save Tier-3 single-halo diagnostics in density, slope, Sigma, and DeltaSigma."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)

    axes[0, 0].loglog(r_kpc, rho_tier2, color="tab:blue", lw=2.0, label="Tier-2")
    axes[0, 0].loglog(r_kpc, rho_tier3, color="tab:red", lw=2.0, ls="--", label="Tier-3")
    axes[0, 0].set_xlabel("r [kpc]")
    axes[0, 0].set_ylabel(r"$\rho$ [Msun / kpc$^3$]")
    axes[0, 0].set_title(f"{title}: 3D Density")
    axes[0, 0].legend()

    axes[0, 1].semilogx(r_kpc, log_slope(r_kpc, rho_tier2), color="tab:blue", lw=2.0, label="Tier-2")
    axes[0, 1].semilogx(r_kpc, log_slope(r_kpc, rho_tier3), color="tab:red", lw=2.0, ls="--", label="Tier-3")
    axes[0, 1].set_xlabel("r [kpc]")
    axes[0, 1].set_ylabel(r"$d\ln \rho / d\ln r$")
    axes[0, 1].set_title("Slope")
    axes[0, 1].axhline(-3.0, color="gray", lw=0.9, ls=":")

    axes[1, 0].loglog(r_projected_kpc, sigma_tier2, color="tab:blue", lw=2.0, label="Tier-2")
    axes[1, 0].loglog(r_projected_kpc, sigma_tier3, color="tab:red", lw=2.0, ls="--", label="Tier-3")
    axes[1, 0].set_xlabel("R [kpc]")
    axes[1, 0].set_ylabel(r"$\Sigma$ [Msun / kpc$^2$]")
    axes[1, 0].set_title("Projected Surface Density")
    axes[1, 0].legend()

    axes[1, 1].loglog(r_projected_kpc, delta_sigma_tier2, color="tab:blue", lw=2.0, label="Tier-2")
    axes[1, 1].loglog(r_projected_kpc, delta_sigma_tier3, color="tab:red", lw=2.0, ls="--", label="Tier-3")
    axes[1, 1].set_xlabel("R [kpc]")
    axes[1, 1].set_ylabel(r"$\Delta\Sigma$ [Msun / kpc$^2$]")
    axes[1, 1].set_title("Excess Surface Density")
    axes[1, 1].legend()

    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def plot_tier2_tier3_stacked_comparison(
    r_projected_kpc: np.ndarray,
    cdm_tier2: np.ndarray,
    sidm_tier2_by_sigma: dict[float, np.ndarray],
    cdm_tier3: np.ndarray,
    sidm_tier3_by_sigma: dict[float, np.ndarray],
    output_path: Path,
    title: str,
) -> None:
    """Save stacked Tier-2 vs Tier-3 comparison with ratio panels."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)

    axes[0].loglog(r_projected_kpc, cdm_tier2, color="black", lw=2.2, label="CDM Tier-2")
    axes[0].loglog(r_projected_kpc, cdm_tier3, color="black", lw=1.8, ls="--", label="CDM Tier-3")
    for sigma_over_m, values in sidm_tier2_by_sigma.items():
        axes[0].loglog(r_projected_kpc, values, lw=1.4, alpha=0.8, label=f"SIDM {sigma_over_m:g} T2")
        axes[0].loglog(
            r_projected_kpc,
            sidm_tier3_by_sigma[sigma_over_m],
            lw=1.8,
            ls="--",
            label=f"SIDM {sigma_over_m:g} T3",
        )
    axes[0].set_xlabel("R [kpc]")
    axes[0].set_ylabel(r"Stacked $\Delta\Sigma$ [Msun / kpc$^2$]")
    axes[0].set_title(f"{title}: Tier-2 vs Tier-3")
    axes[0].legend(fontsize=7, ncol=2)

    axes[1].axhline(1.0, color="black", lw=1.0, ls="--")
    axes[1].semilogx(r_projected_kpc, cdm_tier2 / cdm_tier2, color="black", lw=2.0, label="CDM T2")
    axes[1].semilogx(r_projected_kpc, cdm_tier3 / cdm_tier2, color="black", lw=1.8, ls="--", label="CDM T3")
    for sigma_over_m, values in sidm_tier2_by_sigma.items():
        axes[1].semilogx(r_projected_kpc, values / cdm_tier2, lw=1.5, label=f"SIDM {sigma_over_m:g} / CDM T2")
        axes[1].semilogx(
            r_projected_kpc,
            sidm_tier3_by_sigma[sigma_over_m] / cdm_tier3,
            lw=1.8,
            ls="--",
            label=f"SIDM {sigma_over_m:g} / CDM T3",
        )
    axes[1].set_xlabel("R [kpc]")
    axes[1].set_ylabel(r"SIDM / CDM")
    axes[1].set_title("SIDM-to-CDM Ratios")
    axes[1].legend(fontsize=7, ncol=1)

    axes[2].axhline(1.0, color="black", lw=1.0, ls="--")
    axes[2].semilogx(r_projected_kpc, cdm_tier3 / cdm_tier2, color="black", lw=2.0, label="CDM T3/T2")
    for sigma_over_m, values in sidm_tier2_by_sigma.items():
        axes[2].semilogx(
            r_projected_kpc,
            sidm_tier3_by_sigma[sigma_over_m] / values,
            lw=1.8,
            label=f"SIDM {sigma_over_m:g} T3/T2",
        )
    axes[2].set_xlabel("R [kpc]")
    axes[2].set_ylabel(r"Tier-3 / Tier-2")
    axes[2].set_title("Empirical Outer Correction Impact")
    axes[2].legend(fontsize=7)

    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def plot_tier3_sensitivity_panel(
    r_projected_kpc: np.ndarray,
    baseline_ratio: np.ndarray,
    varied_ratios: dict[str, np.ndarray],
    output_path: Path,
    title: str,
) -> None:
    """Plot one-at-a-time Tier-3 nuisance sensitivity in outskirts ratios."""
    fig, axis = plt.subplots(1, 1, figsize=(6.8, 4.5), constrained_layout=True)
    axis.axhline(1.0, color="black", lw=1.0, ls="--")
    axis.semilogx(r_projected_kpc, baseline_ratio, lw=2.2, color="black", label="baseline")
    for label, values in varied_ratios.items():
        axis.semilogx(r_projected_kpc, values, lw=1.8, label=label)
    axis.set_xlabel("R [kpc]")
    axis.set_ylabel(r"$\Delta\Sigma_{\mathrm{Tier3}} / \Delta\Sigma_{\mathrm{Tier2}}$")
    axis.set_title(title)
    axis.legend(fontsize=8)
    axis.grid(alpha=0.2)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def plot_tier3_redshift_overlay_summary(
    redshift_datasets: list[dict[str, object]],
    sigma_over_m_grid_cm2_g: tuple[float, ...],
    output_path: Path,
    annotation_text: str,
    h: float,
) -> None:
    """Save a 4-panel ensemble summary with multiple redshifts overplotted."""
    if len(redshift_datasets) == 0:
        raise ValueError("redshift_datasets cannot be empty.")

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    axis_rho = axes[0, 0]
    axis_delta_sigma = axes[0, 1]
    axis_ratio = axes[1, 0]
    axis_chi2 = axes[1, 1]
    mass_inset = inset_axes(
        axis_rho,
        width="54%",
        height="34%",
        loc="lower center",
        borderpad=0.8,
    )
    concentration_inset = inset_axes(
        axis_delta_sigma,
        width="54%",
        height="34%",
        loc="lower center",
        borderpad=0.8,
    )

    style_cycle = ["-", "--", "-.", ":"]
    sigma_colors = {sigma: plt.cm.viridis(index / max(len(sigma_over_m_grid_cm2_g) - 1, 1)) for index, sigma in enumerate(sigma_over_m_grid_cm2_g)}
    redshift_colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]

    for index, dataset in enumerate(redshift_datasets):
        label = str(dataset["label"])
        linestyle = style_cycle[index % len(style_cycle)]
        redshift_color = redshift_colors[index % len(redshift_colors)]
        r_3d_kpc = np.asarray(dataset["r_3d_kpc"], dtype=float)
        rho_cdm_h2_msun_pc3 = np.asarray(dataset["rho_cdm_msun_kpc3"], dtype=float) * ((h**2) / 1.0e9)
        r_projected_kpc = np.asarray(dataset["r_projected_kpc"], dtype=float)
        ds_cdm_msun_kpc2 = np.asarray(dataset["delta_sigma_cdm_msun_kpc2"], dtype=float)
        ds_cdm_h_msun_pc2 = ds_cdm_msun_kpc2 * (h / 1.0e6)

        axis_rho.loglog(r_3d_kpc, rho_cdm_h2_msun_pc3, color="black", lw=2.0, ls=linestyle)
        axis_delta_sigma.loglog(r_projected_kpc, ds_cdm_h_msun_pc2, color="black", lw=2.0, ls=linestyle)

        sidm_rho_profiles = dataset["rho_sidm_profiles"]
        sidm_delta_sigma_profiles = dataset["delta_sigma_sidm_profiles"]
        for sigma in sigma_over_m_grid_cm2_g:
            color = sigma_colors[sigma]
            sidm_rho_h2_msun_pc3 = np.asarray(sidm_rho_profiles[sigma], dtype=float) * ((h**2) / 1.0e9)
            sidm_ds_h_msun_pc2 = np.asarray(sidm_delta_sigma_profiles[sigma], dtype=float) * (h / 1.0e6)

            axis_rho.loglog(r_3d_kpc, sidm_rho_h2_msun_pc3, color=color, lw=1.4, ls=linestyle)
            axis_delta_sigma.loglog(r_projected_kpc, sidm_ds_h_msun_pc2, color=color, lw=1.4, ls=linestyle)
            axis_ratio.semilogx(
                r_projected_kpc,
                np.asarray(sidm_delta_sigma_profiles[sigma], dtype=float) / ds_cdm_msun_kpc2,
                color=color,
                lw=1.4,
                ls=linestyle,
            )

        delta_chi2_by_sigma = dataset["delta_chi2_by_sigma"]
        chi2_x = np.asarray(sigma_over_m_grid_cm2_g, dtype=float)
        chi2_y = np.asarray([delta_chi2_by_sigma[sigma] for sigma in sigma_over_m_grid_cm2_g], dtype=float)
        axis_chi2.plot(chi2_x, chi2_y, lw=2.0, ls=linestyle, marker="o", label=label)

        if (
            "m200_msun_samples" in dataset
            and "c200_samples" in dataset
            and "weights" in dataset
        ):
            masses = np.asarray(dataset["m200_msun_samples"], dtype=float)
            concentrations = np.asarray(dataset["c200_samples"], dtype=float)
            weights = np.asarray(dataset["weights"], dtype=float)
            log_mass = np.log10(np.clip(masses, 1.0, np.inf))

            mass_inset.hist(
                log_mass,
                bins=20,
                weights=weights,
                histtype="step",
                color=redshift_color,
                ls=linestyle,
                lw=1.2,
                alpha=0.9,
            )
            concentration_inset.hist(
                concentrations,
                bins=20,
                weights=weights,
                histtype="step",
                color=redshift_color,
                ls=linestyle,
                lw=1.2,
                alpha=0.9,
            )
            mass_mean, mass_p16, mass_p84 = _weighted_summary(log_mass, weights)
            concentration_mean, concentration_p16, concentration_p84 = _weighted_summary(
                concentrations, weights
            )
            mass_inset.axvspan(mass_p16, mass_p84, color=redshift_color, alpha=0.08)
            mass_inset.axvline(mass_mean, color=redshift_color, lw=1.0, ls=linestyle)
            concentration_inset.axvspan(concentration_p16, concentration_p84, color=redshift_color, alpha=0.08)
            concentration_inset.axvline(concentration_mean, color=redshift_color, lw=1.0, ls=linestyle)

    axis_rho.set_xlabel("r [kpc]")
    axis_rho.set_ylabel(r"Stacked $\rho$ [$h^2$ Msun / pc$^3$]")
    axis_rho.set_title("Ensemble 3D Density")

    axis_delta_sigma.set_xlabel("R [kpc]")
    axis_delta_sigma.set_ylabel(r"Stacked $\Delta\Sigma$ [$h$ Msun / pc$^2$]")
    axis_delta_sigma.set_title("Ensemble Lensing Signal")

    axis_ratio.axhline(1.0, color="black", lw=1.0, ls="--")
    axis_ratio.set_xlabel("R [kpc]")
    axis_ratio.set_ylabel(r"$\Delta\Sigma_{\mathrm{SIDM}} / \Delta\Sigma_{\mathrm{CDM}}$")
    axis_ratio.set_title("Stacked Ratio")

    axis_chi2.set_xlabel(r"$\sigma/m$ [cm$^2$/g]")
    axis_chi2.set_ylabel(r"$\Delta\chi^2$ (5% fractional error)")
    axis_chi2.set_title("Distinguishability vs SIDM Strength")
    axis_chi2.grid(alpha=0.25)
    axis_chi2.legend(fontsize=8)

    mass_inset.set_title(r"$\log_{10} M_{200c}$", fontsize=7)
    mass_inset.tick_params(labelsize=6)
    mass_inset.grid(alpha=0.2, ls=":")

    concentration_inset.set_title(r"$c_{200c}$", fontsize=7)
    concentration_inset.tick_params(labelsize=6)
    concentration_inset.grid(alpha=0.2, ls=":")

    sigma_handles = [
        Line2D([0], [0], color=sigma_colors[sigma], lw=2.0, label=f"SIDM {sigma:g}")
        for sigma in sigma_over_m_grid_cm2_g
    ]
    redshift_handles = [
        Line2D([0], [0], color="gray", lw=2.0, ls=style_cycle[index % len(style_cycle)], label=str(dataset["label"]))
        for index, dataset in enumerate(redshift_datasets)
    ]
    cdm_handle = [Line2D([0], [0], color="black", lw=2.0, label="CDM")]
    axis_rho.legend(handles=cdm_handle + sigma_handles + redshift_handles, fontsize=7, loc="best")

    fig.text(
        0.02,
        0.02,
        annotation_text,
        ha="left",
        va="bottom",
        fontsize=8,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.85, "edgecolor": "gray"},
    )
    fig.savefig(output_path, dpi=220)
    plt.close(fig)
