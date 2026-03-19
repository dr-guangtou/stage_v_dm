"""Plotting helpers for benchmark figures."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


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
