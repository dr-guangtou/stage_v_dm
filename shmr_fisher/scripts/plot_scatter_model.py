"""
Generate a visualization of the mass-dependent SHMR scatter model.

Plots sigma(log M* | Mh) as a function of halo mass for the Cao & Tinker
(2020) parameterization, compared to the constant-scatter baseline.

Output: outputs/phase1/scatter_vs_log_Mh.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Add parent directory to path so we can import the package
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import shmr_fisher  # noqa: F401 — sets colossus cosmology
from shmr_fisher.config import SHMRParams
from shmr_fisher.shmr_model import scatter_at_Mh


def plot_scatter_vs_log_Mh() -> None:
    """
    Plot sigma(log M* | Mh) vs. log Mh for constant and mass-dependent models.

    Shows:
    - Constant scatter (sigma_logMs = 0.15 dex), the current baseline.
    - Cao & Tinker (2020) fiducial mass-dependent scatter.
    - Shaded band showing +/- 20% variation in scatter_sigma_rise.
    - Annotations for the physical regimes (dwarf vs. massive galaxies).
    """
    log_Mh = np.linspace(10.0, 15.0, 500)

    # Constant scatter model
    params_const = SHMRParams(use_mass_dependent_scatter=False)
    sigma_const = scatter_at_Mh(log_Mh, params_const)

    # Mass-dependent scatter: fiducial Cao & Tinker (2020)
    params_md = SHMRParams(use_mass_dependent_scatter=True)
    sigma_md = scatter_at_Mh(log_Mh, params_md)

    # Variation band: +/- 30% in scatter_sigma_rise
    params_md_lo = SHMRParams(
        use_mass_dependent_scatter=True,
        scatter_sigma_rise=0.07,
    )
    params_md_hi = SHMRParams(
        use_mass_dependent_scatter=True,
        scatter_sigma_rise=0.13,
    )
    sigma_lo = scatter_at_Mh(log_Mh, params_md_lo)
    sigma_hi = scatter_at_Mh(log_Mh, params_md_hi)

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(8, 5))

    # Constant scatter
    ax.axhline(
        params_const.sigma_logMs, color='gray', ls='--', lw=1.5,
        label=f'Constant ($\\sigma = {params_const.sigma_logMs}$ dex)',
    )

    # Mass-dependent scatter with uncertainty band
    ax.fill_between(
        log_Mh, sigma_lo, sigma_hi,
        alpha=0.25, color='C0',
        label=r'$\sigma_{\rm rise}$ variation ($\pm 30\%$)',
    )
    ax.plot(
        log_Mh, sigma_md, color='C0', lw=2.5,
        label='Cao & Tinker (2020) fiducial',
    )

    # Mark the transition mass
    ax.axvline(
        params_md.scatter_log_Mh_break, color='C3', ls=':', lw=1.2, alpha=0.7,
    )
    ax.annotate(
        r'$\log M_{h,\mathrm{break}}$',
        xy=(params_md.scatter_log_Mh_break, 0.145),
        xytext=(params_md.scatter_log_Mh_break + 0.3, 0.14),
        fontsize=10, color='C3',
        arrowprops=dict(arrowstyle='->', color='C3', lw=1.2),
    )

    # Annotate asymptotic values
    ax.annotate(
        rf'$\sigma_{{\rm high}} = {params_md.scatter_sigma_high}$ dex',
        xy=(14.5, params_md.scatter_sigma_high),
        xytext=(13.8, params_md.scatter_sigma_high + 0.04),
        fontsize=10, color='C0',
        arrowprops=dict(arrowstyle='->', color='C0', lw=1.0),
    )
    low_mh_sigma = params_md.scatter_sigma_high + 2 * params_md.scatter_sigma_rise
    ax.annotate(
        rf'$\sigma_{{\rm high}} + 2\sigma_{{\rm rise}} = {low_mh_sigma}$ dex',
        xy=(10.3, low_mh_sigma),
        xytext=(10.8, low_mh_sigma + 0.03),
        fontsize=10, color='C0',
        arrowprops=dict(arrowstyle='->', color='C0', lw=1.0),
    )

    # Physical regime annotations
    ax.text(
        10.5, 0.10, 'Dwarf galaxies\n(low-mass halos)',
        fontsize=9, ha='center', style='italic', color='gray',
    )
    ax.text(
        14.0, 0.10, 'Clusters\n(high-mass halos)',
        fontsize=9, ha='center', style='italic', color='gray',
    )

    ax.set_xlabel(r'$\log_{10}(M_h / M_\odot)$', fontsize=13)
    ax.set_ylabel(
        r'$\sigma(\log M_* \,|\, M_h)$ [dex]', fontsize=13,
    )
    ax.set_title(
        'SHMR Scatter: Constant vs. Mass-Dependent (Cao & Tinker 2020)',
        fontsize=12,
    )
    ax.set_xlim(10.0, 15.0)
    ax.set_ylim(0.05, 0.50)
    ax.legend(loc='upper right', fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()

    out_path = Path(__file__).resolve().parent.parent / 'outputs' / 'phase1' / 'scatter_vs_log_Mh.png'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {out_path}")
    plt.close(fig)


if __name__ == '__main__':
    # Print a quick diagnostic table
    print("--- Scatter model diagnostic ---")
    params = SHMRParams(use_mass_dependent_scatter=True)
    test_masses = [10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 14.0]
    print(f"{'log Mh':>8s}  {'sigma [dex]':>12s}")
    for lmh in test_masses:
        s = scatter_at_Mh(lmh, params)
        print(f"{lmh:8.1f}  {s:12.4f}")

    # Generate the figure
    plot_scatter_vs_log_Mh()
