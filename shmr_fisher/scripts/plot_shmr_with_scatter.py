"""
Generate a visualization of the fiducial SHMR at z=0, 0.5, 1.0
with mass-dependent scatter (Cao & Tinker 2020).

The main panel shows mean log M*(Mh) with ±1σ bands where σ varies
with halo mass. A secondary panel shows M*/Mh (star formation efficiency).

Output: outputs/phase1/shmr_with_scatter.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import shmr_fisher  # noqa: F401 — sets colossus cosmology
from shmr_fisher.config import SHMRParams
from shmr_fisher.shmr_model import mean_log_Mstar, scatter_at_Mh


def plot_shmr_with_scatter() -> None:
    """
    Plot the SHMR at three redshifts with mass-dependent scatter envelopes.

    Top panel: log M* vs log Mh with ±1σ(Mh) shaded bands.
    Bottom panel: M*/Mh (star formation efficiency) vs log Mh.
    """
    params = SHMRParams(use_mass_dependent_scatter=True)
    log_Mh = np.linspace(10.0, 15.0, 500)

    redshifts = [0.0, 0.5, 1.0]
    colors = ['C0', 'C1', 'C3']

    # Mass-dependent scatter (same at all z)
    sigma = scatter_at_Mh(log_Mh, params)

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(8, 8), sharex=True,
        gridspec_kw={'height_ratios': [3, 2], 'hspace': 0.05},
    )

    # --- Top panel: log M* vs log Mh ---
    for z, color in zip(redshifts, colors):
        log_ms = mean_log_Mstar(log_Mh, params, z)

        # ±1σ envelope from mass-dependent scatter
        ax1.fill_between(
            log_Mh, log_ms - sigma, log_ms + sigma,
            alpha=0.15, color=color,
        )
        ax1.plot(
            log_Mh, log_ms, color=color, lw=2.2,
            label=f'$z = {z:.1f}$',
        )

    # 1:1 line for reference (M* = Mh, i.e., 100% efficiency)
    ax1.plot(
        [10, 15], [10, 15], 'k:', lw=0.8, alpha=0.4,
        label=r'$M_* = M_h$ (100% eff.)',
    )

    ax1.set_ylabel(r'$\log_{10}(M_* / M_\odot)$', fontsize=13)
    ax1.set_ylim(6.5, 12.5)
    ax1.legend(loc='upper left', fontsize=11, framealpha=0.9)
    ax1.set_title(
        'Fiducial SHMR (Moster+2013) with Mass-Dependent Scatter (Cao & Tinker 2020)',
        fontsize=12,
    )
    ax1.grid(True, alpha=0.3)

    # Annotate that scatter bands are ±1σ(Mh)
    ax1.annotate(
        r'shaded: $\pm 1\sigma(M_h)$',
        xy=(10.5, 7.5), fontsize=10, style='italic', color='gray',
    )

    # --- Bottom panel: star formation efficiency M*/Mh ---
    for z, color in zip(redshifts, colors):
        log_ms = mean_log_Mstar(log_Mh, params, z)
        # M*/Mh = 10^(log_ms - log_Mh)
        efficiency = 10.0 ** (log_ms - log_Mh)

        # ±1σ envelope in linear efficiency space
        eff_hi = 10.0 ** (log_ms + sigma - log_Mh)
        eff_lo = 10.0 ** (log_ms - sigma - log_Mh)

        ax2.fill_between(log_Mh, eff_lo, eff_hi, alpha=0.15, color=color)
        ax2.plot(log_Mh, efficiency, color=color, lw=2.2)

    ax2.set_xlabel(r'$\log_{10}(M_h / M_\odot)$', fontsize=13)
    ax2.set_ylabel(r'$M_* / M_h$', fontsize=13)
    ax2.set_yscale('log')
    ax2.set_ylim(1e-5, 0.2)
    ax2.set_xlim(10.0, 15.0)
    ax2.grid(True, alpha=0.3, which='both')

    # Mark the peak efficiency band
    ax2.axhspan(0.02, 0.06, color='gray', alpha=0.08)
    ax2.text(
        14.5, 0.035, 'Peak\nefficiency\n(2–6%)',
        fontsize=8, ha='center', va='center', color='gray', style='italic',
    )

    fig.tight_layout()

    out_path = Path(__file__).resolve().parent.parent / 'outputs' / 'phase1' / 'shmr_with_scatter.png'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {out_path}")
    plt.close(fig)


if __name__ == '__main__':
    plot_shmr_with_scatter()
