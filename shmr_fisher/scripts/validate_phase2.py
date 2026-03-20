"""
Phase 2 QA: generate validation plots for halo model predictions.

Produces:
    outputs/phase2/hod_occupation.png    — N_cen and N_sat vs Mh for 4 M* bins
    outputs/phase2/delta_sigma_model.png — DeltaSigma(R) for 4 M* bins at z=0.3
    outputs/phase2/phase2_summary.png    — 4-panel diagnostic: HOD, DS, b_eff, n_gal vs z
"""

import sys
sys.path.insert(0, '.')

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from shmr_fisher.config import SHMRParams
from shmr_fisher.halo_model import (
    n_cen, n_sat, delta_sigma_nfw, delta_sigma_bin,
    effective_bias, galaxy_number_density,
)

params = SHMRParams()
outdir = Path('outputs') / 'phase2'
outdir.mkdir(parents=True, exist_ok=True)

bins = [(9.5, 10.0), (10.0, 10.5), (10.5, 11.0), (11.0, 11.5)]
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
log_Mh = np.linspace(10.0, 15.0, 300)

# -----------------------------------------------------------------------
# Plot 1: HOD occupation functions
# -----------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax = axes[0]
for (lo, hi), color in zip(bins, colors):
    nc = n_cen(log_Mh, lo, hi, params, 0.3)
    ax.plot(log_Mh, nc, color=color, lw=2,
            label=rf'$\log M_* \in [{lo},{hi}]$')
ax.set_xlabel(r'$\log_{10}(M_h / M_\odot)$', fontsize=13)
ax.set_ylabel(r'$\langle N_\mathrm{cen} \rangle$', fontsize=13)
ax.set_title('Central occupation at z=0.3', fontsize=13)
ax.legend(fontsize=10)
ax.set_xlim(10, 15)
ax.set_ylim(-0.02, 1.05)
ax.grid(True, alpha=0.3)

ax = axes[1]
for (lo, hi), color in zip(bins, colors):
    ns = n_sat(log_Mh, lo, hi, params, 0.3)
    ax.plot(log_Mh, ns, color=color, lw=2,
            label=rf'$\log M_* \in [{lo},{hi}]$')
ax.set_xlabel(r'$\log_{10}(M_h / M_\odot)$', fontsize=13)
ax.set_ylabel(r'$\langle N_\mathrm{sat} \rangle$', fontsize=13)
ax.set_title('Satellite occupation at z=0.3', fontsize=13)
ax.legend(fontsize=10)
ax.set_xlim(10, 15)
ax.set_yscale('log')
ax.set_ylim(1e-5, 100)
ax.grid(True, alpha=0.3)

plt.tight_layout()
fig.savefig(outdir / 'hod_occupation.png', dpi=300, bbox_inches='tight')
print(f'Saved: {outdir / "hod_occupation.png"}')
plt.close(fig)

# -----------------------------------------------------------------------
# Plot 2: DeltaSigma model — SPEC.md Section 2.3
# -----------------------------------------------------------------------
R_Mpc = np.logspace(np.log10(0.1), np.log10(30.0), 30)

fig, ax = plt.subplots(figsize=(8, 6))

for (lo, hi), color in zip(bins, colors):
    ds = delta_sigma_bin(R_Mpc, lo, hi, params, 0.3)
    ds_pc2 = ds / 1e12  # Msun/Mpc^2 -> Msun/pc^2
    ax.plot(R_Mpc, ds_pc2, color=color, lw=2,
            label=rf'$\log M_* \in [{lo},{hi}]$')

# Also show single NFW for reference
for Mh_ref, ls, lbl in [(1e12, '--', r'$M_h=10^{12}$'),
                          (1e13, ':', r'$M_h=10^{13}$')]:
    ds_single = delta_sigma_nfw(R_Mpc, Mh_ref, 0.3) / 1e12
    ax.plot(R_Mpc, ds_single, ls=ls, color='gray', lw=1.5, alpha=0.7,
            label=f'Single NFW {lbl}')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$R$ [physical Mpc]', fontsize=13)
ax.set_ylabel(r'$\Delta\Sigma$ [$M_\odot/\mathrm{pc}^2$]', fontsize=13)
ax.set_title(r'Galaxy-galaxy lensing signal at $z=0.3$', fontsize=13)
ax.set_xlim(0.1, 30)
ax.set_ylim(0.01, 500)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3, which='both')

plt.tight_layout()
fig.savefig(outdir / 'delta_sigma_model.png', dpi=300, bbox_inches='tight')
print(f'Saved: {outdir / "delta_sigma_model.png"}')
plt.close(fig)

# -----------------------------------------------------------------------
# Plot 3: 4-panel summary diagnostic
# -----------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(13, 10))

# Panel A: DeltaSigma at R=1 Mpc vs M* bin, at multiple z
ax = axes[0, 0]
R_1Mpc = np.array([1.0])
z_list = [0.1, 0.3, 0.5, 1.0]
z_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
bin_centers = [0.5 * (lo + hi) for lo, hi in bins]

for z_val, zc in zip(z_list, z_colors):
    ds_vals = []
    for lo, hi in bins:
        ds = delta_sigma_bin(R_1Mpc, lo, hi, params, z_val)
        ds_vals.append(ds[0] / 1e12)
    ax.plot(bin_centers, ds_vals, 'o-', color=zc, lw=2, ms=6,
            label=f'z = {z_val}')

ax.set_xlabel(r'$\log_{10}(M_* / M_\odot)$ bin center', fontsize=12)
ax.set_ylabel(r'$\Delta\Sigma(R=1\,\mathrm{Mpc})$ [$M_\odot/\mathrm{pc}^2$]', fontsize=12)
ax.set_title('Lensing signal vs stellar mass', fontsize=12)
ax.set_yscale('log')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Panel B: Effective bias vs M* bin at multiple z
ax = axes[0, 1]
for z_val, zc in zip(z_list, z_colors):
    b_vals = [effective_bias(lo, hi, params, z_val) for lo, hi in bins]
    ax.plot(bin_centers, b_vals, 'o-', color=zc, lw=2, ms=6,
            label=f'z = {z_val}')

ax.set_xlabel(r'$\log_{10}(M_* / M_\odot)$ bin center', fontsize=12)
ax.set_ylabel(r'$b_\mathrm{eff}$', fontsize=12)
ax.set_title('Effective halo bias vs stellar mass', fontsize=12)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Panel C: Galaxy number density vs M* bin at multiple z
ax = axes[1, 0]
for z_val, zc in zip(z_list, z_colors):
    ng_vals = [galaxy_number_density(lo, hi, params, z_val) for lo, hi in bins]
    ax.plot(bin_centers, ng_vals, 'o-', color=zc, lw=2, ms=6,
            label=f'z = {z_val}')

ax.set_xlabel(r'$\log_{10}(M_* / M_\odot)$ bin center', fontsize=12)
ax.set_ylabel(r'$n_\mathrm{gal}$ [Mpc$^{-3}$]', fontsize=12)
ax.set_title('Galaxy number density vs stellar mass', fontsize=12)
ax.set_yscale('log')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Panel D: DeltaSigma(R) profile comparison: bin-averaged vs single NFW
ax = axes[1, 1]
R_plot = np.logspace(np.log10(0.1), np.log10(30.0), 25)
# Show bin [10.5, 11.0] at z=0.3 with the corresponding mean halo mass
from shmr_fisher.shmr_model import mean_log_Mh
log_Mstar_mid = 10.75
log_Mh_mean = mean_log_Mh(log_Mstar_mid, params, 0.3)
Mh_mean = 10**log_Mh_mean

ds_bin = delta_sigma_bin(R_plot, 10.5, 11.0, params, 0.3) / 1e12
ds_single = delta_sigma_nfw(R_plot, Mh_mean, 0.3) / 1e12

ax.plot(R_plot, ds_bin, 'b-', lw=2.5, label='Bin-averaged (HOD)')
ax.plot(R_plot, ds_single, 'r--', lw=2,
        label=rf'Single NFW ($M_h=10^{{{log_Mh_mean:.1f}}}$)')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$R$ [physical Mpc]', fontsize=12)
ax.set_ylabel(r'$\Delta\Sigma$ [$M_\odot/\mathrm{pc}^2$]', fontsize=12)
ax.set_title(r'$\log M_* \in [10.5, 11.0]$ at $z=0.3$', fontsize=12)
ax.set_xlim(0.1, 30)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3, which='both')

plt.tight_layout()
fig.savefig(outdir / 'phase2_summary.png', dpi=300, bbox_inches='tight')
print(f'Saved: {outdir / "phase2_summary.png"}')
plt.close(fig)

print('\nAll Phase 2 QA plots generated.')
