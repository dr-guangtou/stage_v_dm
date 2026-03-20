"""Regenerate Phase 2 QA figures with mass-dependent scatter."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from shmr_fisher.config import SHMRParams, ForecastConfig, LensingConfig
from shmr_fisher.halo_model import (
    n_cen, n_sat, galaxy_number_density, effective_bias, delta_sigma_bin,
)
from shmr_fisher.covariance import lensing_covariance
from shmr_fisher.survey_configs import surveys
from shmr_fisher.shmr_model import mean_log_Mstar

params = SHMRParams(use_mass_dependent_scatter=True)
survey = surveys["stage4_low_z"]
fc = ForecastConfig()
outdir = Path("outputs") / "phase2"
outdir.mkdir(parents=True, exist_ok=True)
z = 0.3
R_Mpc = np.logspace(np.log10(fc.R_min_Mpc), np.log10(fc.R_max_Mpc), 30)
bins = [(9.5, 10.0), (10.0, 10.5), (10.5, 11.0)]

# --- hod_occupation.png ---
log_Mh = np.linspace(10, 15, 300)
fig, ax = plt.subplots(figsize=(8, 6))
for lo, hi in bins:
    ncen = n_cen(log_Mh, lo, hi, params, z)
    nsat = n_sat(log_Mh, lo, hi, params, z)
    ax.plot(log_Mh, ncen, label=f"N_cen [{lo:.1f},{hi:.1f}]")
    ax.plot(log_Mh, nsat, "--", label=f"N_sat [{lo:.1f},{hi:.1f}]")
ax.set_xlabel(r"log$_{10}$(M$_h$/M$_\odot$)")
ax.set_ylabel("Occupation number")
ax.set_yscale("log")
ax.set_ylim(1e-4, 100)
ax.legend(fontsize=8)
ax.set_title(f"HOD occupation at z={z}")
fig.tight_layout()
fig.savefig(outdir / "hod_occupation.png", dpi=300)
plt.close(fig)
print("hod_occupation.png saved")

# --- covariance_diagnostic.png ---
R_cov = np.logspace(np.log10(fc.R_min_Mpc), np.log10(fc.R_max_Mpc), 10)
lo, hi = 10.0, 10.5
ds_cov = delta_sigma_bin(R_cov, lo, hi, params, z)
lensing = LensingConfig()
n_gal = galaxy_number_density(lo, hi, params, z)
# N_lens = n_gal * survey volume (approximate for illustration)
from colossus.cosmology import cosmology
cosmo = cosmology.getCurrent()
# Comoving volume for the survey area and a thin redshift shell
area_sr = survey.area_deg2 * (np.pi / 180) ** 2
d_c = cosmo.comovingDistance(z_max=z) / cosmo.H0 * 100  # Mpc
dz = 0.1
d_c_lo = cosmo.comovingDistance(z_max=z - dz / 2) / cosmo.H0 * 100
d_c_hi = cosmo.comovingDistance(z_max=z + dz / 2) / cosmo.H0 * 100
vol = area_sr / 3 * (d_c_hi**3 - d_c_lo**3)  # Mpc^3
N_lens = n_gal * vol
# R_cov should be bin edges for lensing_covariance
R_edges = np.logspace(np.log10(fc.R_min_Mpc), np.log10(fc.R_max_Mpc), 11)
R_mid = np.sqrt(R_edges[:-1] * R_edges[1:])  # geometric mean
ds_cov = delta_sigma_bin(R_mid, lo, hi, params, z)
var_ds = lensing_covariance(R_edges, z, N_lens, survey.area_deg2, lensing)
diag_err = np.sqrt(var_ds)  # lensing_covariance returns 1D variance array

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
ax1.loglog(R_Mpc, delta_sigma_bin(R_Mpc, lo, hi, params, z) / 1e12, "b-")
ax1.errorbar(
    R_mid, ds_cov / 1e12, yerr=diag_err / 1e12,
    fmt="ro", capsize=3, label="with errors",
)
ax1.set_xlabel("R [Mpc]")
ax1.set_ylabel(r"$\Delta\Sigma$ [M$_\odot$/pc$^2$]")
ax1.set_title(f"Signal + errors, z={z}, log M*=[{lo},{hi}]")
ax1.legend()

# Diagonal covariance -> identity correlation matrix
# Instead, show S/N per radial bin
snr = ds_cov / diag_err
ax2.semilogx(R_mid, snr, "ko-")
ax2.set_xlabel("R [Mpc]")
ax2.set_ylabel("S/N per bin")
ax2.set_title(f"S/N per radial bin (N_lens = {N_lens:.1e})")
ax2.axhline(1, color="r", ls="--", alpha=0.5)
fig.tight_layout()
fig.savefig(outdir / "covariance_diagnostic.png", dpi=300)
plt.close(fig)
print("covariance_diagnostic.png saved")

# --- phase2_summary.png ---
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

log_Mh_arr = np.linspace(10, 15, 200)
ms = np.array([mean_log_Mstar(m, params, z) for m in log_Mh_arr])
axes[0, 0].plot(log_Mh_arr, ms)
axes[0, 0].set_xlabel("log Mh")
axes[0, 0].set_ylabel("log M*")
axes[0, 0].set_title(f"SHMR at z={z}")

for lo, hi in bins:
    ds = delta_sigma_bin(R_Mpc, lo, hi, params, z)
    axes[0, 1].loglog(R_Mpc, ds / 1e12, label=f"[{lo},{hi}]")
axes[0, 1].set_xlabel("R [Mpc]")
axes[0, 1].set_ylabel(r"$\Delta\Sigma$ [M$_\odot$/pc$^2$]")
axes[0, 1].legend()
axes[0, 1].set_title("Lensing signals")

ngs = [galaxy_number_density(lo, hi, params, z) for lo, hi in bins]
labels = [f"[{lo},{hi}]" for lo, hi in bins]
axes[1, 0].bar(labels, ngs)
axes[1, 0].set_ylabel("n_gal [Mpc$^{-3}$]")
axes[1, 0].set_title("Number densities")
axes[1, 0].set_yscale("log")

bs = [effective_bias(lo, hi, params, z) for lo, hi in bins]
axes[1, 1].bar(labels, bs)
axes[1, 1].set_ylabel("b_eff")
axes[1, 1].set_title("Effective bias")

fig.suptitle(f"Phase 2 Summary - z={z}", fontsize=14)
fig.tight_layout()
fig.savefig(outdir / "phase2_summary.png", dpi=300)
plt.close(fig)
print("phase2_summary.png saved")
