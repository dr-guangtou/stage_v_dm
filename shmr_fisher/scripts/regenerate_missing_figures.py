"""Generate missing sweep figures and systematics comparison."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from shmr_fisher.config import (
    SHMRParams, ForecastConfig, LensingConfig, NuisanceConfig,
)
from shmr_fisher.fisher import (
    compute_fisher_matrix, marginalized_errors, extract_shmr_constraints,
)
from shmr_fisher.survey_configs import surveys
from scripts.run_sweep import parameter_sweep, plot_sweep

params = SHMRParams(use_mass_dependent_scatter=True)
fc_sweep = ForecastConfig(n_R_bins=8, n_Mh_bins=150)
outdir_phase3 = Path("outputs") / "phase3"
outdir_phase3.mkdir(parents=True, exist_ok=True)

# --- sweep_n_gal_total_stage4_low_z.png ---
print("--- Sweep: n_gal_total ---")
r_ngal = parameter_sweep(
    "stage4_low_z", "n_gal_total",
    [1e6, 5e6, 1e7, 5e7, 1e8],
    forecast_config=fc_sweep, shmr_params=params,
)
plot_sweep(r_ngal)
print("sweep_n_gal_total done")

# --- sweep_dlog_Mstar_stage5_wide.png ---
print("\n--- Sweep: dlog_Mstar ---")
r_dlog = parameter_sweep(
    "stage5_wide", "dlog_Mstar",
    [1.0, 0.5, 0.25],
    forecast_config=fc_sweep, shmr_params=params,
)
plot_sweep(r_dlog)
print("sweep_dlog_Mstar done")

# --- systematics_comparison.png ---
print("\n--- Systematics comparison ---")
lensing = LensingConfig()
nc = NuisanceConfig()
survey = surveys["stage4_low_z"]

configs = {
    "No systematics": ForecastConfig(
        n_R_bins=10, n_Mh_bins=200,
        systematic_floor_fraction=0.0, include_nuisance_params=False,
    ),
    "Floor only (5%)": ForecastConfig(
        n_R_bins=10, n_Mh_bins=200,
        systematic_floor_fraction=0.05, include_nuisance_params=False,
    ),
    "Floor + nuisance": ForecastConfig(
        n_R_bins=10, n_Mh_bins=200,
        systematic_floor_fraction=0.05, include_nuisance_params=True,
    ),
}

all_errors = {}
shmr_pnames = None

for label, fc_sys in configs.items():
    print(f"  {label}...")
    nuis = nc if fc_sys.include_nuisance_params else None
    fisher, pnames, meta = compute_fisher_matrix(params, survey, lensing, fc_sys, nuis)
    eigv = np.linalg.eigvalsh(fisher)
    if np.all(eigv > 0):
        shmr_errs, shmr_names = extract_shmr_constraints(fisher, pnames)
        all_errors[label] = dict(zip(shmr_names, shmr_errs))
        if shmr_pnames is None:
            shmr_pnames = shmr_names
    else:
        print(f"    WARNING: not positive definite")

fig, ax = plt.subplots(figsize=(12, 6))
x = np.arange(len(shmr_pnames))
width = 0.8 / len(all_errors)
colors = ["#2ecc71", "#e67e22", "#e74c3c"]

for i, (label, errs) in enumerate(all_errors.items()):
    frac = [errs[p] / abs(getattr(params, p)) * 100 for p in shmr_pnames]
    ax.bar(x + i * width, frac, width, label=label, color=colors[i])

ax.set_xticks(x + width * (len(all_errors) - 1) / 2)
ax.set_xticklabels(shmr_pnames, rotation=45, ha="right")
ax.set_ylabel("Fractional error [%]")
ax.set_title("Impact of systematics on SHMR constraints (Stage-IV)")
ax.legend()
ax.set_yscale("log")
fig.tight_layout()
fig.savefig(outdir_phase3 / "systematics_comparison.png", dpi=300)
plt.close(fig)
print("systematics_comparison.png saved")

print("\nAll missing figures generated.")
