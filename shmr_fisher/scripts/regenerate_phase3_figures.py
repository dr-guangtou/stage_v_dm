"""Regenerate Phase 3 Fisher figures with mass-dependent scatter."""

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

params = SHMRParams(use_mass_dependent_scatter=True)
outdir = Path("outputs") / "phase3"
outdir.mkdir(parents=True, exist_ok=True)

# Run Fisher for all surveys with systematics
fc = ForecastConfig(
    n_R_bins=10, n_Mh_bins=200,
    systematic_floor_fraction=0.05,
    include_nuisance_params=True,
)
nc = NuisanceConfig()
lensing = LensingConfig()

results = {}
for name, survey in surveys.items():
    print(f"\nRunning: {name}")
    nuis = nc if fc.include_nuisance_params else None
    fisher, param_names, metadata = compute_fisher_matrix(
        params, survey, lensing, fc, nuis,
    )
    eigvals = np.linalg.eigvalsh(fisher)
    pos_def = bool(np.all(eigvals > 0))
    if pos_def:
        errors = marginalized_errors(fisher)
        shmr_errors, shmr_names = extract_shmr_constraints(fisher, param_names)
    else:
        errors = np.full(len(param_names), np.inf)
        shmr_errors = errors
        shmr_names = param_names
        print(f"  WARNING: Not positive definite!")

    results[name] = {
        "survey_name": survey.name,
        "param_names": param_names,
        "shmr_param_names": shmr_names,
        "errors": dict(zip(param_names, errors.tolist())),
        "shmr_errors": dict(zip(shmr_names, shmr_errors.tolist())),
        "condition_number": metadata.get("condition_number", np.inf),
        "positive_definite": pos_def,
        "fisher": fisher,
    }
    print(f"  cond = {results[name]['condition_number']:.2e}, pos_def = {pos_def}")
    for pn, err in zip(shmr_names, shmr_errors):
        fid = abs(getattr(params, pn))
        frac = err / fid * 100 if fid > 1e-10 else 0
        print(f"    {pn:15s}: {err:.6f} ({frac:.2f}%)")

# --- fisher_errors_comparison.png ---
fig, ax = plt.subplots(figsize=(12, 6))
survey_names = list(results.keys())

# Use the union of all SHMR param names (sorted by first occurrence)
all_pnames = []
for sname in survey_names:
    for pn in results[sname]["shmr_param_names"]:
        if pn not in all_pnames:
            all_pnames.append(pn)

x = np.arange(len(all_pnames))
width = 0.8 / len(survey_names)

for i, sname in enumerate(survey_names):
    r = results[sname]
    frac_errors = []
    for pn in all_pnames:
        if pn in r["shmr_errors"]:
            fid = abs(getattr(params, pn))
            err = r["shmr_errors"][pn]
            frac_errors.append(err / fid * 100 if fid > 1e-10 else 0)
        else:
            frac_errors.append(0)  # parameter not constrained by this survey
    ax.bar(x + i * width, frac_errors, width, label=r["survey_name"])

ax.set_xticks(x + width * (len(survey_names) - 1) / 2)
ax.set_xticklabels(all_pnames, rotation=45, ha="right")
ax.set_ylabel("Fractional error [%]")
ax.set_title("SHMR parameter constraints (with systematics)")
ax.legend()
ax.set_yscale("log")
fig.tight_layout()
fig.savefig(outdir / "fisher_errors_comparison.png", dpi=300)
plt.close(fig)
print("\nfisher_errors_comparison.png saved")

# --- fisher_ellipses.png ---
# 2D confidence ellipses for a pair of interesting parameters
from matplotlib.patches import Ellipse


def plot_ellipse(ax, fisher_full, param_names, p1, p2, color, label):
    """Plot 1-sigma ellipse for parameters p1, p2."""
    i1 = param_names.index(p1)
    i2 = param_names.index(p2)
    cov_full = np.linalg.inv(fisher_full)
    cov_2d = np.array([
        [cov_full[i1, i1], cov_full[i1, i2]],
        [cov_full[i2, i1], cov_full[i2, i2]],
    ])
    eigvals, eigvecs = np.linalg.eigh(cov_2d)
    angle = np.degrees(np.arctan2(eigvecs[1, 0], eigvecs[0, 0]))
    w, h = 2 * np.sqrt(eigvals)  # 1-sigma
    ell = Ellipse(
        (getattr(params, p1), getattr(params, p2)),
        w, h, angle=angle, fill=False, color=color, lw=2, label=label,
    )
    ax.add_patch(ell)


# Pick two well-constrained SHMR parameters
p1, p2 = "log_M1", "log_epsilon"
fig, ax = plt.subplots(figsize=(8, 6))
colors = plt.cm.Set1(np.linspace(0, 1, len(survey_names)))

for i, sname in enumerate(survey_names):
    r = results[sname]
    if not r["positive_definite"]:
        continue
    try:
        plot_ellipse(ax, r["fisher"], r["param_names"], p1, p2, colors[i], r["survey_name"])
    except (ValueError, np.linalg.LinAlgError) as e:
        print(f"  Skipping ellipse for {sname}: {e}")

ax.set_xlabel(p1)
ax.set_ylabel(p2)
ax.set_title(f"Fisher ellipses: {p1} vs {p2}")
ax.legend()
ax.autoscale()
fig.tight_layout()
fig.savefig(outdir / "fisher_ellipses.png", dpi=300)
plt.close(fig)
print("fisher_ellipses.png saved")

# --- fisher_ngal_scaling.png ---
# Show how errors scale with N_gal (shape-noise dominated -> 1/sqrt(N))
fc_scaling = ForecastConfig(n_R_bins=8, n_Mh_bins=150, vary_z_evolution=False)
base_survey = surveys["stage4_low_z"]
from dataclasses import replace

ngal_factors = [0.25, 0.5, 1.0, 2.0, 4.0]
errors_vs_ngal = []

for factor in ngal_factors:
    modified = replace(base_survey, n_gal_total=base_survey.n_gal_total * factor)
    fisher_m, pnames_m, _ = compute_fisher_matrix(params, modified, lensing, fc_scaling, None)
    eigv = np.linalg.eigvalsh(fisher_m)
    if np.all(eigv > 0):
        errs = marginalized_errors(fisher_m)
        errors_vs_ngal.append(errs)
    else:
        errors_vs_ngal.append(np.full(len(pnames_m), np.inf))

errors_arr = np.array(errors_vs_ngal)  # (n_factors, n_params)
fig, ax = plt.subplots(figsize=(8, 6))
for j, pname in enumerate(pnames_m):
    ax.loglog(ngal_factors, errors_arr[:, j] / errors_arr[2, j], "o-", label=pname)
ax.loglog(ngal_factors, 1.0 / np.sqrt(ngal_factors), "k--", lw=2, label=r"$1/\sqrt{N}$")
ax.set_xlabel("N_gal multiplier")
ax.set_ylabel("Relative error")
ax.set_title("Error scaling with galaxy count")
ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig(outdir / "fisher_ngal_scaling.png", dpi=300)
plt.close(fig)
print("fisher_ngal_scaling.png saved")

print("\nAll Phase 3 figures generated.")
