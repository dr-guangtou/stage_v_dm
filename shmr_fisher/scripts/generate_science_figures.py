"""
Generate all Phase 4 science figures.

Produces:
    outputs/delta_sigma_with_errors.png
    outputs/two_regime_summary.png
    outputs/improvement_factor.png
    outputs/sweep_area_deg2_stage4_low_z.png
    outputs/sweep_log_Mstar_min_stage5_wide.png
    outputs/sweep_z_max_stage4_low_z.png
"""

import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np

from shmr_fisher.config import (
    ForecastConfig, LensingConfig, NuisanceConfig, SHMRParams,
)
from shmr_fisher.plot_results import (
    plot_delta_sigma_with_errors,
    plot_two_regime_summary,
    plot_improvement_factor,
)
from scripts.run_forecast import run_forecast
from scripts.run_sweep import parameter_sweep, plot_sweep

params = SHMRParams()
outdir = Path("outputs")

# Run full forecast with systematics for science figures
print("="*60)
print("Running full forecast (with systematics)...")
print("="*60)

fc = ForecastConfig(
    n_R_bins=10, n_Mh_bins=200,
    systematic_floor_fraction=0.05, include_nuisance_params=True,
)
nc = NuisanceConfig()

results = run_forecast(
    forecast_config=fc, nuisance_config=nc, output_dir="outputs",
)

# -------------------------------------------------------------------
# Figure: DeltaSigma with error bars
# -------------------------------------------------------------------
print("\n--- Generating DeltaSigma with error bars ---")
plot_delta_sigma_with_errors(
    results, params, z=0.3, mstar_bin=(10.5, 11.0),
    save_path=outdir / "delta_sigma_with_errors.png",
)

# -------------------------------------------------------------------
# Figure: Two-regime summary
# -------------------------------------------------------------------
print("\n--- Generating two-regime summary ---")
plot_two_regime_summary(
    results, params,
    save_path=outdir / "two_regime_summary.png",
)

# -------------------------------------------------------------------
# Figure: Improvement factor (Stage-V Wide over Stage-IV Low-z)
# -------------------------------------------------------------------
print("\n--- Generating improvement factor ---")
# Use 5-param Fisher for fair comparison
fc_5p = ForecastConfig(
    n_R_bins=8, n_Mh_bins=150, vary_z_evolution=False,
    systematic_floor_fraction=0.05, include_nuisance_params=True,
)
results_5p = run_forecast(
    survey_names=["stage4_low_z", "stage5_wide"],
    forecast_config=fc_5p, nuisance_config=nc, output_dir="outputs",
)
plot_improvement_factor(
    results_5p["stage4_low_z"], results_5p["stage5_wide"],
    baseline_label="Stage-IV", target_label="Stage-V",
    save_path=outdir / "improvement_factor.png",
)

# -------------------------------------------------------------------
# Remaining sweeps
# -------------------------------------------------------------------
fc_sweep = ForecastConfig(n_R_bins=8, n_Mh_bins=150)

print("\n--- Sweep: area_deg2 ---")
r_area = parameter_sweep(
    "stage4_low_z", "area_deg2",
    [1000, 3000, 7000, 10000, 14000],
    forecast_config=fc_sweep,
)
plot_sweep(r_area)

print("\n--- Sweep: log_Mstar_min ---")
r_mmin = parameter_sweep(
    "stage5_wide", "log_Mstar_min",
    [9.5, 9.0, 8.5, 8.0, 7.5],
    forecast_config=fc_sweep,
)
plot_sweep(r_mmin)

print("\n--- Sweep: z_max ---")
r_zmax = parameter_sweep(
    "stage4_low_z", "z_max",
    [0.2, 0.3, 0.4, 0.6, 0.8, 1.0],
    forecast_config=fc_sweep,
)
plot_sweep(r_zmax)

# -------------------------------------------------------------------
# Append captions
# -------------------------------------------------------------------
timestamp = datetime.now().strftime("%Y-%m-%d %H:%M PDT")

captions = f"""
---

### `delta_sigma_with_errors.png`
**Created:** {timestamp}

Galaxy-galaxy lensing signal DeltaSigma(R) for stellar mass bin [10.5, 11.0] at z=0.3, with measurement error bars for three survey tiers: Stage-III (blue), Stage-IV (orange), Stage-V (green). The black solid line shows the halo-model prediction. Error bars include shape noise and a 5% systematic floor. Radial bins are log-spaced from 0.1 to 30 physical Mpc.

The error bars shrink from Stage-III to Stage-V primarily due to larger lens samples (more lens-source pairs reduce shape noise). At small R (< 0.5 Mpc), errors are smallest because the annular area is small but the source density per unit Mpc^2 is high. At large R (> 10 Mpc), both signal and S/N drop rapidly.

**Purpose:** Visualizes the measurement precision achievable for a representative stellar mass bin across survey generations. This is the raw input to the Fisher matrix — the error bars directly determine the constraining power.

---

### `two_regime_summary.png`
**Created:** {timestamp}

Two-panel side-by-side comparison highlighting the complementarity of different survey strategies.

**Left panel — z=0 SHMR shape + scatter (all 5 surveys):** Bar chart of marginalized fractional errors on log M1, N_0, beta_0, gamma_0, sigma_logMs. All five predefined surveys are shown. Low-z surveys (Stage-III, Stage-IV Low-z) constrain the z=0 shape well because their galaxy samples sit on the SHMR at z~0. Stage-V Wide has the tightest z=0 constraints because it combines large volume with deep mass completeness.

**Right panel — Redshift evolution (wide-z surveys only):** Bar chart of fractional errors on nu_M1, nu_N, nu_beta, nu_gamma. Only surveys that span enough redshift range to constrain evolution are shown (Stage-IV High-z, Stage-V Wide, Stage-V Deep). Stage-V Wide provides the tightest evolution constraints thanks to its wide z-range (0.05–1.0) and large galaxy count. Stage-IV High-z has the weakest evolution constraints due to its narrow mass range (log M*_min=10.8), which limits the number of informative M* bins.

**Purpose:** Main science summary figure demonstrating that low-z depth constrains the z=0 SHMR shape, while wide z-coverage constrains its evolution — motivating the design of next-generation spectroscopic surveys.

---

### `improvement_factor.png`
**Created:** {timestamp}

Bar chart showing the improvement factor sigma(Stage-IV) / sigma(Stage-V) for each of the 5 z=0 SHMR parameters, using 5-parameter Fisher matrices with systematics (5% floor + nuisance marginalization) for fair comparison. Values > 1 mean Stage-V is better.

Each bar is annotated with the numerical improvement factor. Parameters most improved by Stage-V are those most sensitive to the survey's advantages (larger galaxy count, deeper mass completeness, wider z-range). The systematic floor limits the maximum achievable improvement (compared to stat-only, where Stage-V gains would be larger).

**Purpose:** The main "science case" figure — quantifies how much better a Stage-V survey constrains the galaxy-halo connection compared to current Stage-IV capabilities.
"""

with open(outdir / "CAPTION.md", "a") as f:
    f.write(captions)
print(f"\nCaptions appended to {outdir / 'CAPTION.md'}")

print("\n" + "="*60)
print("All science figures generated.")
print("="*60)
