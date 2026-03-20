"""
Generate all Phase 4 science figures.

Produces:
    outputs/phase4/delta_sigma_with_errors.png
    outputs/phase4/two_regime_summary.png
    outputs/phase4/improvement_factor.png
    outputs/phase4/sweep_area_deg2_stage4_low_z.png
    outputs/phase4/sweep_log_Mstar_min_stage5_wide.png
    outputs/phase4/sweep_z_max_stage4_low_z.png
"""

import sys
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

# Use mass-dependent scatter (Cao & Tinker 2020)
params = SHMRParams(use_mass_dependent_scatter=True)
outdir = Path("outputs") / "phase4"
outdir.mkdir(parents=True, exist_ok=True)

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
    forecast_config=fc, nuisance_config=nc, shmr_params=params,
    output_dir="outputs",
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
    forecast_config=fc_5p, nuisance_config=nc, shmr_params=params,
    output_dir="outputs",
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
    forecast_config=fc_sweep, shmr_params=params,
)
plot_sweep(r_area)

print("\n--- Sweep: log_Mstar_min ---")
r_mmin = parameter_sweep(
    "stage5_wide", "log_Mstar_min",
    [9.5, 9.0, 8.5, 8.0, 7.5],
    forecast_config=fc_sweep, shmr_params=params,
)
plot_sweep(r_mmin)

print("\n--- Sweep: z_max ---")
r_zmax = parameter_sweep(
    "stage4_low_z", "z_max",
    [0.2, 0.3, 0.4, 0.6, 0.8, 1.0],
    forecast_config=fc_sweep, shmr_params=params,
)
plot_sweep(r_zmax)

# Captions are maintained manually in outputs/CAPTION.md

print("\n" + "="*60)
print("All science figures generated.")
print("="*60)
