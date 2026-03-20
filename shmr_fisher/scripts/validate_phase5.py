"""
Phase 5 validation: sigma_obs and SMF covariance impact on forecasts.

Runs the dwarf regime forecast with and without observational scatter
(sigma_log_Mstar_obs), and with SMF (Poisson + cosmic variance) vs.
Poisson-only n_gal covariance. Prints comparison tables showing the
impact of each new feature on parameter constraints.

Usage:
    uv run python scripts/validate_phase5.py
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shmr_fisher.config import (
    ForecastConfig,
    LensingConfig,
    NuisanceConfig,
    SHMRParams,
    SurveyConfig,
)
from shmr_fisher.fisher import (
    compute_fisher_matrix,
    extract_shmr_constraints,
    marginalized_errors,
)


def make_dwarf_survey() -> SurveyConfig:
    """Create a simple dwarf-regime survey for testing."""
    return SurveyConfig(
        name="Dwarf Test",
        area_deg2=5000,
        z_min=0.01,
        z_max=0.15,
        n_gal_total=1_000_000,
        log_Mstar_min=7.5,
        log_Mstar_max=9.0,
        dlog_Mstar=0.5,
    )


def run_single(
    survey: SurveyConfig,
    shmr_params: SHMRParams,
    forecast_config: ForecastConfig,
    label: str,
) -> tuple[dict[str, float], list[str]]:
    """Run a single forecast and return SHMR errors."""
    lensing = LensingConfig()
    nuisance = NuisanceConfig() if forecast_config.include_nuisance_params else None

    fisher, param_names, metadata = compute_fisher_matrix(
        shmr_params, survey, lensing, forecast_config, nuisance,
    )

    eigvals = np.linalg.eigvalsh(fisher)
    if np.all(eigvals > 0):
        shmr_errors, shmr_names = extract_shmr_constraints(fisher, param_names)
        return dict(zip(shmr_names, shmr_errors)), shmr_names
    else:
        print(f"  WARNING: Fisher matrix not positive definite for '{label}'")
        return {n: np.inf for n in param_names}, list(param_names)


def print_comparison(
    label_a: str,
    errors_a: dict[str, float],
    label_b: str,
    errors_b: dict[str, float],
    shmr_params: SHMRParams,
) -> None:
    """Print a side-by-side comparison table."""
    print(f"\n{'='*75}")
    print(f"  {label_a}  vs.  {label_b}")
    print(f"{'='*75}")
    print(f"  {'Parameter':<20s} {'Fiducial':>10s} "
          f"{'sigma_A':>12s} {'sigma_B':>12s} {'Ratio B/A':>10s}")
    print(f"  {'-'*20} {'-'*10} {'-'*12} {'-'*12} {'-'*10}")

    for name in errors_a:
        if name not in errors_b:
            continue
        fid = getattr(shmr_params, name, 0.0)
        sig_a = errors_a[name]
        sig_b = errors_b[name]
        ratio = sig_b / sig_a if sig_a > 0 else np.inf
        print(f"  {name:<20s} {fid:>10.4f} "
              f"{sig_a:>12.6f} {sig_b:>12.6f} {ratio:>10.3f}")


def main():
    survey = make_dwarf_survey()
    shmr_params = SHMRParams(use_mass_dependent_scatter=True)

    # Base forecast config: dwarf regime, no z-evolution, fix high-mass params
    base_kwargs = dict(
        n_R_bins=10,
        n_Mh_bins=200,
        systematic_floor_fraction=0.05,
        include_nuisance_params=True,
        vary_z_evolution=False,
        fixed_params=["gamma_0", "log_M1_0"],
    )

    # ---------------------------------------------------------------
    # Test 1: sigma_obs = 0 vs. sigma_obs = 0.1
    # ---------------------------------------------------------------
    print("\n" + "#" * 75)
    print("# TEST 1: Impact of observational stellar mass scatter")
    print("#" * 75)

    fc_no_obs = ForecastConfig(**base_kwargs, sigma_log_Mstar_obs=0.0)
    fc_obs = ForecastConfig(**base_kwargs, sigma_log_Mstar_obs=0.1)

    errors_no_obs, _ = run_single(survey, shmr_params, fc_no_obs,
                                  "sigma_obs=0.0")
    errors_obs, _ = run_single(survey, shmr_params, fc_obs,
                               "sigma_obs=0.1")

    print_comparison(
        "sigma_obs=0.0", errors_no_obs,
        "sigma_obs=0.1", errors_obs,
        shmr_params,
    )
    print("\n  Ratio > 1 means sigma_obs DEGRADES constraints (expected).")

    # ---------------------------------------------------------------
    # Test 2: SMF covariance (Poisson+CV) vs. Poisson-only
    # ---------------------------------------------------------------
    print("\n" + "#" * 75)
    print("# TEST 2: SMF covariance (Poisson + cosmic variance) vs. Poisson-only")
    print("#" * 75)

    fc_smf = ForecastConfig(**base_kwargs, include_smf=True)
    fc_poisson = ForecastConfig(**base_kwargs, include_smf=False)

    errors_smf, _ = run_single(survey, shmr_params, fc_smf,
                               "include_smf=True")
    errors_poisson, _ = run_single(survey, shmr_params, fc_poisson,
                                   "include_smf=False")

    print_comparison(
        "SMF (Poisson+CV)", errors_smf,
        "Poisson-only", errors_poisson,
        shmr_params,
    )
    print("\n  Ratio < 1 means Poisson-only is more optimistic (smaller errors)")
    print("  because it ignores cosmic variance.")

    # ---------------------------------------------------------------
    # Test 3: Combined effect
    # ---------------------------------------------------------------
    print("\n" + "#" * 75)
    print("# TEST 3: Combined: sigma_obs=0.1 + SMF vs. baseline (sigma_obs=0, SMF)")
    print("#" * 75)

    fc_combined = ForecastConfig(
        **base_kwargs,
        sigma_log_Mstar_obs=0.1,
        include_smf=True,
    )
    errors_combined, _ = run_single(survey, shmr_params, fc_combined,
                                    "combined")

    print_comparison(
        "Baseline (s_obs=0, SMF)", errors_smf,
        "Combined (s_obs=0.1, SMF)", errors_combined,
        shmr_params,
    )

    print("\n" + "=" * 75)
    print("Phase 5 validation complete.")
    print("=" * 75)


if __name__ == "__main__":
    main()
