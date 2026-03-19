"""
Parameter sweep driver.

Varies a single survey parameter while holding all others fixed,
computing Fisher constraints at each point. Produces a scaling plot
and saves results as JSON.

Usage:
    uv run python scripts/run_sweep.py --base stage5_wide --param n_gal_total --values 1e6,5e6,1e7,5e7,1e8
    uv run python scripts/run_sweep.py --base stage5_wide --param log_Mstar_min --values 9.5,9.0,8.5,8.0
    uv run python scripts/run_sweep.py --base stage4_low_z --param area_deg2 --values 1000,5000,10000,14000
    uv run python scripts/run_sweep.py --base stage4_low_z --param dlog_Mstar --values 1.0,0.5,0.25

Outputs:
    outputs/sweep_{param}_{base}.png
    outputs/sweep_{param}_{base}.json
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shmr_fisher.config import (
    ForecastConfig,
    LensingConfig,
    NuisanceConfig,
    SHMRParams,
)
from shmr_fisher.fisher import (
    compute_fisher_matrix,
    extract_shmr_constraints,
    marginalized_errors,
)
from shmr_fisher.survey_configs import make_custom_survey, surveys as all_surveys


def parameter_sweep(
    base_survey_name: str,
    sweep_param: str,
    sweep_values: list[float],
    forecast_config: ForecastConfig | None = None,
    lensing_config: LensingConfig | None = None,
    nuisance_config: NuisanceConfig | None = None,
    shmr_params: SHMRParams | None = None,
) -> dict:
    """
    Sweep a single survey parameter and compute Fisher constraints at each value.

    Parameters
    ----------
    base_survey_name : str
        Key in survey_configs.surveys for the baseline survey.
    sweep_param : str
        SurveyConfig field to vary (e.g., 'n_gal_total', 'area_deg2').
    sweep_values : list of float
        Values to evaluate.
    forecast_config, lensing_config, nuisance_config, shmr_params :
        Optional overrides. Defaults used if None.

    Returns
    -------
    results : dict
        Contains sweep_values, param_names, and errors arrays.
    """
    if shmr_params is None:
        shmr_params = SHMRParams()
    if forecast_config is None:
        forecast_config = ForecastConfig(n_R_bins=8, n_Mh_bins=150)
    if lensing_config is None:
        lensing_config = LensingConfig()

    base = all_surveys[base_survey_name]

    shmr_names = None
    all_errors = []

    # When sweeping z_max, the number of varied parameters can change
    # (5 vs 9). Fix evolution params off for consistent output dimension.
    fc_local = forecast_config
    if sweep_param == "z_max":
        from dataclasses import replace
        fc_local = replace(forecast_config, vary_z_evolution=False)

    nuis = nuisance_config if fc_local.include_nuisance_params else None

    for val in sweep_values:
        # Build modified survey config
        kwargs = {
            "name": f"{base.name} ({sweep_param}={val})",
            "area_deg2": base.area_deg2,
            "z_min": base.z_min,
            "z_max": base.z_max,
            "n_gal_total": base.n_gal_total,
            "log_Mstar_min": base.log_Mstar_min,
            "dlog_Mstar": base.dlog_Mstar,
        }
        kwargs[sweep_param] = val

        # Preserve z-dependent completeness if not sweeping log_Mstar_min
        if sweep_param != "log_Mstar_min" and base.log_Mstar_min_func is not None:
            kwargs["log_Mstar_min_func"] = base.log_Mstar_min_func

        survey = make_custom_survey(**kwargs)

        print(f"  {sweep_param}={val}...", end="", flush=True)

        fisher, param_names, meta = compute_fisher_matrix(
            shmr_params, survey, lensing_config, fc_local, nuis,
        )

        eigvals = np.linalg.eigvalsh(fisher)
        if np.all(eigvals > 0):
            errs, names = extract_shmr_constraints(fisher, param_names)
            if shmr_names is None:
                shmr_names = names
            all_errors.append(errs)
            print(f" done (cond={meta.get('condition_number', 0):.1e})")
        else:
            all_errors.append(np.full(len(shmr_names or param_names), np.inf))
            print(" SINGULAR")

    return {
        "base_survey": base_survey_name,
        "sweep_param": sweep_param,
        "sweep_values": sweep_values,
        "shmr_param_names": shmr_names,
        "errors": np.array(all_errors),  # shape (n_values, n_params)
        "fiducial_params": shmr_params.to_dict(),
    }


def plot_sweep(results: dict, output_dir: str = "outputs") -> Path:
    """
    Generate and save a scaling plot from sweep results.

    Parameters
    ----------
    results : dict
        Output of parameter_sweep().
    output_dir : str

    Returns
    -------
    save_path : Path
    """
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    sweep_param = results["sweep_param"]
    sweep_values = results["sweep_values"]
    shmr_names = results["shmr_param_names"]
    errors = results["errors"]
    fid = results["fiducial_params"]

    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
              "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22"]

    fig, ax = plt.subplots(figsize=(9, 6))

    for i, (pname, color) in enumerate(zip(shmr_names, colors)):
        fid_val = abs(fid.get(pname, 1.0))
        if fid_val > 1e-10:
            frac_errors = errors[:, i] / fid_val * 100
        else:
            frac_errors = errors[:, i] * 100
        ax.plot(sweep_values, frac_errors, "o-", color=color, lw=2, ms=5,
                label=pname)

    ax.set_xlabel(sweep_param, fontsize=13)
    ax.set_ylabel(r"$\sigma / |\theta_\mathrm{fid}|$ [%]", fontsize=13)
    ax.set_title(
        f"Constraint scaling: {sweep_param} "
        f"(base: {results['base_survey']})",
        fontsize=13,
    )
    ax.legend(fontsize=9, ncol=2)
    ax.grid(True, alpha=0.3, which="both")

    # Use log scale for params that span large ranges
    if max(sweep_values) / min(sweep_values) > 10:
        ax.set_xscale("log")
    ax.set_yscale("log")

    plt.tight_layout()
    save_path = outdir / f"sweep_{sweep_param}_{results['base_survey']}.png"
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {save_path}")
    plt.close(fig)

    # Save JSON
    json_path = save_path.with_suffix(".json")
    json_data = {
        "base_survey": results["base_survey"],
        "sweep_param": sweep_param,
        "sweep_values": [float(v) for v in sweep_values],
        "shmr_param_names": shmr_names,
        "fractional_errors_percent": {
            pname: [
                float(errors[j, i] / abs(fid.get(pname, 1.0)) * 100)
                if abs(fid.get(pname, 1.0)) > 1e-10
                else float(errors[j, i] * 100)
                for j in range(len(sweep_values))
            ]
            for i, pname in enumerate(shmr_names)
        },
    }
    with open(json_path, "w") as f:
        json.dump(json_data, f, indent=2)
    print(f"Saved: {json_path}")

    # Append caption
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M PDT")
    caption = f"""
---

### `{save_path.name}`
**Created:** {timestamp}

Constraint scaling plot showing how marginalized fractional errors on SHMR parameters change as {sweep_param} is varied from {min(sweep_values)} to {max(sweep_values)}, using the {results['base_survey']} survey as the baseline. Each line corresponds to one SHMR parameter. Log-log axes if the sweep range spans more than one decade.

**Purpose:** Validates AT-7 (monotonic improvement with increasing constraining power) and identifies which parameters benefit most from changes to {sweep_param}.
"""
    with open(outdir / "CAPTION.md", "a") as f:
        f.write(caption)

    return save_path


def main():
    parser = argparse.ArgumentParser(description="SHMR Fisher Parameter Sweep")
    parser.add_argument("--base", required=True, help="Base survey name")
    parser.add_argument("--param", required=True, help="Parameter to sweep")
    parser.add_argument(
        "--values", required=True,
        help="Comma-separated values to evaluate",
    )
    parser.add_argument("--floor", type=float, default=0.0)
    parser.add_argument("--nuisance", action="store_true")
    parser.add_argument("--output-dir", default="outputs")
    args = parser.parse_args()

    values = [float(v) for v in args.values.split(",")]

    fc = ForecastConfig(
        n_R_bins=8,
        n_Mh_bins=150,
        systematic_floor_fraction=args.floor,
        include_nuisance_params=args.nuisance,
    )
    nc = NuisanceConfig() if args.nuisance else None

    print(f"Sweep: {args.param} = {values}")
    print(f"Base survey: {args.base}")

    results = parameter_sweep(
        args.base, args.param, values,
        forecast_config=fc, nuisance_config=nc,
    )
    plot_sweep(results, output_dir=args.output_dir)


if __name__ == "__main__":
    main()
