"""Compare CDM baselines: parametricSIDM tiny-sigma vs colossus DK14.

This script generates a convergence diagnostic figure for dwarf and cluster regimes.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sidm_stagev_forecast.config import DEFAULT_COSMOLOGY
from sidm_stagev_forecast.cosmology import rdelta
from sidm_stagev_forecast.io import append_figure_caption_entries
from sidm_stagev_forecast.profiles import sidm_profile_from_parametric_model


def _configure_colossus(root: Path) -> None:
    from colossus import settings
    from colossus.cosmology import cosmology

    cache_base = root / ".colossus_runtime"
    cache_base.mkdir(parents=True, exist_ok=True)
    settings.BASE_DIR = str(cache_base)
    settings.PERSISTENCE = "rw"

    cosmology_name = "stagev_planck18"
    if cosmology_name not in cosmology.cosmologies:
        cosmology.addCosmology(
            cosmology_name,
            {
                "flat": True,
                "H0": DEFAULT_COSMOLOGY.h * 100.0,
                "Om0": DEFAULT_COSMOLOGY.omega_m,
                "Ob0": 0.049,
                "sigma8": 0.811,
                "ns": 0.965,
            },
        )
    cosmology.setCosmology(cosmology_name)


def _dk14_cdm_profiles(
    r_kpc: np.ndarray,
    m200_msun: float,
    c200: float,
    z: float,
) -> tuple[np.ndarray, np.ndarray]:
    from colossus.halo.profile_dk14 import DK14Profile

    h = DEFAULT_COSMOLOGY.h
    profile = DK14Profile(M=m200_msun * h, c=c200, z=z, mdef="200c")

    r_kpc_h = r_kpc * h
    rho_h2 = profile.density(r_kpc_h)
    delta_sigma_h = profile.deltaSigma(r_kpc_h)

    rho_msun_kpc3 = np.asarray(rho_h2, dtype=float) * (h**2)
    delta_sigma_msun_kpc2 = np.asarray(delta_sigma_h, dtype=float) * h
    return rho_msun_kpc3, delta_sigma_msun_kpc2


def _parametric_tiny_sigma_profiles(
    r_kpc: np.ndarray,
    m200_msun: float,
    c200: float,
    z: float,
    sigma_tiny_cm2_g: float,
) -> tuple[np.ndarray, np.ndarray]:
    from sidm_stagev_forecast.projection import delta_sigma_of_R

    profile = sidm_profile_from_parametric_model(
        r_kpc=r_kpc,
        m200_msun=m200_msun,
        c200=c200,
        z=z,
        sigma_over_m=sigma_tiny_cm2_g,
    )
    ds = delta_sigma_of_R(
        r_projected_kpc=r_kpc,
        r_kpc=r_kpc,
        rho_msun_kpc3=profile["rho_msun_kpc3"],
        n_z=700,
    )["delta_sigma_msun_kpc2"]
    return np.asarray(profile["rho_msun_kpc3"], dtype=float), ds


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    _configure_colossus(root)

    sigma_tiny = 1.0e-4
    regimes = [
        {"label": "dwarf", "m200_msun": 1.0e10, "c200": 15.0, "z": 0.2, "r_min": 0.2, "r_max": 300.0},
        {"label": "cluster", "m200_msun": 1.0e14, "c200": 4.0, "z": 0.3, "r_min": 5.0, "r_max": 3000.0},
    ]

    fig, axes = plt.subplots(2, 4, figsize=(18, 8), constrained_layout=True)

    rows = []
    for col, regime in enumerate(regimes):
        r_kpc = np.geomspace(regime["r_min"], regime["r_max"], 160)

        rho_param, ds_param = _parametric_tiny_sigma_profiles(
            r_kpc=r_kpc,
            m200_msun=regime["m200_msun"],
            c200=regime["c200"],
            z=regime["z"],
            sigma_tiny_cm2_g=sigma_tiny,
        )
        rho_dk14, ds_dk14 = _dk14_cdm_profiles(
            r_kpc=r_kpc,
            m200_msun=regime["m200_msun"],
            c200=regime["c200"],
            z=regime["z"],
        )

        rho_ratio = rho_param / rho_dk14
        ds_ratio = ds_param / ds_dk14
        r200_kpc = rdelta(regime["m200_msun"], regime["z"], DEFAULT_COSMOLOGY, definition="200c")
        convergence_mask = (
            (r_kpc >= 0.05 * r200_kpc)
            & (r_kpc <= r200_kpc)
            & np.isfinite(rho_ratio)
            & np.isfinite(ds_ratio)
            & (rho_dk14 > 0.0)
            & (ds_dk14 > 0.0)
        )

        axes[0, 2 * col].loglog(r_kpc, rho_param, lw=2.0, label=r"parametricSIDM, $\sigma/m=10^{-4}$")
        axes[0, 2 * col].loglog(r_kpc, rho_dk14, lw=2.0, ls="--", label="colossus DK14")
        axes[0, 2 * col].set_xlabel("r [kpc]")
        axes[0, 2 * col].set_ylabel(r"$\rho$ [Msun / kpc$^3$]")
        axes[0, 2 * col].set_title(f"{regime['label']} density")
        axes[0, 2 * col].legend(fontsize=7)

        axes[0, 2 * col + 1].axhline(1.0, color="black", ls="--", lw=1.0)
        axes[0, 2 * col + 1].semilogx(r_kpc[convergence_mask], rho_ratio[convergence_mask], lw=2.0)
        axes[0, 2 * col + 1].set_xlabel("r [kpc]")
        axes[0, 2 * col + 1].set_ylabel(r"$\rho_{\rm param}/\rho_{\rm DK14}$")
        axes[0, 2 * col + 1].set_title(f"{regime['label']} density ratio (0.05-1 R200c)")
        axes[0, 2 * col + 1].set_ylim(0.5, 1.5)

        axes[1, 2 * col].loglog(r_kpc, ds_param, lw=2.0, label=r"parametricSIDM, $\sigma/m=10^{-4}$")
        axes[1, 2 * col].loglog(r_kpc, ds_dk14, lw=2.0, ls="--", label="colossus DK14")
        axes[1, 2 * col].set_xlabel("R [kpc]")
        axes[1, 2 * col].set_ylabel(r"$\Delta\Sigma$ [Msun / kpc$^2$]")
        axes[1, 2 * col].set_title(f"{regime['label']} DeltaSigma")

        axes[1, 2 * col + 1].axhline(1.0, color="black", ls="--", lw=1.0)
        axes[1, 2 * col + 1].semilogx(r_kpc[convergence_mask], ds_ratio[convergence_mask], lw=2.0)
        axes[1, 2 * col + 1].set_xlabel("R [kpc]")
        axes[1, 2 * col + 1].set_ylabel(r"$\Delta\Sigma_{\rm param}/\Delta\Sigma_{\rm DK14}$")
        axes[1, 2 * col + 1].set_title(f"{regime['label']} DeltaSigma ratio (0.05-1 R200c)")
        axes[1, 2 * col + 1].set_ylim(0.5, 1.5)

        rows.append(
            {
                "regime": regime["label"],
                "z": regime["z"],
                "sigma_tiny_cm2_g": sigma_tiny,
                "r200c_kpc": float(r200_kpc),
                "convergence_r_min_kpc": float(0.05 * r200_kpc),
                "convergence_r_max_kpc": float(r200_kpc),
                "rho_ratio_median": float(np.median(rho_ratio[convergence_mask])),
                "rho_ratio_max_abs_dev": float(np.max(np.abs(rho_ratio[convergence_mask] - 1.0))),
                "ds_ratio_median": float(np.median(ds_ratio[convergence_mask])),
                "ds_ratio_max_abs_dev": float(np.max(np.abs(ds_ratio[convergence_mask] - 1.0))),
            }
        )

    figure_path = root / "outputs" / "figures" / "cdm_engine_convergence_parametric_vs_colossus_dk14.png"
    fig.savefig(figure_path, dpi=220)
    plt.close(fig)

    table_path = root / "outputs" / "tables" / "cdm_engine_convergence_parametric_vs_colossus_dk14.csv"
    pd.DataFrame(rows).to_csv(table_path, index=False)

    created_at = pd.Timestamp(figure_path.stat().st_mtime, unit="s").tz_localize("UTC").tz_convert("America/Phoenix").isoformat()
    append_figure_caption_entries(
        caption_path=root / "outputs" / "figures" / "CAPTION.md",
        entries=[
            {
                "filename": figure_path.name,
                "created_at": created_at,
                "caption": (
                    "Convergence diagnostic comparing CDM profiles from parametricSIDM at tiny cross section "
                    "(sigma/m=1e-4 cm^2/g) versus colossus DK14 CDM. Panels show density and DeltaSigma curves "
                    "and their ratios for dwarf and cluster benchmarks. Ratio panels use the convergence window "
                    "0.05*R200c to 1.0*R200c."
                ),
            }
        ],
    )


if __name__ == "__main__":
    main()
