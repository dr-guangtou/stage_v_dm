"""Tests for YAML-driven ensemble configuration parsing and validation."""

from pathlib import Path

import numpy as np
import yaml

from sidm_stagev_forecast.config import log_radius_grid
from sidm_stagev_forecast.ensemble import generate_ensemble
from sidm_stagev_forecast.ensemble_yaml import load_ensemble_yaml_config


def _docs_path(filename: str) -> Path:
    return Path(__file__).resolve().parents[1] / "docs" / filename


def test_cluster_yaml_parses_and_generates_expected_median_mass() -> None:
    """Cluster YAML should route to HMF and produce high-mass ensemble."""
    parsed = load_ensemble_yaml_config(_docs_path("cluster_ensemble_config.yaml"))
    assert parsed.mode == "HMF"
    assert parsed.tier2_config["enabled"] is True
    assert parsed.tier2_config["outer_profile_model"] == "dk14_like"

    halos = generate_ensemble(parsed.mode, parsed.ensemble_config)
    median_mass = float(np.median([halo["m200_msun"] for halo in halos]))

    assert 1.5e14 <= median_mass <= 6.0e14
    assert np.isclose(np.sum([halo["weight"] for halo in halos]), 1.0)


def test_dwarf_yaml_parses_and_generates_expected_median_mass() -> None:
    """Dwarf YAML should route to SHMR and produce dwarf-scale halo masses."""
    parsed = load_ensemble_yaml_config(_docs_path("dwarf_ensemble_config.yaml"))
    assert parsed.mode == "SHMR"
    assert parsed.tier2_config["enabled"] is True
    assert parsed.tier2_config["r_match_mode"] in {"fraction_r200c", "fraction_r200m", "fixed_kpc"}

    halos = generate_ensemble(parsed.mode, parsed.ensemble_config)
    median_mass = float(np.median([halo["m200_msun"] for halo in halos]))

    assert 5.0e9 <= median_mass <= 2.0e10
    assert np.isclose(np.sum([halo["weight"] for halo in halos]), 1.0)


def test_yaml_projection_bins_match_config() -> None:
    """Projection bins should map exactly from YAML into runtime grids."""
    cluster = load_ensemble_yaml_config(_docs_path("cluster_ensemble_config.yaml"))
    dwarf = load_ensemble_yaml_config(_docs_path("dwarf_ensemble_config.yaml"))

    r_cluster = log_radius_grid(
        cluster.projection_config["r_min_kpc"],
        cluster.projection_config["r_max_kpc"],
        cluster.projection_config["n_r_bins"],
    )
    r_dwarf = log_radius_grid(
        dwarf.projection_config["r_min_kpc"],
        dwarf.projection_config["r_max_kpc"],
        dwarf.projection_config["n_r_bins"],
    )

    assert len(r_cluster) == cluster.projection_config["n_r_bins"]
    assert len(r_dwarf) == dwarf.projection_config["n_r_bins"]


def test_yaml_reproducibility_and_single_halo_limit() -> None:
    """Fixed seed should reproduce; N=1 should preserve single-halo limit."""
    parsed = load_ensemble_yaml_config(_docs_path("cluster_ensemble_config.yaml"))

    ensemble_first = generate_ensemble(parsed.mode, parsed.ensemble_config)
    ensemble_second = generate_ensemble(parsed.mode, parsed.ensemble_config)
    assert ensemble_first == ensemble_second

    single_halo_config = dict(parsed.ensemble_config)
    single_halo_config["n_halos"] = 1
    single_halo = generate_ensemble(parsed.mode, single_halo_config)
    assert len(single_halo) == 1
    assert np.isclose(single_halo[0]["weight"], 1.0)


def test_yaml_tier2_defaults_to_disabled_when_block_absent(tmp_path: Path) -> None:
    """Tier-1 YAML files without tier2 block should preserve original behavior."""
    configuration = {
        "mode": "HMF",
        "ensemble": {"N_halos": 8, "random_seed": 9},
        "redshift": {"model": "fixed", "value": 0.4},
        "mass_function": {"model": "power_law", "M_min": 1.0e14, "M_max": 1.0e15},
        "selection": {"type": "none"},
        "concentration": {"model": "maccio", "scatter_logc": 0.1},
        "stacking": {"weight_scheme": "equal"},
        "projection": {"R_min_kpc": 50, "R_max_kpc": 2000, "N_R_bins": 20},
        "sidm": {"sigma_over_m_grid": [0.2, 0.5]},
    }
    path = tmp_path / "tier1_only.yaml"
    with path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(configuration, handle)

    parsed = load_ensemble_yaml_config(path)
    assert parsed.tier2_config["enabled"] is False
