"""Tests for YAML-driven ensemble configuration parsing and validation."""

from pathlib import Path

import numpy as np

from sidm_stagev_forecast.config import log_radius_grid
from sidm_stagev_forecast.ensemble import generate_ensemble
from sidm_stagev_forecast.ensemble_yaml import load_ensemble_yaml_config


def _docs_path(filename: str) -> Path:
    return Path(__file__).resolve().parents[1] / "docs" / filename


def test_cluster_yaml_parses_and_generates_expected_median_mass() -> None:
    """Cluster YAML should route to HMF and produce high-mass ensemble."""
    parsed = load_ensemble_yaml_config(_docs_path("cluster_ensemble_config.yaml"))
    assert parsed.mode == "HMF"

    halos = generate_ensemble(parsed.mode, parsed.ensemble_config)
    median_mass = float(np.median([halo["m200_msun"] for halo in halos]))

    assert 1.5e14 <= median_mass <= 6.0e14
    assert np.isclose(np.sum([halo["weight"] for halo in halos]), 1.0)


def test_dwarf_yaml_parses_and_generates_expected_median_mass() -> None:
    """Dwarf YAML should route to SHMR and produce dwarf-scale halo masses."""
    parsed = load_ensemble_yaml_config(_docs_path("dwarf_ensemble_config.yaml"))
    assert parsed.mode == "SHMR"

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
