"""Tests for Tier-1 halo ensemble utilities."""

import numpy as np

from sidm_stagev_forecast.config import (
    CLUSTER_HMF_ENSEMBLE_CONFIG_EXAMPLE,
    DWARF_SHMR_ENSEMBLE_CONFIG_EXAMPLE,
)
from sidm_stagev_forecast.ensemble import (
    generate_ensemble,
    sample_halo_ensemble,
    summarize_ensemble,
)


def test_sample_halo_ensemble_reproducible_for_fixed_seed() -> None:
    """Sampling should be deterministic for a fixed random seed."""
    first = sample_halo_ensemble(
        n_halos=64,
        mean_mass_msun=3.0e14,
        mass_scatter_dex=0.2,
        z=0.4,
        concentration_scatter_dex=0.1,
        seed=42,
    )
    second = sample_halo_ensemble(
        n_halos=64,
        mean_mass_msun=3.0e14,
        mass_scatter_dex=0.2,
        z=0.4,
        concentration_scatter_dex=0.1,
        seed=42,
    )

    assert first == second


def test_sample_halo_ensemble_equal_weights_and_positive_values() -> None:
    """Sampled ensemble should satisfy basic positivity and normalization."""
    halos = sample_halo_ensemble(
        n_halos=50,
        mean_mass_msun=3.0e14,
        mass_scatter_dex=0.2,
        z=0.4,
        concentration_scatter_dex=0.1,
        seed=7,
    )

    weights = np.asarray([halo["weight"] for halo in halos])
    masses = np.asarray([halo["m200_msun"] for halo in halos])
    concentrations = np.asarray([halo["c200"] for halo in halos])

    assert np.all(weights > 0.0)
    assert np.isclose(np.sum(weights), 1.0)
    assert np.all(masses > 0.0)
    assert np.all(concentrations > 0.0)


def test_summarize_ensemble_recovers_requested_scatter_scale() -> None:
    """Sample standard deviation in log10(M) should track configured scatter."""
    requested_scatter = 0.2
    halos = sample_halo_ensemble(
        n_halos=20_000,
        mean_mass_msun=3.0e14,
        mass_scatter_dex=requested_scatter,
        z=0.4,
        concentration_scatter_dex=0.1,
        seed=3,
    )
    summary = summarize_ensemble(halos)

    assert np.isclose(summary["std_log10_mass_dex"], requested_scatter, atol=0.01)
    assert np.isclose(summary["weight_sum"], 1.0, atol=1.0e-12)


def test_documented_config_examples_generate_valid_ensembles() -> None:
    """Config examples should run and produce positive masses and normalized weights."""
    cluster_halos = generate_ensemble("HMF", CLUSTER_HMF_ENSEMBLE_CONFIG_EXAMPLE)
    dwarf_halos = generate_ensemble("SHMR", DWARF_SHMR_ENSEMBLE_CONFIG_EXAMPLE)

    for halos in (cluster_halos, dwarf_halos):
        masses = np.asarray([halo["m200_msun"] for halo in halos])
        weights = np.asarray([halo["weight"] for halo in halos])
        concentrations = np.asarray([halo["c200"] for halo in halos])

        assert np.all(masses > 0.0)
        assert np.all(concentrations > 0.0)
        assert np.isclose(np.sum(weights), 1.0)
