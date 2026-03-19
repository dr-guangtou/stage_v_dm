"""Validation tests for SHMR-based ensemble generation."""

import numpy as np

from sidm_stagev_forecast.ensemble import generate_ensemble


def test_shmr_mass_distribution_matches_analytic_logspace_mapping() -> None:
    """With zero halo scatter, SHMR halo-mass moments should match analytic mapping."""
    mean_log10_mstar = 7.5
    sigma_log10_mstar = 0.5
    a = 10.5
    b = 1.4
    pivot = 7.5

    halos = generate_ensemble(
        mode="SHMR",
        config_dict={
            "n_halos": 25_000,
            "seed": 19,
            "redshift": 0.2,
            "stellar_mass_distribution": {
                "type": "lognormal",
                "mean_log10_mstar": mean_log10_mstar,
                "sigma_log10_mstar": sigma_log10_mstar,
            },
            "shmr_model": {"type": "power_law", "a": a, "b": b, "pivot_log10_mstar": pivot},
            "halo_scatter_dex": 0.0,
            "mhalo_min_msun": 1.0e8,
            "mhalo_max_msun": 1.0e14,
            "concentration_scatter_dex": 0.0,
            "weight_mode": "equal",
        },
    )

    expected_mean = a + b * (mean_log10_mstar - pivot)
    expected_std = b * sigma_log10_mstar

    sampled_log10_halo_mass = np.log10([halo["m200_msun"] for halo in halos])
    assert np.isclose(np.mean(sampled_log10_halo_mass), expected_mean, atol=0.02)
    assert np.isclose(np.std(sampled_log10_halo_mass), expected_std, atol=0.02)


def test_shmr_convergence_with_halo_count() -> None:
    """Larger SHMR samples should better recover expected log halo-mass mean."""
    mean_log10_mstar = 7.5
    a = 10.4
    b = 1.5
    pivot = 7.5
    expected_mean = a + b * (mean_log10_mstar - pivot)

    common_config = {
        "redshift": 0.2,
        "stellar_mass_distribution": {
            "type": "lognormal",
            "mean_log10_mstar": mean_log10_mstar,
            "sigma_log10_mstar": 0.5,
        },
        "shmr_model": {"type": "power_law", "a": a, "b": b, "pivot_log10_mstar": pivot},
        "halo_scatter_dex": 0.3,
        "mhalo_min_msun": 1.0e6,
        "mhalo_max_msun": 1.0e16,
        "concentration_scatter_dex": 0.1,
        "weight_mode": "equal",
    }
    small_errors: list[float] = []
    large_errors: list[float] = []
    for seed in range(8):
        halos_small = generate_ensemble(
            mode="SHMR",
            config_dict={"n_halos": 250, "seed": seed, **common_config},
        )
        halos_large = generate_ensemble(
            mode="SHMR",
            config_dict={"n_halos": 5000, "seed": seed, **common_config},
        )
        small_errors.append(
            abs(np.mean(np.log10([halo["m200_msun"] for halo in halos_small])) - expected_mean)
        )
        large_errors.append(
            abs(np.mean(np.log10([halo["m200_msun"] for halo in halos_large])) - expected_mean)
        )

    assert np.mean(large_errors) < np.mean(small_errors)
    assert np.mean(large_errors) < 0.03


def test_shmr_single_halo_limit_is_deterministic_without_scatter() -> None:
    """Single-halo SHMR output should be deterministic when scatter is disabled."""
    halos = generate_ensemble(
        mode="SHMR",
        config_dict={
            "n_halos": 1,
            "seed": 4,
            "redshift": 0.2,
            "stellar_mass_distribution": {
                "type": "lognormal",
                "mean_log10_mstar": 7.5,
                "sigma_log10_mstar": 0.0,
            },
            "shmr_model": {"type": "power_law", "a": 10.5, "b": 1.4, "pivot_log10_mstar": 7.5},
            "halo_scatter_dex": 0.0,
            "concentration_scatter_dex": 0.0,
            "weight_mode": "equal",
        },
    )
    expected_mass = 10.0**10.5

    assert len(halos) == 1
    assert np.isclose(halos[0]["m200_msun"], expected_mass)
    assert np.isclose(halos[0]["weight"], 1.0)
