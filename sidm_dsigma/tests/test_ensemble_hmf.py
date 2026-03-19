"""Validation tests for HMF-based ensemble generation."""

import numpy as np

from sidm_stagev_forecast.ensemble import generate_ensemble


def _expected_log10_mass_mean_power_law(
    alpha: float,
    mass_min_msun: float,
    mass_max_msun: float,
) -> float:
    mass_grid = np.geomspace(mass_min_msun, mass_max_msun, 20_000)
    integrand = mass_grid ** (-alpha) * mass_grid
    cumulative = np.zeros_like(integrand)
    log_mass_grid = np.log(mass_grid)
    cumulative[1:] = np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(log_mass_grid))
    normalization = cumulative[-1]

    weighted_log_mass = np.zeros_like(integrand)
    weighted_log_mass[1:] = np.cumsum(
        0.5
        * (
            np.log10(mass_grid[1:]) * integrand[1:]
            + np.log10(mass_grid[:-1]) * integrand[:-1]
        )
        * np.diff(log_mass_grid)
    )
    return float(weighted_log_mass[-1] / normalization)


def test_hmf_mass_distribution_matches_power_law_expectation() -> None:
    """Sampled HMF masses should reproduce expected mean log-mass."""
    alpha = 1.9
    mass_min_msun = 1.0e14
    mass_max_msun = 1.0e15
    expected_mean_log10_mass = _expected_log10_mass_mean_power_law(
        alpha,
        mass_min_msun,
        mass_max_msun,
    )

    halos = generate_ensemble(
        mode="HMF",
        config_dict={
            "n_halos": 30_000,
            "seed": 13,
            "redshift": 0.4,
            "mass_min_msun": mass_min_msun,
            "mass_max_msun": mass_max_msun,
            "hmf_model": {"type": "power_law", "alpha": alpha},
            "selection_model": {"type": "none"},
            "concentration_scatter_dex": 0.0,
            "weight_mode": "equal",
        },
    )
    sampled_mean_log10_mass = np.mean(np.log10([halo["m200_msun"] for halo in halos]))
    assert np.isclose(sampled_mean_log10_mass, expected_mean_log10_mass, atol=0.015)


def test_hmf_convergence_with_halo_count() -> None:
    """Larger HMF samples should approach expected mean-mass more tightly."""
    alpha = 2.0
    mass_min_msun = 1.0e14
    mass_max_msun = 1.0e15
    expected_mean_log10_mass = _expected_log10_mass_mean_power_law(
        alpha,
        mass_min_msun,
        mass_max_msun,
    )

    small_errors: list[float] = []
    large_errors: list[float] = []
    for seed in range(8):
        halos_small = generate_ensemble(
            mode="HMF",
            config_dict={
                "n_halos": 200,
                "seed": seed,
                "redshift": 0.4,
                "mass_min_msun": mass_min_msun,
                "mass_max_msun": mass_max_msun,
                "hmf_model": {"type": "power_law", "alpha": alpha},
                "selection_model": {"type": "none"},
                "concentration_scatter_dex": 0.1,
                "weight_mode": "equal",
            },
        )
        halos_large = generate_ensemble(
            mode="HMF",
            config_dict={
                "n_halos": 5000,
                "seed": seed,
                "redshift": 0.4,
                "mass_min_msun": mass_min_msun,
                "mass_max_msun": mass_max_msun,
                "hmf_model": {"type": "power_law", "alpha": alpha},
                "selection_model": {"type": "none"},
                "concentration_scatter_dex": 0.1,
                "weight_mode": "equal",
            },
        )
        small_errors.append(
            abs(
                np.mean(np.log10([halo["m200_msun"] for halo in halos_small]))
                - expected_mean_log10_mass
            )
        )
        large_errors.append(
            abs(
                np.mean(np.log10([halo["m200_msun"] for halo in halos_large]))
                - expected_mean_log10_mass
            )
        )

    assert np.mean(large_errors) < np.mean(small_errors)
    assert np.mean(large_errors) < 0.03


def test_hmf_single_halo_limit_has_unit_weight() -> None:
    """Single-halo HMF ensemble should return one halo with weight 1."""
    halos = generate_ensemble(
        mode="HMF",
        config_dict={
            "n_halos": 1,
            "seed": 11,
            "redshift": 0.4,
            "mass_min_msun": 3.0e14,
            "mass_max_msun": 3.1e14,
            "hmf_model": {"type": "power_law", "alpha": 1.9},
            "selection_model": {"type": "threshold", "log10_m_cut": 14.0},
            "concentration_scatter_dex": 0.0,
            "weight_mode": "equal",
        },
    )

    assert len(halos) == 1
    assert np.isclose(halos[0]["weight"], 1.0)
    assert halos[0]["m200_msun"] > 0.0
