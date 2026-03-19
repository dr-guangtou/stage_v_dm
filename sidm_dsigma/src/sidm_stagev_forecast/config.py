"""Configuration objects for benchmark SIDM/CDM forecasts."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class CosmologyConfig:
    """Flat LCDM configuration used for all computations."""

    omega_m: float = 0.3
    omega_lambda: float = 0.7
    h: float = 0.7


@dataclass(frozen=True)
class BenchmarkConfig:
    """Single benchmark halo configuration."""

    label: str
    m200_msun: float
    c200: float
    redshift: float
    sigma_over_m_grid_cm2_g: tuple[float, ...]
    r_3d_min_kpc: float
    r_3d_max_kpc: float
    r_lensing_min_kpc: float
    r_lensing_max_kpc: float


@dataclass(frozen=True)
class ForecastConfig:
    """Global forecast settings."""

    n_r_3d: int = 220
    n_r_lensing: int = 40
    default_mass_definition: str = "200c"


@dataclass(frozen=True)
class EnsembleBenchmarkConfig:
    """Tier-1 ensemble benchmark settings for stacked cluster-like forecasts."""

    n_halos: int = 100
    mean_mass_msun: float = 3.0e14
    mass_scatter_dex: float = 0.2
    concentration_scatter_dex: float = 0.1
    redshift: float = 0.4
    sigma_over_m_grid_cm2_g: tuple[float, ...] = (0.05, 0.1, 0.2, 0.3)
    n_r_3d: int = 220
    n_r_lensing_per_halo: int = 56
    r_common_min_kpc: float = 30.0
    r_common_max_kpc: float = 3000.0
    n_r_common: int = 44
    projection_n_z: int = 700
    seed: int = 7


@dataclass(frozen=True)
class Tier2Config:
    """Tier-2 hybrid (SIDM inner + DK14-like outskirts) controls."""

    enabled: bool = False
    outer_profile_model: str = "dk14_like"
    stitch_method: str = "logistic_logrho_blend"
    r_match_mode: str = "fraction_r200m"
    r_match_value: float = 0.8
    smooth_width_dex: float = 0.15
    continuity: str = "density"
    regime: str = "cluster"


@dataclass(frozen=True)
class Tier3Config:
    """Tier-3 empirical outer-correction controls."""

    enabled: bool = False
    correction_model: str = "rt_gamma_shift"
    sigma_pivot: float = 1.0
    apply_to_regimes: tuple[str, ...] = ("cluster",)
    calibration_mode: str = "manual_preset"
    preset: str = "none"
    regime: str = "cluster"


DEFAULT_COSMOLOGY = CosmologyConfig()
DEFAULT_FORECAST_CONFIG = ForecastConfig()
DEFAULT_ENSEMBLE_BENCHMARK = EnsembleBenchmarkConfig()
DEFAULT_TIER2_CONFIG = Tier2Config()
DEFAULT_TIER3_CONFIG = Tier3Config()

CLUSTER_HMF_ENSEMBLE_CONFIG_EXAMPLE: dict[str, object] = {
    "n_halos": 100,
    "seed": 7,
    "redshift": 0.4,
    "mass_min_msun": 1.0e14,
    "mass_max_msun": 1.0e15,
    "n_mass_grid": 512,
    "hmf_model": {"type": "power_law", "alpha": 1.9},
    "selection_model": {"type": "threshold", "log10_m_cut": 14.0},
    "concentration_model": {"type": "maccio"},
    "concentration_scatter_dex": 0.15,
    "weight_mode": "equal",
    "tier2": {
        "enabled": True,
        "outer_profile_model": "dk14_like",
        "stitch_method": "logistic_logrho_blend",
        "r_match_mode": "fraction_r200m",
        "r_match_value": 0.8,
        "smooth_width_dex": 0.15,
        "continuity": "density",
        "regime": "cluster",
    },
    "tier3": {
        "enabled": True,
        "correction_model": "rt_gamma_shift",
        "sigma_pivot": 1.0,
        "apply_to_regimes": ("cluster",),
        "calibration_mode": "manual_preset",
        "preset": "moderate",
        "regime": "cluster",
        "rt_shift": {"A_rt": -0.08, "use_concentration_dependence": True, "A_c": 0.5},
        "gamma_shift": {"A_gamma": -0.2},
    },
}

DWARF_SHMR_ENSEMBLE_CONFIG_EXAMPLE: dict[str, object] = {
    "n_halos": 200,
    "seed": 17,
    "redshift": 0.2,
    "stellar_mass_distribution": {
        "type": "lognormal",
        "mean_log10_mstar": 7.5,
        "sigma_log10_mstar": 0.5,
    },
    "shmr_model": {"type": "power_law", "a": 10.5, "b": 1.4, "pivot_log10_mstar": 7.5},
    "halo_scatter_dex": 0.3,
    "mhalo_min_msun": 1.0e9,
    "mhalo_max_msun": 3.0e11,
    "concentration_model": {"type": "maccio"},
    "concentration_scatter_dex": 0.15,
    "weight_mode": "equal",
    "tier2": {
        "enabled": False,
        "outer_profile_model": "dk14_like",
        "stitch_method": "logistic_logrho_blend",
        "r_match_mode": "fraction_r200c",
        "r_match_value": 1.0,
        "smooth_width_dex": 0.18,
        "continuity": "density",
        "regime": "dwarf",
    },
    "tier3": {
        "enabled": False,
        "correction_model": "rt_gamma_shift",
        "sigma_pivot": 1.0,
        "apply_to_regimes": ("cluster",),
        "calibration_mode": "manual_preset",
        "preset": "none",
        "regime": "dwarf",
    },
}

DWARF_BENCHMARK = BenchmarkConfig(
    label="dwarf",
    m200_msun=1.0e10,
    c200=15.0,
    redshift=0.3,
    sigma_over_m_grid_cm2_g=(0.0, 10.0, 20.0, 50.0, 100.0),
    r_3d_min_kpc=0.1,
    r_3d_max_kpc=300.0,
    r_lensing_min_kpc=3.0,
    r_lensing_max_kpc=300.0,
)

CLUSTER_BENCHMARK = BenchmarkConfig(
    label="cluster",
    m200_msun=1.0e14,
    c200=4.0,
    redshift=0.3,
    sigma_over_m_grid_cm2_g=(0.0, 0.05, 0.1, 0.2, 0.3),
    r_3d_min_kpc=5.0,
    r_3d_max_kpc=5000.0,
    r_lensing_min_kpc=30.0,
    r_lensing_max_kpc=3000.0,
)

BENCHMARKS = (DWARF_BENCHMARK, CLUSTER_BENCHMARK)


def log_radius_grid(r_min_kpc: float, r_max_kpc: float, n_points: int) -> np.ndarray:
    """Return a logarithmically spaced radius grid in kpc."""
    return np.geomspace(r_min_kpc, r_max_kpc, n_points)
