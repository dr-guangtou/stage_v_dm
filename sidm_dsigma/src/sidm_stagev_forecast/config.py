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


DEFAULT_COSMOLOGY = CosmologyConfig()
DEFAULT_FORECAST_CONFIG = ForecastConfig()

DWARF_BENCHMARK = BenchmarkConfig(
    label="dwarf",
    m200_msun=1.0e10,
    c200=15.0,
    redshift=0.3,
    sigma_over_m_grid_cm2_g=(0.0, 0.2, 0.5, 1.0, 2.0),
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
    sigma_over_m_grid_cm2_g=(0.0, 0.2, 0.5, 1.0, 2.0),
    r_3d_min_kpc=5.0,
    r_3d_max_kpc=5000.0,
    r_lensing_min_kpc=30.0,
    r_lensing_max_kpc=3000.0,
)

BENCHMARKS = (DWARF_BENCHMARK, CLUSTER_BENCHMARK)


def log_radius_grid(r_min_kpc: float, r_max_kpc: float, n_points: int) -> np.ndarray:
    """Return a logarithmically spaced radius grid in kpc."""
    return np.geomspace(r_min_kpc, r_max_kpc, n_points)
