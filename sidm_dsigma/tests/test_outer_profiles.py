"""Tests for DK14-like outer-profile utilities."""

import os
from pathlib import Path

import numpy as np
import pytest

from sidm_stagev_forecast.config import DEFAULT_COSMOLOGY
from sidm_stagev_forecast.outer_profiles import build_dk14_outer_profile, build_outer_profile
from sidm_stagev_forecast.plotting import log_slope


def test_dk14_like_profile_positive_and_finite() -> None:
    """DK14-like profile components should remain finite and positive."""
    r_kpc = np.geomspace(5.0, 8000.0, 240)
    profile = build_dk14_outer_profile(
        r_kpc=r_kpc,
        mass_msun=3.0e14,
        z=0.4,
        concentration=4.0,
        mass_def="200c",
        cosmology=DEFAULT_COSMOLOGY,
        dk14_params={"regime": "cluster"},
    )

    assert np.all(np.isfinite(profile["rho_total_msun_kpc3"]))
    assert np.all(np.isfinite(profile["rho_inner_baseline_msun_kpc3"]))
    assert np.all(np.isfinite(profile["rho_outer_term_msun_kpc3"]))
    assert np.all(profile["rho_total_msun_kpc3"] > 0.0)
    assert np.all(profile["rho_outer_term_msun_kpc3"] > 0.0)


def test_dk14_like_profile_has_outer_steepening_signature() -> None:
    """A DK14-like profile should show steepening relative to inner halo behavior."""
    r_kpc = np.geomspace(5.0, 8000.0, 320)
    profile = build_dk14_outer_profile(
        r_kpc=r_kpc,
        mass_msun=3.0e14,
        z=0.4,
        concentration=4.0,
        mass_def="200c",
        cosmology=DEFAULT_COSMOLOGY,
        dk14_params={"regime": "cluster"},
    )
    slope = log_slope(r_kpc, profile["rho_total_msun_kpc3"])

    assert np.min(slope) < -2.8


def test_outer_profile_selector_rejects_unknown_model() -> None:
    """Outer profile selector should fail fast for unsupported model names."""
    with pytest.raises(ValueError):
        build_outer_profile(
            r_kpc=np.geomspace(5.0, 8000.0, 40),
            mass_msun=3.0e14,
            z=0.4,
            concentration=4.0,
            mass_def="200c",
            cosmology=DEFAULT_COSMOLOGY,
            outer_profile_model="unknown_model",
            outer_params={"regime": "cluster"},
        )


def test_colossus_diemer23_outer_profile_optional(tmp_path: Path) -> None:
    """Diemer23 colossus backend should return a finite positive profile when available."""
    pytest.importorskip("colossus")
    os.environ["COLOSSUS_CACHE_DIR"] = str(tmp_path / "colossus_cache")
    profile = build_outer_profile(
        r_kpc=np.geomspace(5.0, 8000.0, 80),
        mass_msun=3.0e14,
        z=0.4,
        concentration=4.0,
        mass_def="200c",
        cosmology=DEFAULT_COSMOLOGY,
        outer_profile_model="colossus_diemer23",
        outer_params={"regime": "cluster"},
    )

    assert np.all(np.isfinite(profile["rho_total_msun_kpc3"]))
    assert np.all(profile["rho_total_msun_kpc3"] > 0.0)
    assert profile["metadata"]["model"] in {"colossus_diemer23", "dk14_like"}
