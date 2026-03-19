"""Tests for Tier-2 profile stitching behavior."""

import numpy as np

from sidm_stagev_forecast.stitch import stitch_inner_outer_profile


def test_logistic_logrho_stitch_respects_inner_outer_limits() -> None:
    """Hybrid profile should follow inner profile at small r and outer at large r."""
    r_kpc = np.geomspace(1.0, 5000.0, 500)
    rho_inner = 1.0e7 * (r_kpc / 10.0) ** (-1.2)
    rho_outer = 1.0e5 * (r_kpc / 1000.0) ** (-2.4)

    stitched = stitch_inner_outer_profile(
        r_kpc=r_kpc,
        rho_inner_sidm_msun_kpc3=rho_inner,
        rho_outer_reference_msun_kpc3=rho_outer,
        method="logistic_logrho_blend",
        r_match_kpc=900.0,
        smooth_width_dex=0.1,
        continuity="density",
    )

    small_radius_mask = r_kpc < 40.0
    large_radius_mask = r_kpc > 3500.0

    assert np.allclose(stitched[small_radius_mask], rho_inner[small_radius_mask], rtol=0.05)
    assert np.allclose(stitched[large_radius_mask], rho_outer[large_radius_mask], rtol=0.05)
    assert np.all(np.isfinite(stitched))
    assert np.all(stitched > 0.0)
