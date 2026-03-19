"""Tests for optional NFW reference cross-check utilities."""

import numpy as np

from sidm_stagev_forecast.reference_validation import (
    detect_reference_backends,
    run_nfw_reference_crosscheck,
)


def test_detect_reference_backends_has_expected_keys() -> None:
    availability = detect_reference_backends()
    assert set(availability.keys()) == {"colossus", "pyccl"}
    assert isinstance(availability["colossus"], bool)
    assert isinstance(availability["pyccl"], bool)


def test_run_nfw_reference_crosscheck_returns_required_columns() -> None:
    table, status = run_nfw_reference_crosscheck()

    required_columns = {
        "r_projected_kpc",
        "delta_sigma_local_msun_kpc2",
        "delta_sigma_analytic_msun_kpc2",
        "frac_diff_local_vs_analytic",
    }
    assert required_columns.issubset(set(table.columns))
    assert len(table) > 0
    assert np.all(np.isfinite(table["frac_diff_local_vs_analytic"]))
    assert len(status) == 2
