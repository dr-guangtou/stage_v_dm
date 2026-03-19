"""Tests for ensemble-oriented forecast metrics."""

import numpy as np

from sidm_stagev_forecast.forecast import (
    delta_chi2_with_fractional_error,
    evaluate_stacked_distinguishability,
    sigma_separation_from_delta_chi2,
)


def test_delta_chi2_with_fractional_error_matches_manual_expression() -> None:
    reference = np.asarray([10.0, 20.0, 40.0])
    model = np.asarray([11.0, 18.0, 44.0])
    fractional_error = 0.1

    manual = np.sum(((model - reference) / (fractional_error * np.abs(reference))) ** 2)
    measured = delta_chi2_with_fractional_error(model, reference, fractional_error)

    assert np.isclose(measured, manual, rtol=0.0, atol=1.0e-14)


def test_sigma_separation_is_sqrt_delta_chi2() -> None:
    assert np.isclose(sigma_separation_from_delta_chi2(9.0), 3.0)


def test_evaluate_stacked_distinguishability_contains_expected_fields() -> None:
    reference = np.asarray([5.0, 8.0, 12.0, 20.0])
    model = np.asarray([5.5, 7.4, 11.0, 19.0])
    scenarios = {"optimistic": 0.05, "conservative": 0.10}

    results = evaluate_stacked_distinguishability(
        model=model,
        reference=reference,
        fractional_error_scenarios=scenarios,
    )

    assert "delta_chi2_optimistic" in results
    assert "delta_chi2_conservative" in results
    assert "sigma_separation_optimistic" in results
    assert "sigma_separation_conservative" in results
    assert results["delta_chi2_optimistic"] > results["delta_chi2_conservative"]
