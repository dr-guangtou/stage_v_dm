"""Toy distinguishability metrics for SIDM vs CDM forecasts."""

from __future__ import annotations

import numpy as np


def fractional_error_model(
    r_projected_kpc: np.ndarray,
    regime: str,
    scenario: str = "baseline",
) -> np.ndarray:
    """Return toy per-bin fractional uncertainties for DeltaSigma."""
    radii = np.asarray(r_projected_kpc)
    regime = regime.lower()
    scenario = scenario.lower()

    if regime not in {"dwarf", "cluster"}:
        raise ValueError("regime must be 'dwarf' or 'cluster'.")
    if scenario not in {"baseline", "conservative"}:
        raise ValueError("scenario must be 'baseline' or 'conservative'.")

    pivot = np.median(radii)
    radial_scale = np.abs(np.log10(radii / pivot))

    if regime == "dwarf":
        if scenario == "baseline":
            return np.clip(0.12 + 0.04 * radial_scale, 0.10, 0.20)
        return np.clip(0.20 + 0.06 * radial_scale, 0.18, 0.35)

    if scenario == "baseline":
        return np.clip(0.05 + 0.02 * radial_scale, 0.03, 0.10)
    return np.clip(0.10 + 0.03 * radial_scale, 0.08, 0.18)


def delta_chi2(model: np.ndarray, reference: np.ndarray, sigma: np.ndarray) -> float:
    """Return Delta chi^2 = sum((model-reference)^2 / sigma^2)."""
    model = np.asarray(model)
    reference = np.asarray(reference)
    sigma = np.asarray(sigma)

    if np.any(sigma <= 0.0):
        raise ValueError("All uncertainty values must be positive.")

    residual = model - reference
    return float(np.sum((residual / sigma) ** 2))


def delta_chi2_with_fractional_error(
    model: np.ndarray,
    reference: np.ndarray,
    fractional_error: float | np.ndarray,
    floor: float = 1.0e-12,
) -> float:
    """Return Delta chi^2 using sigma_i = frac_i * |reference_i|."""
    reference = np.asarray(reference)
    sigma_abs = np.maximum(np.asarray(fractional_error) * np.abs(reference), floor)
    return delta_chi2(model=np.asarray(model), reference=reference, sigma=sigma_abs)


def sigma_separation_from_delta_chi2(delta_chi2_value: float) -> float:
    """Return effective Gaussian sigma-equivalent separation."""
    if delta_chi2_value < 0.0:
        raise ValueError("delta_chi2_value must be non-negative.")
    return float(np.sqrt(delta_chi2_value))


def evaluate_stacked_distinguishability(
    model: np.ndarray,
    reference: np.ndarray,
    fractional_error_scenarios: dict[str, float],
) -> dict[str, float]:
    """Evaluate Delta chi^2 and sigma-separation over named error scenarios."""
    results: dict[str, float] = {}
    for scenario_name, scenario_fraction in fractional_error_scenarios.items():
        metric_delta_chi2 = delta_chi2_with_fractional_error(
            model=model,
            reference=reference,
            fractional_error=scenario_fraction,
        )
        results[f"delta_chi2_{scenario_name}"] = metric_delta_chi2
        results[f"sigma_separation_{scenario_name}"] = sigma_separation_from_delta_chi2(
            metric_delta_chi2
        )
    return results


def required_uniform_fractional_precision(
    model: np.ndarray,
    reference: np.ndarray,
    target_sigma_significance: float,
    floor: float = 1.0e-12,
) -> float:
    """Estimate required constant fractional precision for target significance.

    Assumes sigma_i = frac * |reference_i| and significance^2 = Delta chi^2.
    """
    model = np.asarray(model)
    reference = np.asarray(reference)

    denominator = np.maximum(np.abs(reference), floor)
    relative_delta = (model - reference) / denominator

    target_delta_chi2 = target_sigma_significance**2
    return float(np.sqrt(np.sum(relative_delta**2) / target_delta_chi2))
