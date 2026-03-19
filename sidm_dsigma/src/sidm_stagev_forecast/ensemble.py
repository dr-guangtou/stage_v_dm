"""Flexible halo-ensemble generation for HMF and SHMR selection modes."""

from __future__ import annotations

import warnings
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True)
class EnsembleConfig:
    """Backward-compatible Tier-1 log-normal ensemble settings."""

    n_halos: int = 100
    mean_mass_msun: float = 3.0e14
    mass_scatter_dex: float = 0.2
    redshift: float = 0.4
    concentration_scatter_dex: float = 0.1
    seed: int = 7


def concentration_from_mass_maccio(
    m200_msun: np.ndarray | float,
    z: np.ndarray | float,
) -> np.ndarray:
    """Return a lightweight Dutton-Maccio-like concentration relation."""
    mass = np.asarray(m200_msun, dtype=float)
    redshift = np.asarray(z, dtype=float)

    if np.any(mass <= 0.0):
        raise ValueError("All halo masses must be positive.")
    if np.any(redshift < 0.0):
        raise ValueError("Redshift must be non-negative.")

    pivot_mass_msun = 2.0e12
    return 5.71 * (mass / pivot_mass_msun) ** (-0.084) * (1.0 + redshift) ** (-0.47)


def _validate_n_halos(config_dict: dict[str, Any]) -> int:
    n_halos = int(config_dict["n_halos"])
    if n_halos <= 0:
        raise ValueError("n_halos must be positive.")
    return n_halos


def _rng_from_config(config_dict: dict[str, Any]) -> np.random.Generator:
    return np.random.default_rng(config_dict.get("seed"))


def _sample_redshifts(
    n_halos: int,
    config_dict: dict[str, Any],
    random_generator: np.random.Generator,
) -> np.ndarray:
    if "redshift_distribution" not in config_dict:
        redshift = float(config_dict.get("redshift", 0.4))
        if redshift < 0.0:
            raise ValueError("redshift must be non-negative.")
        return np.full(n_halos, redshift)

    distribution = dict(config_dict["redshift_distribution"])
    distribution_type = str(distribution.get("type", "uniform")).lower()
    if distribution_type == "uniform":
        z_min = float(distribution["z_min"])
        z_max = float(distribution["z_max"])
        if z_min < 0.0 or z_max < z_min:
            raise ValueError("Invalid uniform redshift bounds.")
        return random_generator.uniform(z_min, z_max, size=n_halos)
    if distribution_type == "discrete":
        z_values = np.asarray(distribution["z_values"], dtype=float)
        probabilities = np.asarray(distribution["probabilities"], dtype=float)
        if len(z_values) != len(probabilities):
            raise ValueError("z_values and probabilities must have equal length.")
        if np.any(z_values < 0.0) or np.any(probabilities < 0.0):
            raise ValueError("Discrete redshift distribution values must be non-negative.")
        probabilities = probabilities / np.sum(probabilities)
        return random_generator.choice(z_values, size=n_halos, p=probabilities)
    if distribution_type == "gaussian":
        mean = float(distribution["mean"])
        sigma = float(distribution["sigma"])
        z_min = float(distribution.get("z_min", 0.0))
        z_max = float(distribution.get("z_max", mean + 5.0 * sigma))
        if sigma <= 0.0:
            raise ValueError("Gaussian redshift sigma must be positive.")
        if z_min < 0.0 or z_max <= z_min:
            raise ValueError("Invalid gaussian redshift truncation bounds.")

        accepted: list[float] = []
        while len(accepted) < n_halos:
            draws = random_generator.normal(mean, sigma, size=n_halos)
            valid_draws = draws[(draws >= z_min) & (draws <= z_max)]
            accepted.extend(valid_draws.tolist())
        return np.asarray(accepted[:n_halos], dtype=float)

    raise ValueError("Unsupported redshift_distribution type.")


def _mass_grid_from_config(config_dict: dict[str, Any]) -> np.ndarray:
    if "mass_grid_msun" in config_dict:
        mass_grid = np.asarray(config_dict["mass_grid_msun"], dtype=float)
    else:
        mass_min_msun = float(config_dict.get("mass_min_msun", 1.0e14))
        mass_max_msun = float(config_dict.get("mass_max_msun", 1.0e15))
        n_mass = int(config_dict.get("n_mass_grid", 512))
        mass_grid = np.geomspace(mass_min_msun, mass_max_msun, n_mass)

    if np.any(mass_grid <= 0.0):
        raise ValueError("mass_grid_msun values must be positive.")
    if not np.all(np.diff(mass_grid) > 0.0):
        raise ValueError("mass_grid_msun must be strictly increasing.")
    return mass_grid


def _resolve_hmf_callable(config_dict: dict[str, Any]) -> Callable[[np.ndarray, float], np.ndarray]:
    if "hmf_callable" in config_dict:
        return config_dict["hmf_callable"]

    hmf_model = dict(config_dict.get("hmf_model", {"type": "power_law", "alpha": 1.9}))
    model_type = str(hmf_model.get("type", "power_law")).lower()

    if model_type == "power_law":
        alpha = float(hmf_model.get("alpha", 1.9))
        normalization = float(hmf_model.get("normalization", 1.0))

        def power_law_hmf(mass_msun: np.ndarray, z: float) -> np.ndarray:
            del z
            return normalization * np.asarray(mass_msun, dtype=float) ** (-alpha)

        return power_law_hmf

    if model_type == "schechter":
        alpha = float(hmf_model.get("alpha", 1.9))
        m_cut_msun = float(hmf_model.get("m_cut_msun", 1.0e15))
        normalization = float(hmf_model.get("normalization", 1.0))

        def schechter_hmf(mass_msun: np.ndarray, z: float) -> np.ndarray:
            del z
            mass = np.asarray(mass_msun, dtype=float)
            return normalization * mass ** (-alpha) * np.exp(-mass / m_cut_msun)

        return schechter_hmf
    if model_type == "tinker08":
        # Prefer optional external implementations; fallback to a smooth proxy.
        try:
            from colossus.cosmology import cosmology as colossus_cosmology
            from colossus.lss import mass_function
        except Exception:
            warnings.warn(
                "Tinker08 requested but colossus not available; using power-law proxy alpha=1.9.",
                RuntimeWarning,
                stacklevel=2,
            )
            alpha = float(hmf_model.get("fallback_alpha", 1.9))
            normalization = float(hmf_model.get("normalization", 1.0))

            def fallback_tinker08(mass_msun: np.ndarray, z: float) -> np.ndarray:
                del z
                return normalization * np.asarray(mass_msun, dtype=float) ** (-alpha)

            return fallback_tinker08

        def tinker08_hmf(mass_msun: np.ndarray, z: float) -> np.ndarray:
            if "planck18" not in colossus_cosmology.cosmologies:
                colossus_cosmology.addCosmology(
                    "planck18",
                    {
                        "flat": True,
                        "H0": 67.4,
                        "Om0": 0.315,
                        "Ob0": 0.049,
                        "sigma8": 0.811,
                        "ns": 0.965,
                    },
                )
            colossus_cosmology.setCosmology("planck18")

            mass = np.asarray(mass_msun, dtype=float)
            dndlnm = mass_function.massFunction(
                x=mass,
                z=z,
                mdef="200c",
                model="tinker08",
                q_out="dndlnM",
            )
            dndm = dndlnm / np.maximum(mass, 1.0e-30)
            dndm[~np.isfinite(dndm)] = 0.0
            return np.maximum(dndm, 0.0)

        return tinker08_hmf

    raise ValueError("Unsupported hmf_model type.")


def _resolve_selection_callable(
    config_dict: dict[str, Any],
) -> Callable[[np.ndarray, np.ndarray], np.ndarray]:
    if "selection_callable" in config_dict:
        return config_dict["selection_callable"]

    selection_model = dict(config_dict.get("selection_model", {"type": "none"}))
    selection_type = str(selection_model.get("type", "none")).lower()

    if selection_type == "none":

        def no_selection(mass_msun: np.ndarray, z: np.ndarray) -> np.ndarray:
            del z
            return np.ones_like(np.asarray(mass_msun, dtype=float))

        return no_selection

    if selection_type == "threshold":
        log10_m_cut = float(selection_model.get("log10_m_cut", 14.0))

        def threshold_selection(mass_msun: np.ndarray, z: np.ndarray) -> np.ndarray:
            del z
            return (np.log10(np.asarray(mass_msun, dtype=float)) >= log10_m_cut).astype(float)

        return threshold_selection

    if selection_type == "logistic":
        log10_m_cut = float(selection_model.get("log10_m_cut", 14.0))
        sigma_log10_m = float(selection_model.get("sigma_log10_m", 0.2))
        if sigma_log10_m <= 0.0:
            raise ValueError("sigma_log10_m must be positive for logistic selection.")

        def logistic_selection(mass_msun: np.ndarray, z: np.ndarray) -> np.ndarray:
            del z
            log_mass = np.log10(np.asarray(mass_msun, dtype=float))
            return 1.0 / (1.0 + np.exp(-(log_mass - log10_m_cut) / sigma_log10_m))

        return logistic_selection

    raise ValueError("Unsupported selection_model type.")


def _sample_mass_from_hmf(
    n_halos: int,
    mass_grid_msun: np.ndarray,
    z_samples: np.ndarray,
    hmf_callable: Callable[[np.ndarray, float], np.ndarray],
    selection_callable: Callable[[np.ndarray, np.ndarray], np.ndarray],
    random_generator: np.random.Generator,
) -> np.ndarray:
    samples = np.empty(n_halos, dtype=float)
    log_mass_grid = np.log(mass_grid_msun)

    for index in range(n_halos):
        redshift = float(z_samples[index])
        hmf_values = np.asarray(hmf_callable(mass_grid_msun, redshift), dtype=float)
        selection_values = np.asarray(
            selection_callable(mass_grid_msun, np.full_like(mass_grid_msun, redshift)),
            dtype=float,
        )

        if hmf_values.shape != mass_grid_msun.shape:
            raise ValueError("hmf_callable must return same shape as mass grid.")
        if selection_values.shape != mass_grid_msun.shape:
            raise ValueError("selection callable must return same shape as mass grid.")
        if np.any(hmf_values < 0.0) or np.any(selection_values < 0.0):
            raise ValueError("hmf and selection values must be non-negative.")

        integrand = hmf_values * selection_values * mass_grid_msun
        if not np.any(integrand > 0.0):
            raise ValueError("Effective HMF x selection has zero support on mass grid.")

        cumulative = np.zeros_like(integrand)
        cumulative[1:] = np.cumsum(
            0.5 * (integrand[1:] + integrand[:-1]) * np.diff(log_mass_grid)
        )
        cumulative /= cumulative[-1]

        draw = random_generator.uniform(0.0, 1.0)
        sampled_log_mass = np.interp(draw, cumulative, log_mass_grid)
        samples[index] = np.exp(sampled_log_mass)

    return samples


def _sample_stellar_masses(
    n_halos: int,
    config_dict: dict[str, Any],
    random_generator: np.random.Generator,
) -> np.ndarray:
    stellar_mass_distribution = dict(
        config_dict.get(
            "stellar_mass_distribution",
            {"type": "lognormal", "mean_log10_mstar": 7.5, "sigma_log10_mstar": 0.5},
        )
    )
    distribution_type = str(stellar_mass_distribution.get("type", "lognormal")).lower()
    if distribution_type != "lognormal":
        raise ValueError("Only lognormal stellar_mass_distribution is supported in Tier-1.")

    mean_log10_mstar = float(stellar_mass_distribution["mean_log10_mstar"])
    sigma_log10_mstar = float(stellar_mass_distribution["sigma_log10_mstar"])
    log10_mstar_min = float(stellar_mass_distribution.get("log10_mstar_min", -np.inf))
    log10_mstar_max = float(stellar_mass_distribution.get("log10_mstar_max", np.inf))
    if sigma_log10_mstar < 0.0:
        raise ValueError("sigma_log10_mstar must be non-negative.")
    if log10_mstar_max <= log10_mstar_min:
        raise ValueError("log10_mstar_max must be greater than log10_mstar_min.")

    if sigma_log10_mstar == 0.0:
        log10_stellar_mass = np.full(n_halos, mean_log10_mstar)
    else:
        accepted: list[float] = []
        while len(accepted) < n_halos:
            draws = random_generator.normal(mean_log10_mstar, sigma_log10_mstar, size=n_halos)
            valid_draws = draws[(draws >= log10_mstar_min) & (draws <= log10_mstar_max)]
            accepted.extend(valid_draws.tolist())
        log10_stellar_mass = np.asarray(accepted[:n_halos], dtype=float)
    return 10.0**log10_stellar_mass


def _stellar_to_halo_mass(
    stellar_mass_msun: np.ndarray,
    config_dict: dict[str, Any],
    random_generator: np.random.Generator,
) -> np.ndarray:
    shmr_model = dict(
        config_dict.get(
            "shmr_model",
            {"type": "power_law", "a": 11.0, "b": 1.8, "pivot_log10_mstar": 8.0},
        )
    )
    model_type = str(shmr_model.get("type", "power_law")).lower()
    if model_type != "power_law":
        if model_type != "behroozi13":
            raise ValueError("Supported SHMR models are power_law and behroozi13.")

    log10_stellar_mass = np.log10(np.asarray(stellar_mass_msun, dtype=float))
    if model_type == "power_law":
        a = float(shmr_model["a"])
        b = float(shmr_model["b"])
        pivot_log10_mstar = float(shmr_model.get("pivot_log10_mstar", 8.0))
        mean_log10_halo_mass = a + b * (log10_stellar_mass - pivot_log10_mstar)
    else:
        # Lightweight Behroozi-inspired proxy tuned for dwarf-regime medians.
        pivot_log10_mstar = float(shmr_model.get("pivot_log10_mstar", 7.5))
        x_value = log10_stellar_mass - pivot_log10_mstar
        mean_log10_halo_mass = 10.0 + 1.25 * x_value + 0.10 * x_value**2

    halo_scatter_dex = float(config_dict.get("halo_scatter_dex", 0.3))
    if halo_scatter_dex < 0.0:
        raise ValueError("halo_scatter_dex must be non-negative.")

    sampled_log10_halo_mass = random_generator.normal(mean_log10_halo_mass, halo_scatter_dex)

    halo_mass_msun = 10.0**sampled_log10_halo_mass
    mhalo_min_msun = float(config_dict.get("mhalo_min_msun", 1.0e9))
    mhalo_max_msun = float(config_dict.get("mhalo_max_msun", 1.0e13))
    return np.clip(halo_mass_msun, mhalo_min_msun, mhalo_max_msun)


def _sample_concentrations(
    halo_mass_msun: np.ndarray,
    z_samples: np.ndarray,
    config_dict: dict[str, Any],
    random_generator: np.random.Generator,
) -> np.ndarray:
    concentration_scatter_dex = float(config_dict.get("concentration_scatter_dex", 0.15))
    if concentration_scatter_dex < 0.0:
        raise ValueError("concentration_scatter_dex must be non-negative.")

    concentration_model = dict(config_dict.get("concentration_model", {"type": "maccio"}))
    concentration_type = str(concentration_model.get("type", "maccio")).lower()

    if concentration_type == "maccio":
        concentration_mean = concentration_from_mass_maccio(halo_mass_msun, z_samples)
    elif concentration_type == "diemerjoyce19":
        try:
            from colossus.cosmology import cosmology as colossus_cosmology
            from colossus.halo import concentration
        except Exception:
            warnings.warn(
                (
                    "DiemerJoyce19 requested but colossus not available; "
                    "falling back to maccio relation."
                ),
                RuntimeWarning,
                stacklevel=2,
            )
            concentration_mean = concentration_from_mass_maccio(halo_mass_msun, z_samples)
        else:
            if "planck18" not in colossus_cosmology.cosmologies:
                colossus_cosmology.addCosmology(
                    "planck18",
                    {
                        "flat": True,
                        "H0": 67.4,
                        "Om0": 0.315,
                        "Ob0": 0.049,
                        "sigma8": 0.811,
                        "ns": 0.965,
                    },
                )
            colossus_cosmology.setCosmology("planck18")
            concentration_mean = concentration.concentration(
                M=np.asarray(halo_mass_msun, dtype=float),
                mdef="200c",
                z=np.asarray(z_samples, dtype=float),
                model="diemer19",
            )
            concentration_mean = np.asarray(concentration_mean, dtype=float)
            invalid = ~np.isfinite(concentration_mean) | (concentration_mean <= 0.0)
            if np.any(invalid):
                concentration_mean[invalid] = concentration_from_mass_maccio(
                    np.asarray(halo_mass_msun, dtype=float)[invalid],
                    np.asarray(z_samples, dtype=float)[invalid],
                )
    elif concentration_type == "power_law":
        c0 = float(concentration_model.get("c0", 4.5))
        m_pivot_msun = float(concentration_model.get("m_pivot_msun", 3.0e14))
        beta = float(concentration_model.get("beta", -0.1))
        gamma = float(concentration_model.get("gamma", -0.5))
        concentration_mean = (
            c0
            * (np.asarray(halo_mass_msun, dtype=float) / m_pivot_msun) ** beta
            * (1.0 + np.asarray(z_samples, dtype=float)) ** gamma
        )
    else:
        raise ValueError("Unsupported concentration_model type.")

    sampled_log10_concentration = random_generator.normal(
        np.log10(concentration_mean),
        concentration_scatter_dex,
    )
    return np.maximum(1.5, 10.0**sampled_log10_concentration)


def _compute_weights(
    halo_mass_msun: np.ndarray,
    z_samples: np.ndarray,
    config_dict: dict[str, Any],
) -> np.ndarray:
    if "weight_callable" in config_dict:
        raw_weights = np.asarray(
            config_dict["weight_callable"](halo_mass_msun, z_samples),
            dtype=float,
        )
    else:
        weight_mode = str(config_dict.get("weight_mode", "equal")).lower()
        if weight_mode == "equal":
            raw_weights = np.ones_like(halo_mass_msun, dtype=float)
        elif weight_mode == "mass":
            raw_weights = np.asarray(halo_mass_msun, dtype=float)
        else:
            raise ValueError("Unsupported weight_mode.")

    if np.any(raw_weights < 0.0):
        raise ValueError("Weights must be non-negative.")
    if not np.any(raw_weights > 0.0):
        raise ValueError("At least one weight must be positive.")
    return raw_weights / np.sum(raw_weights)


def generate_ensemble(
    mode: str,
    config_dict: dict[str, Any],
) -> list[dict[str, float]]:
    """Generate a halo ensemble for `mode in {'HMF', 'SHMR'}`."""
    mode_normalized = mode.strip().upper()
    if mode_normalized not in {"HMF", "SHMR"}:
        raise ValueError("mode must be 'HMF' or 'SHMR'.")

    n_halos = _validate_n_halos(config_dict)
    random_generator = _rng_from_config(config_dict)
    z_samples = _sample_redshifts(n_halos, config_dict, random_generator)

    if mode_normalized == "HMF":
        mass_grid_msun = _mass_grid_from_config(config_dict)
        hmf_callable = _resolve_hmf_callable(config_dict)
        selection_callable = _resolve_selection_callable(config_dict)
        halo_mass_msun = _sample_mass_from_hmf(
            n_halos=n_halos,
            mass_grid_msun=mass_grid_msun,
            z_samples=z_samples,
            hmf_callable=hmf_callable,
            selection_callable=selection_callable,
            random_generator=random_generator,
        )
        stellar_mass_msun = None
    else:
        stellar_mass_msun = _sample_stellar_masses(n_halos, config_dict, random_generator)
        halo_mass_msun = _stellar_to_halo_mass(stellar_mass_msun, config_dict, random_generator)

    concentration_samples = _sample_concentrations(
        halo_mass_msun=halo_mass_msun,
        z_samples=z_samples,
        config_dict=config_dict,
        random_generator=random_generator,
    )
    weights = _compute_weights(halo_mass_msun, z_samples, config_dict)

    halos: list[dict[str, float]] = []
    for index in range(n_halos):
        halo_record: dict[str, float] = {
            "M200": float(halo_mass_msun[index]),
            "m200_msun": float(halo_mass_msun[index]),
            "z": float(z_samples[index]),
            "c200": float(concentration_samples[index]),
            "weight": float(weights[index]),
        }
        if stellar_mass_msun is not None:
            halo_record["mstar_msun"] = float(stellar_mass_msun[index])
            halo_record["log10_mstar_msun"] = float(np.log10(stellar_mass_msun[index]))
        halos.append(halo_record)
    return halos


def sample_halo_ensemble(
    n_halos: int,
    mean_mass_msun: float = 3.0e14,
    mass_scatter_dex: float = 0.2,
    z: float = 0.4,
    concentration_scatter_dex: float = 0.1,
    seed: int | None = None,
) -> list[dict[str, float]]:
    """Backward-compatible lognormal-mass sampler built on `generate_ensemble`."""

    def lognormal_hmf(mass_msun: np.ndarray, redshift: float) -> np.ndarray:
        del redshift
        log_mass = np.log10(np.asarray(mass_msun, dtype=float))
        if mass_scatter_dex == 0.0:
            return np.exp(-1.0e6 * (log_mass - np.log10(mean_mass_msun)) ** 2)
        variance = mass_scatter_dex**2
        return np.exp(-0.5 * ((log_mass - np.log10(mean_mass_msun)) ** 2) / variance)

    return generate_ensemble(
        mode="HMF",
        config_dict={
            "n_halos": n_halos,
            "seed": seed,
            "redshift": z,
            "mass_min_msun": mean_mass_msun / 10.0,
            "mass_max_msun": mean_mass_msun * 10.0,
            "n_mass_grid": 512,
            "hmf_callable": lognormal_hmf,
            "selection_model": {"type": "none"},
            "concentration_scatter_dex": concentration_scatter_dex,
            "weight_mode": "equal",
        },
    )


def summarize_ensemble(halo_ensemble: list[dict[str, float]]) -> dict[str, float]:
    """Return compact ensemble summary statistics for validation and logging."""
    if len(halo_ensemble) == 0:
        raise ValueError("halo_ensemble cannot be empty.")

    mass = np.asarray([halo["m200_msun"] for halo in halo_ensemble], dtype=float)
    concentration = np.asarray([halo["c200"] for halo in halo_ensemble], dtype=float)
    weight = np.asarray([halo["weight"] for halo in halo_ensemble], dtype=float)

    if np.any(mass <= 0.0):
        raise ValueError("All sampled masses must be positive.")
    if np.any(concentration <= 0.0):
        raise ValueError("All sampled concentrations must be positive.")
    if np.any(weight < 0.0):
        raise ValueError("All weights must be non-negative.")

    summary = {
        "n_halos": float(len(halo_ensemble)),
        "mean_mass_msun": float(np.mean(mass)),
        "median_mass_msun": float(np.median(mass)),
        "std_log10_mass_dex": float(np.std(np.log10(mass))),
        "mean_c200": float(np.mean(concentration)),
        "median_c200": float(np.median(concentration)),
        "weight_sum": float(np.sum(weight)),
    }

    if "mstar_msun" in halo_ensemble[0]:
        stellar_mass = np.asarray([halo["mstar_msun"] for halo in halo_ensemble], dtype=float)
        summary["mean_stellar_mass_msun"] = float(np.mean(stellar_mass))
        summary["std_log10_stellar_mass_dex"] = float(np.std(np.log10(stellar_mass)))

    return summary
