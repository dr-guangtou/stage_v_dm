"""
Fisher matrix computation for SHMR parameter forecasts.

Computes the Fisher information matrix by accumulating contributions
from all (redshift bin, stellar mass bin) combinations. For each bin,
the observables are DeltaSigma(R), b_eff, and n_gal. Derivatives are
computed via central finite differences.

The Fisher matrix inverse gives the parameter covariance:
    sigma(theta_i) = sqrt((F^{-1})_{ii})

This module also provides utilities for marginalization, conditional
errors, Fisher ellipses, and external priors.

References
----------
Tegmark, Taylor & Heavens (1997), ApJ, 480, 22 — Fisher matrix formalism.
Oguri & Takada (2011), PRD, 83, 023008 — Fisher forecast for lensing.
"""

from __future__ import annotations

import warnings

import numpy as np

from .config import (
    ForecastConfig, LensingConfig, NuisanceConfig, SHMRParams, SurveyConfig,
)
from .covariance import (
    clustering_covariance,
    d_ln_inv_sigma_crit_d_dz,
    lensing_covariance,
    survey_volume,
)
from .halo_model import (
    delta_sigma_bin,
    effective_bias,
    galaxy_number_density,
)


# ---------------------------------------------------------------------------
# Parameter selection
# ---------------------------------------------------------------------------

def get_varied_params(
    params: SHMRParams,
    forecast_config: ForecastConfig,
    z_min: float,
    z_max: float,
) -> list[tuple[str, float]]:
    """
    Auto-select SHMR parameters to vary based on survey z-range.

    Always varies the 5 z=0 shape + scatter parameters. Also varies the
    4 evolution parameters if the survey spans enough redshift range AND
    the user has enabled evolution variation.

    Parameters
    ----------
    params : SHMRParams
        Fiducial SHMR parameters.
    forecast_config : ForecastConfig
        Forecast configuration (contains vary_z_evolution flag).
    z_min : float
        Minimum survey redshift.
    z_max : float
        Maximum survey redshift.

    Returns
    -------
    varied : list of (param_name, fiducial_value) tuples
    """
    # Always vary the z=0 shape parameters
    varied = [
        ('log_M1_0', params.log_M1_0),
        ('N_0', params.N_0),
        ('beta_0', params.beta_0),
        ('gamma_0', params.gamma_0),
    ]

    # Scatter parameters: mass-dependent (Cao & Tinker 2020) or constant
    if params.use_mass_dependent_scatter:
        varied.append(('scatter_sigma_high', params.scatter_sigma_high))
        varied.append(('scatter_sigma_rise', params.scatter_sigma_rise))
    else:
        varied.append(('sigma_logMs', params.sigma_logMs))

    # Vary evolution params only if survey spans enough redshift range
    if (forecast_config.vary_z_evolution
            and z_max > 0.4
            and (z_max - z_min) > 0.3):
        varied += [
            ('nu_M1', params.nu_M1),
            ('nu_N', params.nu_N),
            ('nu_beta', params.nu_beta),
            ('nu_gamma', params.nu_gamma),
        ]

    return varied


# ---------------------------------------------------------------------------
# Numerical derivatives
# ---------------------------------------------------------------------------

def compute_derivatives(
    observable_func: callable,
    params: SHMRParams,
    varied_params: list[tuple[str, float]],
    frac_step: float,
    **kwargs,
) -> np.ndarray:
    """
    Numerical derivatives of an observable w.r.t. SHMR parameters.

    Uses central finite differences:
        df/dtheta_i = [f(theta_i + delta) - f(theta_i - delta)] / (2*delta)
    where delta = frac_step * |theta_i|.

    Special case: if |theta_i| < 1e-10 (parameter near zero, e.g. nu_N),
    uses an absolute step of 0.01 instead.

    Parameters
    ----------
    observable_func : callable
        f(params, **kwargs) -> array of shape (N_obs,).
    params : SHMRParams
        Fiducial parameters.
    varied_params : list of (name, fiducial_value)
        Parameters to differentiate with respect to.
    frac_step : float
        Fractional step size.
    **kwargs
        Additional arguments passed to observable_func.

    Returns
    -------
    derivs : array, shape (N_params, N_obs)
        Jacobian matrix.
    """
    n_params = len(varied_params)

    # Evaluate fiducial to get output shape
    f0 = np.atleast_1d(observable_func(params, **kwargs))
    n_obs = len(f0)

    derivs = np.zeros((n_params, n_obs))

    for i, (name, fid_val) in enumerate(varied_params):
        # Step size: fractional or absolute for near-zero parameters
        if abs(fid_val) < 1e-10:
            delta = 0.01  # Absolute step for parameters near zero
        else:
            delta = frac_step * abs(fid_val)

        # Central finite difference
        params_plus = params.copy(**{name: fid_val + delta})
        params_minus = params.copy(**{name: fid_val - delta})

        f_plus = np.atleast_1d(observable_func(params_plus, **kwargs))
        f_minus = np.atleast_1d(observable_func(params_minus, **kwargs))

        derivs[i, :] = (f_plus - f_minus) / (2.0 * delta)

    return derivs


# ---------------------------------------------------------------------------
# Fisher matrix computation
# ---------------------------------------------------------------------------

def compute_fisher_matrix(
    shmr_params: SHMRParams,
    survey_config: SurveyConfig,
    lensing_config: LensingConfig,
    forecast_config: ForecastConfig,
    nuisance_config: NuisanceConfig | None = None,
) -> tuple[np.ndarray, list[str], dict]:
    """
    Main Fisher matrix computation.

    Accumulates Fisher information from lensing (DeltaSigma) and
    clustering (b_eff, n_gal) across all (z-bin, M*-bin) combinations.

    Supports two optional systematic error features (toggled via
    ForecastConfig):
    1. Fractional systematic floor on DeltaSigma covariance.
    2. Nuisance parameter marginalization (shear calibration m,
       source photo-z bias dz) with Gaussian priors.

    Parameters
    ----------
    shmr_params : SHMRParams
        Fiducial SHMR parameters.
    survey_config : SurveyConfig
        Spectroscopic survey configuration.
    lensing_config : LensingConfig
        Fixed lensing survey parameters.
    forecast_config : ForecastConfig
        Forecast behavior settings.
    nuisance_config : NuisanceConfig or None
        Nuisance parameter priors. Required if
        forecast_config.include_nuisance_params is True.

    Returns
    -------
    fisher : array, shape (N_params, N_params)
        Fisher information matrix.
    param_names : list of str
        Names of varied parameters (SHMR + nuisance if enabled).
    metadata : dict
        Diagnostic information (bins, N_lens, condition number, etc.).
    """
    # 1. Determine which SHMR parameters to vary
    varied = get_varied_params(
        shmr_params, forecast_config,
        survey_config.z_min, survey_config.z_max,
    )
    shmr_param_names = [name for name, _ in varied]
    n_shmr = len(varied)

    # Determine nuisance parameters
    nuisance_names: list[str] = []
    if forecast_config.include_nuisance_params:
        if nuisance_config is None:
            nuisance_config = NuisanceConfig()  # Use defaults
        nuisance_names = ['shear_m', 'photo_dz_source']

    n_nuisance = len(nuisance_names)
    n_total = n_shmr + n_nuisance
    all_param_names = shmr_param_names + nuisance_names

    # 2. Build redshift bins
    z_edges = np.arange(
        survey_config.z_min,
        survey_config.z_max + 0.001,
        forecast_config.dz,
    )
    if len(z_edges) < 2:
        z_edges = np.array([survey_config.z_min, survey_config.z_max])

    z_bins = [(z_edges[i], z_edges[i + 1]) for i in range(len(z_edges) - 1)]

    # 3. Radial bin edges for lensing
    R_bin_edges = np.logspace(
        np.log10(forecast_config.R_min_Mpc),
        np.log10(forecast_config.R_max_Mpc),
        forecast_config.n_R_bins + 1,
    )
    R_bin_centers = np.sqrt(R_bin_edges[:-1] * R_bin_edges[1:])

    # Halo mass grid for integrals
    log_Mh_grid = np.linspace(
        forecast_config.log_Mh_min,
        forecast_config.log_Mh_max,
        forecast_config.n_Mh_bins,
    )

    # 4. Initialize Fisher matrix and metadata
    fisher = np.zeros((n_total, n_total))
    metadata = {
        'z_bins': z_bins,
        'Mstar_bins_per_z': {},
        'N_lens_per_bin': {},
        'fiducial_observables': {},
        'param_names': all_param_names,
        'shmr_param_names': shmr_param_names,
        'nuisance_param_names': nuisance_names,
        'n_shmr': n_shmr,
        'n_nuisance': n_nuisance,
        'n_total': n_total,
        'systematic_floor_fraction': forecast_config.systematic_floor_fraction,
    }

    # Pre-compute total predicted galaxy counts per (z, M*) bin
    # to distribute N_gal_total proportionally
    bin_info = []
    total_predicted = 0.0

    for z_lo, z_hi in z_bins:
        z_mid = 0.5 * (z_lo + z_hi)
        try:
            mstar_bins = survey_config.get_stellar_mass_bins(z_mid)
        except ValueError:
            continue

        vol = survey_volume(z_lo, z_hi, survey_config.area_deg2)

        for ms_lo, ms_hi in mstar_bins:
            n_gal_vol = galaxy_number_density(
                ms_lo, ms_hi, shmr_params, z_mid,
                log_Mh_grid=log_Mh_grid,
            )
            n_predicted = n_gal_vol * vol
            bin_info.append((z_mid, z_lo, z_hi, ms_lo, ms_hi,
                             n_predicted, vol))
            total_predicted += n_predicted

    if total_predicted <= 0:
        warnings.warn(
            f"No galaxies predicted for survey '{survey_config.name}'. "
            "Returning zero Fisher matrix.",
            stacklevel=2,
        )
        return fisher, all_param_names, metadata

    # 5. Accumulate Fisher contributions from each bin
    f_sys = forecast_config.systematic_floor_fraction

    for z_mid, z_lo, z_hi, ms_lo, ms_hi, n_pred, vol in bin_info:
        # Distribute survey galaxies proportional to predicted counts
        N_lens = survey_config.n_gal_total * (n_pred / total_predicted)

        if N_lens < 1:
            continue

        # Store metadata
        metadata['Mstar_bins_per_z'].setdefault(z_mid, []).append(
            (ms_lo, ms_hi))
        metadata['N_lens_per_bin'][(z_mid, ms_lo)] = N_lens

        # --- Lensing observable: DeltaSigma(R) ---
        def ds_func(p, z=z_mid, lo=ms_lo, hi=ms_hi):
            return delta_sigma_bin(
                R_bin_centers, lo, hi, p, z,
                log_Mh_grid=log_Mh_grid,
            )

        ds_fiducial = ds_func(shmr_params)
        var_ds = lensing_covariance(
            R_bin_edges, z_mid, N_lens,
            survey_config.area_deg2, lensing_config,
        )

        # Apply systematic floor: irreducible fractional uncertainty
        # on DeltaSigma from baryonic effects, miscentering, etc.
        if f_sys > 0:
            var_sys = (f_sys * ds_fiducial) ** 2
            var_ds = var_ds + var_sys

        # SHMR derivatives of DeltaSigma
        J_ds_shmr = compute_derivatives(
            ds_func, shmr_params, varied,
            forecast_config.frac_step,
        )  # shape (n_shmr, n_R_bins)

        # Build full Jacobian including nuisance parameters
        if n_nuisance > 0:
            J_ds_nuisance = np.zeros((n_nuisance, len(R_bin_centers)))

            # d(DS_obs)/dm = DS_true: shear calibration bias
            # At fiducial m=0: DS_obs = (1+m) * DS_true
            J_ds_nuisance[0, :] = ds_fiducial

            # d(DS_obs)/d(dz_s) = DS_true * d(ln Sigma_crit_eff^{-1})/d(dz_s)
            # Photo-z bias shifts the source n(z), modifying Sigma_crit_eff
            dln_inv_sc = d_ln_inv_sigma_crit_d_dz(z_mid, lensing_config)
            J_ds_nuisance[1, :] = ds_fiducial * dln_inv_sc

            J_ds_full = np.vstack([J_ds_shmr, J_ds_nuisance])
        else:
            J_ds_full = J_ds_shmr

        # Fisher contribution from lensing: F += J^T C^{-1} J
        valid_ds = (var_ds > 0) & np.isfinite(var_ds) & (ds_fiducial > 0)
        if np.any(valid_ds):
            inv_var_ds = np.zeros_like(var_ds)
            inv_var_ds[valid_ds] = 1.0 / var_ds[valid_ds]
            fisher += (
                J_ds_full[:, valid_ds]
                @ np.diag(inv_var_ds[valid_ds])
                @ J_ds_full[:, valid_ds].T
            )

        # --- Clustering observables: b_eff and n_gal ---
        # Not affected by lensing nuisance parameters
        def beff_func(p, z=z_mid, lo=ms_lo, hi=ms_hi):
            return np.array([effective_bias(
                lo, hi, p, z, log_Mh_grid=log_Mh_grid,
            )])

        def ngal_func(p, z=z_mid, lo=ms_lo, hi=ms_hi):
            return np.array([galaxy_number_density(
                lo, hi, p, z, log_Mh_grid=log_Mh_grid,
            )])

        b_fid = beff_func(shmr_params)[0]
        n_fid = ngal_func(shmr_params)[0]

        # Clustering covariance uses MODEL n_gal (not N_lens/V)
        var_b, var_n = clustering_covariance(n_fid, b_fid, vol, z_mid)

        # SHMR derivatives of clustering observables
        J_b_shmr = compute_derivatives(
            beff_func, shmr_params, varied,
            forecast_config.frac_step,
        )
        J_n_shmr = compute_derivatives(
            ngal_func, shmr_params, varied,
            forecast_config.frac_step,
        )

        # Zero-pad clustering Jacobians for nuisance rows
        # (clustering is unaffected by lensing nuisance params)
        if n_nuisance > 0:
            J_b_full = np.zeros((n_total, 1))
            J_b_full[:n_shmr, :] = J_b_shmr
            J_n_full = np.zeros((n_total, 1))
            J_n_full[:n_shmr, :] = J_n_shmr
        else:
            J_b_full = J_b_shmr
            J_n_full = J_n_shmr

        # Fisher contribution from b_eff
        if var_b > 0 and np.isfinite(var_b) and b_fid > 0:
            fisher += (J_b_full @ J_b_full.T) / var_b

        # Fisher contribution from n_gal
        if var_n > 0 and np.isfinite(var_n) and n_fid > 0:
            fisher += (J_n_full @ J_n_full.T) / var_n

        # Store fiducial observables
        metadata['fiducial_observables'][(z_mid, ms_lo)] = {
            'ds': ds_fiducial,
            'b_eff': b_fid,
            'n_gal': n_fid,
            'var_ds': var_ds,
            'var_b': var_b,
            'var_n': var_n,
            'N_lens': N_lens,
        }

    # 6. Add Gaussian priors on nuisance parameters
    if n_nuisance > 0 and nuisance_config is not None:
        fisher = add_external_prior(
            fisher, all_param_names, 'shear_m', nuisance_config.sigma_m)
        fisher = add_external_prior(
            fisher, all_param_names, 'photo_dz_source',
            nuisance_config.sigma_dz_source)

    # 7. Check condition number
    if np.any(fisher != 0):
        cond = np.linalg.cond(fisher)
        metadata['condition_number'] = cond
        if cond > 1e10:
            warnings.warn(
                f"Fisher matrix condition number is {cond:.2e} (> 1e10). "
                "Consider fixing evolution parameters or adding priors.",
                stacklevel=2,
            )
    else:
        metadata['condition_number'] = np.inf

    return fisher, all_param_names, metadata


# ---------------------------------------------------------------------------
# Post-processing
# ---------------------------------------------------------------------------

def marginalized_errors(fisher_matrix: np.ndarray) -> np.ndarray:
    """
    Compute 1-sigma marginalized errors from the Fisher matrix.

    sigma(theta_i) = sqrt((F^{-1})_{ii})

    Parameters
    ----------
    fisher_matrix : array, shape (N, N)
        Fisher information matrix.

    Returns
    -------
    errors : array, shape (N,)
        Marginalized 1-sigma errors.
    """
    cov = np.linalg.inv(fisher_matrix)
    return np.sqrt(np.diag(cov))


def conditional_errors(fisher_matrix: np.ndarray) -> np.ndarray:
    """
    Compute conditional (unmarginalized) errors.

    sigma_cond(theta_i) = 1 / sqrt(F_{ii})

    These are always smaller than marginalized errors.

    Parameters
    ----------
    fisher_matrix : array, shape (N, N)

    Returns
    -------
    errors : array, shape (N,)
    """
    return 1.0 / np.sqrt(np.diag(fisher_matrix))


def fisher_ellipse(
    fisher_matrix: np.ndarray,
    i: int,
    j: int,
    n_sigma: float = 1.0,
    n_points: int = 100,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute the marginalized Fisher ellipse for parameter pair (i, j).

    Parameters
    ----------
    fisher_matrix : array, shape (N, N)
    i, j : int
        Parameter indices.
    n_sigma : float
        Number of sigma for the ellipse contour.
    n_points : int
        Number of points on the ellipse.

    Returns
    -------
    x, y : arrays, shape (n_points,)
        Ellipse coordinates centered on zero (add fiducial values to plot).
    """
    cov = np.linalg.inv(fisher_matrix)
    # Extract 2x2 submatrix
    sub_cov = np.array([
        [cov[i, i], cov[i, j]],
        [cov[j, i], cov[j, j]],
    ])

    # Eigendecomposition
    eigvals, eigvecs = np.linalg.eigh(sub_cov)

    # Parametric ellipse
    theta = np.linspace(0, 2 * np.pi, n_points)
    ellipse = np.array([
        n_sigma * np.sqrt(eigvals[0]) * np.cos(theta),
        n_sigma * np.sqrt(eigvals[1]) * np.sin(theta),
    ])

    # Rotate to parameter space
    rotated = eigvecs @ ellipse
    return rotated[0], rotated[1]


def add_external_prior(
    fisher: np.ndarray,
    param_names: list[str],
    param_name: str,
    sigma_prior: float,
) -> np.ndarray:
    """
    Add a Gaussian prior on a single parameter.

    F[i,i] += 1/sigma_prior^2

    Parameters
    ----------
    fisher : array, shape (N, N)
        Fisher matrix (modified in place AND returned).
    param_names : list of str
    param_name : str
        Name of parameter to add prior to.
    sigma_prior : float
        1-sigma prior width.

    Returns
    -------
    fisher : array, shape (N, N)
        Modified Fisher matrix.
    """
    idx = param_names.index(param_name)
    fisher[idx, idx] += 1.0 / sigma_prior**2
    return fisher


def extract_shmr_constraints(
    fisher: np.ndarray,
    param_names: list[str],
) -> tuple[np.ndarray, list[str]]:
    """
    Extract SHMR-only marginalized errors after marginalizing over nuisance params.

    Inverts the full Fisher matrix (which includes nuisance parameters
    with priors), then extracts the SHMR sub-block of the covariance.

    Parameters
    ----------
    fisher : array, shape (N_total, N_total)
        Full Fisher matrix including nuisance parameters.
    param_names : list of str
        All parameter names.

    Returns
    -------
    shmr_errors : array, shape (N_shmr,)
        Marginalized 1-sigma errors on SHMR parameters only.
    shmr_names : list of str
        SHMR parameter names.
    """
    nuisance_set = {'shear_m', 'photo_dz_source'}
    shmr_idx = [i for i, n in enumerate(param_names) if n not in nuisance_set]
    shmr_names = [param_names[i] for i in shmr_idx]

    # Full covariance -> extract SHMR sub-block
    cov = np.linalg.inv(fisher)
    shmr_errors = np.sqrt(np.diag(cov)[shmr_idx])

    return shmr_errors, shmr_names
