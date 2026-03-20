"""
Configuration dataclasses for the SHMR forecast pipeline.

This module defines all configurable parameters for the forecast:
- SHMRParams: the 9-parameter Moster+2013 SHMR model
- SurveyConfig: spectroscopic survey specification
- LensingConfig: fixed Stage-IV lensing assumptions
- ForecastConfig: controls for the Fisher matrix computation

All physical quantities carry explicit units in their docstrings.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields, replace
from typing import Callable

import numpy as np


@dataclass
class SHMRParams:
    """
    Parameterized stellar-halo mass relation following Moster+2013.

    The mean SHMR is a double power-law for M*/Mh as a function of Mh,
    with linear redshift evolution in z/(1+z):

        M*(Mh, z) = 2 N(z) Mh [(Mh/M1(z))^(-beta(z)) + (Mh/M1(z))^(gamma(z))]^(-1)

    Redshift evolution (Moster+2013 Eq. 11-14):
        log10(M1(z)) = log_M1_0 + nu_M1 * z/(1+z)
        N(z)         = N_0      + nu_N  * z/(1+z)
        beta(z)      = beta_0   + nu_beta * z/(1+z)
        gamma(z)     = gamma_0  + nu_gamma * z/(1+z)

    Attributes
    ----------
    log_M1_0 : float
        Characteristic halo mass at z=0 [log10(Msun)].
    N_0 : float
        Peak stellar-to-halo mass ratio at z=0 [dimensionless].
    beta_0 : float
        Low-mass power-law slope at z=0 [dimensionless].
    gamma_0 : float
        High-mass power-law slope at z=0 [dimensionless].
    nu_M1 : float
        Redshift evolution of log M1 [dimensionless].
    nu_N : float
        Redshift evolution of N [dimensionless].
    nu_beta : float
        Redshift evolution of beta [dimensionless].
    nu_gamma : float
        Redshift evolution of gamma [dimensionless].
    sigma_logMs : float
        Log-normal scatter in M* at fixed Mh [dex], constant with z.

    Notes
    -----
    Fiducial values from Moster+2013 Table 1.
    """

    # z=0 SHMR shape
    log_M1_0: float = 11.59
    N_0: float = 0.0351
    beta_0: float = 1.376
    gamma_0: float = 0.608

    # Redshift evolution
    nu_M1: float = 1.195
    nu_N: float = -0.0247
    nu_beta: float = -0.826
    nu_gamma: float = 0.329

    # Scatter
    # When use_mass_dependent_scatter is False, a constant scatter is used:
    sigma_logMs: float = 0.15

    # Mass-dependent scatter following Cao & Tinker (2020, MNRAS, 498, 5080).
    # sigma(log M* | Mh) = sigma_high + sigma_rise * [1 - tanh((log_Mh - log_Mh_break) / delta)]
    # At high Mh (>> break): sigma -> sigma_high
    # At low  Mh (<< break): sigma -> sigma_high + 2 * sigma_rise
    use_mass_dependent_scatter: bool = False
    scatter_sigma_high: float = 0.18    # Asymptotic scatter at high Mh [dex]
    scatter_sigma_rise: float = 0.10    # Half-amplitude of the rise toward low Mh [dex]
    scatter_log_Mh_break: float = 12.0  # Transition halo mass [log10(Msun)]
    scatter_delta: float = 0.4          # Transition width [dex]

    def to_dict(self) -> dict[str, float]:
        """Return all parameters as a {name: value} dictionary.

        Returns
        -------
        d : dict
            Parameter names mapped to their values.
        """
        return {f.name: getattr(self, f.name) for f in fields(self)}

    @classmethod
    def from_dict(cls, d: dict[str, float]) -> SHMRParams:
        """Construct an SHMRParams instance from a dictionary.

        Parameters
        ----------
        d : dict
            Dictionary with parameter names as keys.

        Returns
        -------
        params : SHMRParams
        """
        return cls(**{k: v for k, v in d.items() if k in {f.name for f in fields(cls)}})

    def copy(self, **overrides: float) -> SHMRParams:
        """Return a new instance with specified parameters replaced.

        This is used for finite-difference derivatives in the Fisher matrix.

        Parameters
        ----------
        **overrides : float
            Parameter names and their new values.

        Returns
        -------
        params : SHMRParams
            New instance with overrides applied.
        """
        return replace(self, **overrides)


@dataclass
class SurveyConfig:
    """
    Spectroscopic survey parameters -- fully generic.

    A survey is defined by its area, redshift range, and a description of
    the galaxy sample (number density + mass completeness) which CAN vary
    with redshift.

    Attributes
    ----------
    name : str
        Human-readable survey label.
    area_deg2 : float
        Sky area [deg^2].
    z_min : float
        Minimum redshift.
    z_max : float
        Maximum redshift.
    n_gal_total : float
        Total number of spec-z galaxies in the survey.
    log_Mstar_min : float
        Global stellar mass completeness floor [log10(Msun)].
    log_Mstar_max : float
        Upper edge of highest stellar mass bin [log10(Msun)].
    dlog_Mstar : float
        Stellar mass bin width [dex].
    log_Mstar_min_func : callable or None
        Optional z-dependent completeness: f(z) -> log10(Msun).
        If None, uses the constant log_Mstar_min at all redshifts.
    """

    name: str
    area_deg2: float
    z_min: float
    z_max: float
    n_gal_total: float
    log_Mstar_min: float
    log_Mstar_max: float = 12.0
    dlog_Mstar: float = 0.5
    log_Mstar_min_func: Callable[[float], float] | None = field(
        default=None, repr=False
    )

    def get_log_Mstar_min(self, z: float) -> float:
        """Return the stellar mass completeness limit at redshift z.

        Parameters
        ----------
        z : float
            Redshift.

        Returns
        -------
        log_Mstar_min : float
            Completeness limit [log10(Msun)].
        """
        if self.log_Mstar_min_func is not None:
            return self.log_Mstar_min_func(z)
        return self.log_Mstar_min

    def get_n_gal_per_deg2(self) -> float:
        """Return the average galaxy surface density [galaxies/deg^2].

        Returns
        -------
        n : float
            Average number of galaxies per square degree.
        """
        return self.n_gal_total / self.area_deg2

    def get_stellar_mass_bins(self, z: float | None = None) -> list[tuple[float, float]]:
        """Return stellar mass bin edges for a given redshift.

        Bin edges run from get_log_Mstar_min(z) to log_Mstar_max in steps
        of dlog_Mstar.

        Parameters
        ----------
        z : float or None
            Redshift for z-dependent completeness. If None, uses the
            global log_Mstar_min.

        Returns
        -------
        bins : list of (lo, hi) tuples
            Stellar mass bin edges [log10(Msun)].

        Raises
        ------
        ValueError
            If no bins can be formed (completeness limit >= upper bound).
        """
        mmin = self.get_log_Mstar_min(z) if z is not None else self.log_Mstar_min
        edges = np.arange(mmin, self.log_Mstar_max + 0.01, self.dlog_Mstar)
        bins = [(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]
        if len(bins) == 0:
            raise ValueError(
                f"No stellar mass bins: completeness limit {mmin:.1f} >= "
                f"upper bound {self.log_Mstar_max:.1f}"
            )
        return bins


@dataclass
class LensingConfig:
    """
    Fixed lensing survey parameters (LSST-like / Euclid-like).

    These are held constant across all forecast comparisons --
    we only vary the spectroscopic survey properties.

    Attributes
    ----------
    n_source_per_arcmin2 : float
        Effective source galaxy density [arcmin^-2].
    sigma_e : float
        Shape noise per ellipticity component [dimensionless].
    z_source_median : float
        Median source redshift.

    Notes
    -----
    Source redshift distribution follows a Smail-type form:
        n(z) ~ z^2 exp(-(z/z0)^1.5)
    with z0 = z_source_median / 1.41.
    """

    n_source_per_arcmin2: float = 25.0
    sigma_e: float = 0.26
    z_source_median: float = 1.0


@dataclass
class NuisanceConfig:
    """
    Configuration for lensing nuisance parameters.

    These parameters capture systematic uncertainties in the lensing
    measurement that are not part of the SHMR model itself. They are
    marginalized over in the Fisher matrix with Gaussian priors.

    Attributes
    ----------
    sigma_m : float
        Gaussian prior width on multiplicative shear calibration bias m
        [dimensionless]. The observed DeltaSigma is (1+m) * DeltaSigma_true.
        Fiducial m=0. Typical values: 0.01 (optimistic) to 0.05 (conservative).
    sigma_dz_source : float
        Gaussian prior width on source photo-z bias delta_z [dimensionless].
        A positive delta_z shifts the source n(z) to higher redshift,
        decreasing Sigma_crit_eff. Fiducial delta_z=0.
        Typical values: 0.01 (optimistic) to 0.05 (conservative).
    """

    sigma_m: float = 0.02
    sigma_dz_source: float = 0.03


@dataclass
class ForecastConfig:
    """
    Master configuration controlling the Fisher forecast behavior.

    Attributes
    ----------
    R_min_Mpc : float
        Minimum projected radius for lensing radial bins [physical Mpc].
    R_max_Mpc : float
        Maximum projected radius for lensing radial bins [physical Mpc].
    n_R_bins : int
        Number of log-spaced radial bins.
    dz : float
        Default redshift bin width.
    frac_step : float
        Fractional step size for numerical derivatives.
    vary_z_evolution : bool
        If False, fix nu_* parameters (for low-z-only surveys).
    log_Mh_min : float
        Minimum halo mass for integration [log10(Msun)].
    log_Mh_max : float
        Maximum halo mass for integration [log10(Msun)].
    n_Mh_bins : int
        Number of log-spaced halo mass bins for integration.
    systematic_floor_fraction : float
        Fractional systematic floor on DeltaSigma covariance. The systematic
        variance in each radial bin is (f_sys * DeltaSigma_fid)^2, added in
        quadrature with the statistical variance. Set to 0.0 to disable
        (default). Typical values: 0.05-0.10.
    include_nuisance_params : bool
        If True, include lensing nuisance parameters (shear calibration m,
        source photo-z bias dz) in the Fisher matrix and marginalize over
        them with Gaussian priors from NuisanceConfig.
    fixed_params : list of str
        SHMR parameter names to hold fixed (not varied in the Fisher matrix).
        Useful for focused science cases, e.g., fixing gamma_0 and log_M1_0
        for dwarf-galaxy forecasts where the high-mass SHMR is irrelevant.
        Overrides the default parameter selection in get_varied_params().
        Empty list (default) means no additional parameters are fixed beyond
        the standard logic (e.g., nu_* when vary_z_evolution is False).
    sigma_log_Mstar_obs : float
        Statistical uncertainty in observed log M* [dex]. When > 0, this
        observational scatter is added in quadrature to the intrinsic SHMR
        scatter, broadening the HOD occupation functions. This mimics the
        effect of photometric stellar mass errors on bin assignments.
        Default 0.0 (no observational scatter).
    include_smf : bool
        If True, include the stellar mass function (SMF) as a data vector
        in the Fisher matrix, using Poisson + cosmic variance covariance
        instead of Poisson-only for the n_gal contribution. Default True.
    """

    R_min_Mpc: float = 0.1
    R_max_Mpc: float = 30.0
    n_R_bins: int = 10
    dz: float = 0.2
    frac_step: float = 0.01
    vary_z_evolution: bool = True
    log_Mh_min: float = 10.0
    log_Mh_max: float = 15.5
    n_Mh_bins: int = 200
    systematic_floor_fraction: float = 0.0
    include_nuisance_params: bool = False
    fixed_params: list[str] = field(default_factory=list)
    sigma_log_Mstar_obs: float = 0.0
    include_smf: bool = True
