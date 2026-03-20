"""
SHMR Forecast: Fisher matrix forecasts for stellar-halo mass relation constraints.

This package computes how well spectroscopic surveys (combined with Stage-IV
imaging lensing) can constrain the galaxy-halo connection, using an analytic
halo model built on the Moster+2013 SHMR parameterization.

Modules
-------
config : Configuration dataclasses (SHMRParams, SurveyConfig, etc.)
config_io : YAML-based run configuration loader
shmr_model : Parameterized SHMR with redshift evolution
halo_model : Observable predictions (DeltaSigma, b_eff, n_gal)
covariance : Analytic noise models for lensing and clustering
fisher : Fisher matrix computation and marginalization
survey_configs : Predefined Stage-III/IV/V survey configurations
validate : Validation checks for each pipeline phase
plot_results : Diagnostic and science figures
"""

from colossus.cosmology import cosmology

# Set fiducial cosmology at import time — all modules use this
cosmology.setCosmology('planck18')
