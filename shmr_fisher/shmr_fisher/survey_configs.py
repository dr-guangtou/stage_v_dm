"""
Predefined survey configurations spanning Stage-III through Stage-V.

These are generic archetypes representing the constraining power of
different survey generations. They are NOT tied to specific named
surveys — the user can modify any parameter or create custom configs
via make_custom_survey().

Each survey is defined by its footprint (area, z-range), galaxy sample
(total count, mass completeness), and stellar mass binning.

References
----------
PLAN.md Task 10 for the design rationale behind each configuration.
"""

from __future__ import annotations

from .config import SurveyConfig


# ---------------------------------------------------------------------------
# Stage III: existing spectroscopic surveys (e.g., SDSS/GAMA-like)
# ---------------------------------------------------------------------------

stage3_shallow_wide = SurveyConfig(
    name="Stage-III Shallow Wide",
    area_deg2=7500,
    z_min=0.02,
    z_max=0.2,
    n_gal_total=700_000,
    log_Mstar_min=9.5,
)

# ---------------------------------------------------------------------------
# Stage IV: current-generation surveys (e.g., DESI-like)
# ---------------------------------------------------------------------------

stage4_low_z = SurveyConfig(
    name="Stage-IV Low-z",
    area_deg2=14000,
    z_min=0.05,
    z_max=0.4,
    n_gal_total=10_000_000,
    log_Mstar_min=9.0,
)

stage4_high_z = SurveyConfig(
    name="Stage-IV High-z",
    area_deg2=14000,
    z_min=0.4,
    z_max=1.0,
    n_gal_total=5_000_000,
    log_Mstar_min=10.8,
)

# ---------------------------------------------------------------------------
# Stage V: next-generation wide-field spectroscopic surveys
# ---------------------------------------------------------------------------

stage5_wide = SurveyConfig(
    name="Stage-V Wide",
    area_deg2=10000,
    z_min=0.05,
    z_max=1.0,
    n_gal_total=50_000_000,
    log_Mstar_min=8.5,
    # Mass completeness degrades with z: 8.5 at z=0, 9.5 at z=1
    log_Mstar_min_func=lambda z: 8.5 + 1.0 * z,
)

stage5_deep = SurveyConfig(
    name="Stage-V Deep",
    area_deg2=3000,
    z_min=0.1,
    z_max=1.5,
    n_gal_total=20_000_000,
    log_Mstar_min=8.0,
    # 8.0 at z=0, 9.5 at z=1.5
    log_Mstar_min_func=lambda z: 8.0 + 1.0 * z,
)

# ---------------------------------------------------------------------------
# Survey dictionary for batch processing
# ---------------------------------------------------------------------------

surveys: dict[str, SurveyConfig] = {
    "stage3_shallow_wide": stage3_shallow_wide,
    "stage4_low_z": stage4_low_z,
    "stage4_high_z": stage4_high_z,
    "stage5_wide": stage5_wide,
    "stage5_deep": stage5_deep,
}


# ---------------------------------------------------------------------------
# Factory for custom surveys
# ---------------------------------------------------------------------------

def make_custom_survey(
    name: str,
    area_deg2: float,
    z_min: float,
    z_max: float,
    n_gal_total: float,
    log_Mstar_min: float,
    dlog_Mstar: float = 0.5,
    **kwargs,
) -> SurveyConfig:
    """
    Quick factory for creating a custom survey configuration.

    Useful for parameter sweeps where a single property is varied.

    Parameters
    ----------
    name : str
        Human-readable label.
    area_deg2 : float
        Sky area [deg^2].
    z_min, z_max : float
        Redshift range.
    n_gal_total : float
        Total number of spec-z galaxies.
    log_Mstar_min : float
        Stellar mass completeness floor [log10(Msun)].
    dlog_Mstar : float
        Stellar mass bin width [dex]. Default 0.5.
    **kwargs
        Additional SurveyConfig fields (e.g., log_Mstar_min_func).

    Returns
    -------
    survey : SurveyConfig
    """
    return SurveyConfig(
        name=name,
        area_deg2=area_deg2,
        z_min=z_min,
        z_max=z_max,
        n_gal_total=n_gal_total,
        log_Mstar_min=log_Mstar_min,
        dlog_Mstar=dlog_Mstar,
        **kwargs,
    )
