# SIDM vs CDM Stage-V Lensing Forecast

Minimal, back-of-the-envelope forecasting framework for assessing whether a future Stage-V spectroscopic redshift survey, combined with weak lensing, can distinguish cold dark matter (CDM) from self-interacting dark matter (SIDM).

This project is intentionally scoped for a ~1 day implementation and a science-case / workshop-level result, not a publication-grade end-to-end forecast.

## Goal

Build a fast pipeline that:

1. generates 3D halo density profiles for CDM and SIDM,
2. projects them into surface density `Sigma(R)` and excess surface density `DeltaSigma(R)`,
3. compares the signals for a small set of benchmark SIDM cross sections, and
4. reports a simple distinguishability metric such as `Delta chi^2` or the radial precision required for detection.

The two benchmark halo regimes are:

- Dwarf halo: `M200 = 1e10 Msun`
- Cluster halo: `M200 = 1e14 Msun`

## Philosophy

This forecast is designed to answer a narrow question quickly:

> If I assume a specific SIDM parameterization, how large is the induced change in `rho(r)` and `DeltaSigma(R)` at dwarf and cluster scales, and what measurement precision would be needed to separate it from CDM?

The forecast should therefore:

- use a specific SIDM parameterization, not a vague “SIDM vs CDM” label;
- stay close to the observable (`DeltaSigma`), while preserving the underlying 3D profile outputs for sanity checks;
- avoid a full HOD / survey systematics model in version 1;
- be easy to hand off to an autonomous coding agent.

## Recommended modeling stack

### Core SIDM engine

Use Daneng Yang’s `parametricSIDM` repository as the engine for generating SIDM 3D density profiles from a CDM halo description.

Repository: <https://github.com/DanengYang/parametricSIDM>

### Lensing / profile utilities

Do not depend on `cluster_toolkit` in the implementation.

Preferred approach:

- implement the projection integrals directly in NumPy/SciPy;
- optionally use CCL or Colossus as references or cross-checks for conventions and normalization;
- treat halomod as optional future infrastructure for 2-halo or HOD-like extensions.

Useful references:

- CCL halo-profile API and projected / cumulative 2D quantities
- Colossus halo-profile machinery, including `surface density` and `DeltaSigma`
- CLMM for weak-lensing observable conventions and unit handling

## Benchmark parameterization

Use an effective halo-scale SIDM cross section as the first parameterization:

- `sigma_over_m_eff` in `cm^2 / g`

For the first-pass forecast, benchmark grids can be:

### Dwarf halo (`1e10 Msun`)
- `0.0`
- `0.2`
- `0.5`
- `1.0`
- `2.0`

### Cluster halo (`1e14 Msun`)
- `0.0`
- `0.2`
- `0.5`
- `1.0`
- `2.0`

These are not meant to be the final particle-physics parameterization. They are a practical first layer that can later be mapped onto a velocity-dependent SIDM model.

## Minimal deliverables

### 1. Profile generator

A module that returns, for a given halo mass, concentration, redshift, and SIDM cross section:

- `r`
- `rho_3d(r)`
- `M_enclosed(<r)`
- `Vc(r)`

### 2. Projection module

A module that numerically computes:

- `Sigma(R)`
- `Sigma_bar(<R)`
- `DeltaSigma(R) = Sigma_bar(<R) - Sigma(R)`

### 3. Toy forecast module

A module that defines a simple measurement model and computes:

- per-bin uncertainties,
- `Delta chi^2` between CDM and each SIDM benchmark,
- required fractional precision to reach a target significance.

### 4. Output products

At minimum, the code should produce:

- 3D density-profile plots,
- `DeltaSigma(R)` plots,
- ratio plots relative to CDM,
- a summary table of benchmark models and distinguishability.

## Scope control: what version 1 should include

Include:

- isolated halo profiles,
- simple NFW CDM baselines,
- direct numerical projection to `DeltaSigma`,
- toy radial bins and toy errors,
- dwarf and cluster benchmark cases.

Do not include in version 1 unless trivial:

- halo mass scatter,
- central/satellite mixing,
- miscentering,
- 2-halo term,
- baryons beyond a simple optional placeholder,
- source photo-z systematics,
- intrinsic alignments,
- survey masking / selection realism.

## Suggested radial ranges

These are starting points, not hard requirements.

### Dwarf halo
- 3D profile radius: `0.1 - 300 kpc`
- lensing radius: `3 - 300 kpc`

### Cluster halo
- 3D profile radius: `5 - 5000 kpc`
- lensing radius: `30 - 3000 kpc`

Use logarithmic bins.

## Suggested toy observational assumptions

These are workshop-level placeholders.

### Dwarf lensing

Use a somewhat optimistic first-pass fractional uncertainty model, for example:

- baseline: `10-20%` fractional error per radial bin over the best-constrained radial range,
- plus an alternate conservative curve.

### Cluster lensing

Use a more mature weak-lensing precision benchmark, for example:

- baseline: `3-10%` fractional error per radial bin,
- plus an alternate conservative curve.

The main point is not to claim survey realism, but to translate theory-level profile differences into an interpretable measurement requirement.

## Primary science outputs

The most useful summary statements are expected to be of the form:

- “At `M200 = 1e10 Msun`, SIDM with `sigma/m = X` changes `DeltaSigma` by ~Y% over radii `R1-R2`.”
- “At `M200 = 1e14 Msun`, distinguishing CDM from SIDM with `sigma/m = X` requires roughly Z% precision per radial bin over `R1-R2`.”
- “Dwarf halos provide larger intrinsic SIDM signal, while clusters provide a more observationally accessible `DeltaSigma` measurement.”

## Repository layout

Suggested minimal layout:

```text
sidm_stagev_forecast/
  README.md
  SPEC.md
  pyproject.toml
  requirements.txt
  src/
    sidm_stagev/
      __init__.py
      config.py
      cosmology.py
      profiles.py
      projection.py
      forecast.py
      plotting.py
      io.py
  notebooks/
    01_dwarf_cluster_forecast.ipynb
  scripts/
    run_benchmarks.py
  outputs/
    figures/
    tables/
```

## Implementation notes

- Keep all units explicit.
- Prefer physical radii for the forecast products unless there is a strong reason not to.
- Keep the code modular enough that a later version can swap in baryons, the 2-halo term, or a more realistic covariance model.
- Avoid hidden conventions.
- Save intermediate products (e.g. profile tables) for debugging and reproducibility.

## References

1. Daneng Yang et al., A Parametric Model for Self-Interacting Dark Matter Halos — public code and method description.
2. Daneng Yang GitHub note pointing to the separate SIDM_Lensing_Model for lensing-specific applications.
3. CCL documentation for projected and cumulative 2D halo-profile quantities.
4. Colossus documentation for halo density, surface density, and `DeltaSigma` profile utilities.
5. CLMM documentation / paper for weak-lensing observable definitions and unit conventions.

---

## Tier-1: Stacked ΔΣ(R) from a Halo Ensemble

This project implements a first-order forecasting pipeline to assess whether a future Stage-V spectroscopic survey combined with weak lensing can distinguish Cold Dark Matter (CDM) from Self-Interacting Dark Matter (SIDM).

This version upgrades the original single-fiducial-halo approach to a **halo ensemble stacking framework**, motivated by recent simulation-based SIDM analyses where constraints are derived from stacked cluster weak-lensing profiles.

---

# Scientific Motivation

Recent SIDM simulation studies show:

- SIDM modifies inner halo density structure.
- Observational constraints are derived from **stacked ΔΣ(R)**.
- Results depend on full halo mass distributions.

Therefore, we now compute:

Halo Ensemble → ρ(r) → Σ(R) → ΔΣ(R) → Stack → Forecast

This is a controlled, fast, transparent science-case forecast.

---

# Scope of Tier-1

We model:

- A halo ensemble drawn from a mass distribution
- Individual halo ρ(r) using:
  - CDM NFW baseline
  - SIDM profiles via `parametricSIDM`
- Numerical projection to Σ(R) and ΔΣ(R)
- Stacked ΔΣ(R)
- Toy Δχ² distinguishability metric

We explicitly DO NOT include:

- Splashback modeling
- 2-halo term
- Subhalos
- Baryons
- Miscentering
- Full survey covariance
- HOD modeling

This is an inner-halo-focused ensemble forecast.

---

# SIDM Parameterization

We adopt elastic SIDM with effective cross section:

σ/m_eff  (cm²/g)

Cluster benchmark grid:

0.2, 0.5, 1.0, 2.0 cm²/g

---

# Pipeline Architecture

Modules:

config.py
cosmology.py
profiles.py
sidm_wrapper.py
projection.py
ensemble.py
stacking.py
forecast.py

Flow:

1. Sample halo ensemble
2. Generate CDM profile
3. Generate SIDM profile
4. Project to ΔΣ(R)
5. Interpolate to common R grid
6. Stack
7. Compute Δχ²

---

# Ensemble Construction (Tier-1)

Minimum implementation:

- Mass distribution: log-normal
- Mean mass: 3×10^14 Msun (cluster regime)
- Scatter: 0.2 dex
- Fixed redshift (default z=0.4)
- Concentration from mass–concentration relation
- Equal-weight stacking

---

# Outputs

- Stacked ΔΣ(R) curves
- Ratio ΔΣ_SIDM / ΔΣ_CDM
- Δχ² vs σ/m
- Summary precision table

---

# Validation

Must verify:

- Projection matches analytic NFW ΔΣ
- Single-halo limit recovered
- Mass conservation
- Stack convergence with N_halos

---

# Limitations

This Tier-1 model captures:

- Inner halo SIDM modifications
- Ensemble stacking effects

It does NOT capture:

- Self-consistent splashback physics
- Accretion-region dynamics
- Subhalo population changes

These require simulation-calibrated modeling.

---

Tier-1 = Ensemble stacking without splashback.
Future tiers may extend outskirts modeling.
