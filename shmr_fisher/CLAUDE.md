# CLAUDE.md — SHMR Forecast Project

## Project Overview

This is a **Fisher forecast tool for the stellar-halo mass relation (SHMR)**, predicting how well spectroscopic surveys (combined with Stage-IV imaging lensing) can constrain the galaxy-halo connection. The user is a professional astronomer (associate professor) who codes primarily in Python. The project is scientific research code, not a web application.

## Coding Standards

### Language & Style
- **Python 3.10+**. Use type hints throughout. Use `dataclasses` for configuration objects.
- Follow **PEP 8**. Use `snake_case` for functions and variables, `CamelCase` for classes.
- Maximum line length: 100 characters (not 79 — scientific code often has long expressions).
- Use f-strings for string formatting.

### Documentation Is Mandatory
- Every module: a docstring at the top explaining its purpose and relationship to the forecast pipeline.
- Every public function: a numpy-style docstring with `Parameters`, `Returns`, and `Notes` sections.
- Every class: docstring explaining the physical meaning of each attribute, including **units**.
- Inline comments for any non-obvious physics or numerical choices (e.g., why a particular integration range, why a specific step size for derivatives).
- When implementing an equation from a paper, cite it: `# Moster+2013 Eq. 2` or `# Following Leauthaud+2012 Eq. 15`.

### Units — This Is Critical
- **Always document units in variable names or comments.** Ambiguous units are the #1 source of bugs in astrophysics code.
- Halo masses: `Msun` (NOT `Msun/h`). If colossus returns h-scaled quantities, convert immediately and document.
- Distances for ΔΣ radial bins: physical Mpc.
- ΔΣ: compute in `Msun/Mpc^2` internally; convert to `Msun/pc^2` only for plotting.
- Colossus profiles use comoving kpc/h internally. Every call to colossus profile methods must have a comment explaining the unit conversion.
- Angular diameter distances: physical Mpc.
- Survey areas: deg² in configs, convert to steradians for calculations.
- Number densities: per Mpc³ (comoving) for volume densities, per deg² for survey densities.

### Numerical Code
- Use `np.logspace` for halo mass grids (they span 5+ orders of magnitude).
- Use `scipy.integrate.quad` or `np.trapezoid` on log-spaced grids for HMF integrals. Prefer `np.trapezoid` on dense grids for speed; use `quad` only when accuracy matters more than speed.
- Finite differences: central differences with configurable step size. Always check convergence.
- When inverting the Fisher matrix, check the condition number first. Log a warning if > 1e10.

### Error Handling
- Raise `ValueError` with a descriptive message for physically nonsensical inputs (e.g., negative halo mass, z < 0, empty stellar mass bins).
- Use `warnings.warn` for non-fatal issues (e.g., Fisher matrix has high condition number, a stellar mass bin has < 100 predicted galaxies).
- Never silently return NaN or zero. If a calculation fails, raise or warn.

### Testing
- Write validation functions (not full pytest yet) that check:
  - SHMR peak efficiency is 3–5% at Mh ~ 10^12 at z=0
  - ΔΣ for a single NFW halo matches expected order of magnitude
  - Fisher matrix is positive definite
  - Doubling N_gal halves errors (in shape-noise-dominated regime)
- Put these in a `validate.py` module. Run them after each phase.

## Key Dependencies

| Package | Purpose | Notes |
|---------|---------|-------|
| `colossus` | Cosmology, HMF, halo bias, NFW profiles, concentration | Primary engine. Use `cosmology.setCosmology('planck18')` at startup. |
| `numpy` | Arrays, integration | Standard |
| `scipy` | Root-finding (brentq), special functions (erf), integration (quad) | For SHMR inverse and HOD |
| `matplotlib` | Plotting | For validation and science figures |
| `astropy` | Units, constants (G, c) | For Σ_crit calculation. Use `astropy.constants` and `astropy.units`. |
| `pyyaml` | YAML config parsing | For loading run configurations from YAML files |

**Do NOT use** `halotools` or `halomod` — they require simulation data downloads and are overkill for this analytic forecast.

## Project Structure

```
shmr_fisher/
├── CLAUDE.md              # This file
├── SPEC.md                # Detailed specification
├── DESIGN.md              # Pipeline design document (how the forecast works)
├── configs/               # YAML run configurations
│   ├── default.yaml       # 4-survey comparison (Stage-III through Stage-V)
│   ├── stage4_vs_stage5.yaml  # 2-survey head-to-head
│   └── dwarf_regime.yaml  # Single low-z survey for dwarf science
├── shmr_fisher/           # Main package
│   ├── __init__.py
│   ├── config.py          # Dataclasses: SHMRParams, SurveyConfig, LensingConfig, ForecastConfig, NuisanceConfig
│   ├── config_io.py       # YAML config loader: RunConfig, load_run_config()
│   ├── shmr_model.py      # SHMR parameterization (Moster+2013 with z-evolution)
│   ├── halo_model.py      # ΔΣ (1-halo NFW + 2-halo), b_eff, n_gal via HMF integrals
│   ├── covariance.py      # Lensing + clustering covariance matrices
│   ├── fisher.py          # Fisher matrix computation, marginalization, priors
│   ├── systematics.py     # Systematic error helpers (floor, nuisance derivatives)
│   ├── survey_configs.py  # Predefined Stage-III/IV/V survey configurations (legacy)
│   ├── validate.py        # Validation/sanity checks after each phase
│   └── plot_results.py    # All plotting functions (dynamic survey colors)
├── scripts/
│   ├── run_forecast.py    # Main driver: --config YAML or legacy CLI
│   ├── run_sweep.py       # Parameter sweep driver: --config or --base
│   └── generate_science_figures.py  # Full figure generation from YAML config
├── outputs/               # Generated outputs, organized by run name
│   ├── {run_name}/        # Per-run directory (e.g., outputs/default/)
│   │   ├── config.yaml    # Copy of input config for reproducibility
│   │   ├── forecast_results.json
│   │   ├── forecast_results.npz
│   │   ├── CAPTION.md     # Figure captions for this run
│   │   └── *.png          # Figures (300 dpi PNG)
│   └── phase1/...phase4/  # Legacy phase-based outputs (from development)
├── docs/                  # Development journals
└── README.md              # User-facing documentation
```

## Domain Context for Claude Code

### What is the SHMR?
The stellar-halo mass relation maps dark matter halo mass (Mh) to the stellar mass (M*) of the galaxy it hosts. It peaks at M*/Mh ~ 3–5% around Mh ~ 10^12 Msun and falls off as a double power-law at lower and higher masses. The scatter is log-normal at fixed Mh: either constant (σ_logMs ~ 0.15 dex) or mass-dependent following Cao & Tinker (2020), increasing from ~0.18 dex at cluster masses to ~0.38 dex at dwarf masses.

### What are the observables?
- **Galaxy-galaxy lensing ΔΣ(R):** The excess surface mass density around galaxies, measured from weak gravitational lensing. Directly probes the mean halo mass of a galaxy sample.
- **Galaxy clustering w_p(r_p):** The projected two-point correlation function. Probes halo mass via the halo bias. In this project, we simplify to (b_eff, n_gal) per bin.

### What is the Fisher matrix doing?
Given a parameterized SHMR model and survey noise properties, the Fisher matrix tells us how precisely each parameter can be measured. It is the expected curvature of the log-likelihood at the fiducial model. The inverse of the Fisher matrix gives the parameter covariance.

### Two science regimes
- **Dwarf galaxies (M* < 10^9.5, z < 0.4):** Redshift evolution is weak. The science goal is to constrain the low-mass slope (beta) and scatter of the z~0 SHMR.
- **Massive galaxies at higher z (M* > 10^10.5, z up to 1.0):** Redshift evolution matters. The science goal is to constrain how the SHMR evolves (the nu_* parameters).

## Workflow Conventions

### Phase-by-phase development
This project is developed incrementally following the plan in `SPEC.md`. After completing each phase:
1. Run the relevant validation checks in `validate.py`.
2. Generate any diagnostic plots specified in the plan.
3. Commit with a message like `Phase 2: halo model predictions — ΔΣ validated`.

### When stuck on colossus units
The most likely failure mode is unit conversion errors with colossus NFW profiles. If ΔΣ values are off by orders of magnitude:
1. Check whether colossus is returning comoving vs physical quantities.
2. Check h-scaling: colossus often works in Msun/h and kpc/h internally.
3. Write a standalone test: compute ΔΣ for M200m = 10^14 Msun, c=5, z=0. The result should be ~100 Msun/pc^2 at R ~ 0.5 Mpc. If you get ~10^8 or ~10^{-4}, you have a unit error.
4. As a fallback, implement the analytic Wright & Brainerd (2000) NFW ΔΣ directly — it's a well-known closed-form expression.

### When the Fisher matrix is ill-conditioned
If `np.linalg.cond(fisher) > 1e10`:
1. Check which parameters are poorly constrained (look at diagonal elements).
2. Most likely cause: evolution parameters (nu_*) are unconstrained by a low-z-only survey. Fix them and re-run with 5 parameters.
3. If the 5-parameter Fisher is still ill-conditioned, check for parameter degeneracies (e.g., log_M1 and N are correlated). Try reparameterizing.

### Output conventions
- Save Fisher matrices and errors as JSON (for portability) and as numpy `.npz` files (for reloading).
- **Default figure format: PNG** at 300 dpi, with descriptive filenames like `shmr_validation.png`, `fisher_comparison_stage4_vs_stage5.png`.
- All figures should have axis labels with units, legends, and a title or annotation indicating the survey configuration.
- **Per-run output directories:** When using YAML configs, all outputs go to `outputs/{run_name}/` where `run_name` is the config filename stem (e.g., `default.yaml` → `outputs/default/`). The YAML config is copied into the output directory for reproducibility.

### YAML-driven workflow (preferred)
- Survey configurations are defined in YAML files under `configs/`.
- Each config can define 1–4 surveys. The run name is derived from the filename.
- Use `--config configs/default.yaml` with `run_forecast.py`, `run_sweep.py`, or `generate_science_figures.py`.
- Legacy CLI flags (`--surveys`, `--systematics`) still work for backward compatibility.

### Figure captioning (MANDATORY)
- Every time a figure is saved to `outputs/`, a detailed caption **must** be appended to the run's `CAPTION.md` (either `outputs/{run_name}/CAPTION.md` or the legacy `outputs/CAPTION.md`).
- Each caption entry must include: the filename, creation timestamp, a full description of what is plotted (axes, units, data shown, colors/line styles), the physical interpretation of key features, and the figure's purpose in the pipeline.
- This ensures figures remain interpretable months later without re-running the code.
