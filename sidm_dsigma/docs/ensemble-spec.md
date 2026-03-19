# ENSEMBLE SPECIFICATION
## Flexible Halo Ensemble Framework
### For SIDM vs CDM Stacked ΔΣ(R) Forecasts

This document defines two ensemble-generation frameworks:

1. HMF-based (halo-mass-selected)
2. SHMR-based (galaxy-selected)

Both are selection-aware and configurable.

These frameworks feed into:

Halo Ensemble → ρ(r) → Σ(R) → ΔΣ(R) → Stack → Forecast

---

# Design Philosophy

Simulation-based SIDM constraints (e.g., cluster weak lensing) rely on:

- Realistic halo mass distributions
- Selection effects
- Ensemble stacking

Therefore, our ensemble framework must:

- Be flexible
- Explicitly encode selection
- Be survey-configurable
- Remain computationally lightweight

---

# Framework 1: HMF-Based Ensemble
## (Halo-Mass-Selected)

## Scientific Use Case

- SZ-selected clusters
- Richness-selected clusters
- Mass-limited samples
- Cosmology-driven halo selection

Typical regime:
M200 ≳ 1e14 Msun

---

## 1. Mathematical Definition

The ensemble is drawn from:

P(M, z) ∝ S(M, z) × dn/dM(M, z)

Where:

dn/dM = halo mass function
S(M, z) = selection function

---

## 2. Minimal Implementation (Tier-1)

Assumptions:

- Fixed redshift z = z0
- Mass range: [M_min, M_max]
- Selection = sharp mass threshold

P(M) ∝ dn/dM(M, z0),  M_min < M < M_max

Default cluster configuration:

z0 = 0.4
M_min = 1e14 Msun
M_max = 1e15 Msun

---

## 3. Advanced Implementation (Selection-Aware)

Allow:

S(M, z) = logistic or error-function selection

Example:

S(M) = 1 / [1 + exp(-(logM - logM_cut)/σ_sel)]

This mimics:

- Richness selection
- SZ selection threshold
- Observable–mass scatter

---

## 4. Redshift Sampling Options

Option A (Tier-1):
- Fixed z

Option B:
- Sample z from survey distribution
- Then sample M from dn/dM(M, z)

---

## 5. Concentration Model

c(M, z) from mass–concentration relation

Add lognormal scatter:

σ_logc ≈ 0.15–0.2

---

## 6. Weights for Stacking

Allow:

- Equal weight
- Weight ∝ 1/Σ_crit²(z)
- Weight ∝ observable proxy

Default Tier-1: equal weight

---

## 7. Summary: HMF-Based Ensemble

Best for:

- Cluster samples
- Halo-mass-selected analyses
- Direct comparison to simulation literature

Not appropriate for dwarf galaxy lens samples.

---

# Framework 2: SHMR-Based Ensemble
## (Galaxy-Selected)

## Scientific Use Case

- Spectroscopic dwarf galaxies
- Stellar-mass-selected samples
- Galaxy-targeted surveys
- Field galaxy lensing

---

## 1. Mathematical Definition

We draw galaxies from:

P(M_star, z) = galaxy stellar mass function

Then map to halo mass via SHMR:

P(M_halo | M_star)

Total halo probability:

P(M_halo, z) =
∫ P(M_star, z) × P(M_halo | M_star) dM_star

---

## 2. Minimal Implementation (Tier-1)

Assumptions:

- Fixed z = z0
- Stellar mass log-normal distribution
- Deterministic SHMR + scatter

Example dwarf configuration:

z0 = 0.2
log10(M_star/Msun) ~ Normal(7.5, 0.5 dex)
Halo scatter at fixed M_star = 0.3 dex

---

## 3. SHMR Model Options

Option A (Simple Power Law):

log M_halo = A + B log M_star

Option B (Behroozi/Moster parameterization)

Option C:
Tabulated SHMR with interpolation

---

## 4. Central vs Satellite Handling

Tier-1:

Assume all galaxies are centrals.

Advanced option:

Include satellite fraction f_sat
Sample subhalo masses separately.

---

## 5. Concentration

Same as HMF-based:

c(M, z) with scatter.

---

## 6. Weights for Stacking

Options:

- Equal per galaxy
- Σ_crit⁻² weighting
- Luminosity weighting

Default Tier-1: equal per galaxy

---

## 7. Summary: SHMR-Based Ensemble

Best for:

- Dwarf galaxy lensing
- Stellar-mass-selected surveys
- Velocity-scale comparisons

Not appropriate for pure halo-mass-selected cluster samples.

---

# Comparison Table

| Feature | HMF-Based | SHMR-Based |
|----------|------------|-------------|
| Selection | Halo mass | Stellar mass |
| Typical mass | ≥1e14 Msun | 1e9–1e11 Msun |
| Survey analogue | SZ / richness | Spectroscopic galaxies |
| Outer profile realism | Moderate | Low |
| Velocity regime | High | Low |

---

# API Design

Both frameworks must implement:

generate_ensemble(
    mode="HMF" or "SHMR",
    config_dict
)

Returns:

List of halos:

{
  "M200": float,
  "z": float,
  "c200": float,
  "weight": float
}

---

# Explicit Non-Goals (Tier-1)

- Splashback modeling
- Subhalo population modeling
- Environmental dependence
- Assembly bias

---

# Future Extensions

- DK14 outer profile attachment
- Velocity-dependent σ/m mapping
- Observable-weighted selection
- Real survey redshift distributions

---

End of Specification.
