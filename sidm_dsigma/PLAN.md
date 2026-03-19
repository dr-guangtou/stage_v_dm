# PLAN.md

## Development roadmap

### Phase 0 — bootstrap
- Create package skeleton and configuration file.
- Define cosmology, mass definition, and benchmark cases.
- Add a reproducible environment file.

### Phase 1 — baseline profiles
- Implement NFW baseline from `M200`, `c200`, and `z`.
- Wrap `parametricSIDM` behind a stable local function.
- Verify basic outputs: `rho(r)`, `M(<r)`, `Vc(r)`.

### Phase 2 — projection
- Implement numerical `Sigma(R)`.
- Implement `Sigma_bar(<R)` and `DeltaSigma(R)`.
- Validate projection on a standard NFW case against analytic or trusted-reference expectations.

### Phase 3 — benchmark forecast
- Run dwarf benchmark grid.
- Run cluster benchmark grid.
- Plot `rho(r)`, `DeltaSigma(R)`, and ratios relative to CDM.
- Add toy fractional error models and compute `Delta chi^2`.

### Phase 4 — outputs
- Produce one notebook for end-to-end reproduction.
- Save figures to `outputs/figures/`.
- Save summary tables to `outputs/tables/`.
- Write a short interpretation note in the notebook output.

## Suggested implementation order

1. `config.py`
2. `cosmology.py`
3. `profiles.py`
4. `projection.py`
5. `forecast.py`
6. `plotting.py`
7. notebook / CLI script

## Near-term stretch goals

- Optional baryonic placeholder term.
- Optional reference cross-check with `pyccl` or `colossus` for standard profiles.
- Optional precision sweep to report the required measurement accuracy for `2sigma` / `3sigma` separation.

## Explicit non-goals for v1

- full HOD,
- 2-halo term,
- miscentering,
- realistic covariance,
- photo-z systematics,
- central/satellite population modeling,
- publication-grade survey forecasting.

# ============================================
# EXTENSION: Tier-1 Halo Ensemble Development
# ============================================

Assume single-halo pipeline complete.

---

## Phase 6 — Halo Ensemble Module

- Implement halo sampler
- Add configurable mass distribution
- Add concentration scatter
- Add reproducible seed support
- Validate ensemble statistics

---

## Phase 7 — Stacking Module

- Implement interpolation to common R grid
- Implement weighted averaging
- Validate single-halo limit
- Verify convergence with increasing N_halos

---

## Phase 8 — Ensemble Forecast

- Generate stacked ΔΣ(R) for:
    CDM
    σ/m = 0.2, 0.5, 1.0, 2.0
- Produce ratio plots
- Compute Δχ²

---

## Phase 9 — Stability Tests

- Increase halo count
- Test mass scatter sensitivity
- Confirm SIDM signal survives stacking

---

## Milestone

Produce a 4-panel summary figure:

1. Ensemble ρ(r)
2. Ensemble ΔΣ(R)
3. Ratio ΔΣ_SIDM / ΔΣ_CDM
4. Δχ² vs σ/m

Tier-1 complete when:

- Ensemble stacking works
- Results reproducible
- ΔΣ differences interpretable
