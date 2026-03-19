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
