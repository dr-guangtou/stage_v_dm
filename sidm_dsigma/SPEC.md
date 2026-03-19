# SPEC.md

## Project title

Minimal Stage-V SIDM Lensing Forecast

## Objective

Implement a lightweight forecast pipeline that uses a parametric SIDM halo model to compare CDM and SIDM predictions for:

- the 3D density profile `rho(r)`, and
- the galaxy-galaxy / cluster-galaxy lensing observable `DeltaSigma(R)`

for two benchmark halo masses:

- `M200 = 1e10 Msun`
- `M200 = 1e14 Msun`

The pipeline is intended for a 1st-order science-case forecast that can be completed quickly and handed off to a coding agent.

---

## 1. Scientific assumptions

### 1.1 Forecast philosophy

This is not a full survey forecast. It is a controlled, interpretable comparison of halo-profile predictions.

We want to answer:

> For a specified SIDM parameterization, how different are the predicted `rho(r)` and `DeltaSigma(R)` from CDM, and what approximate measurement precision would be needed to tell them apart?

### 1.2 SIDM parameterization

Use a halo-scale effective cross section:

```text
sigma_over_m_eff   [cm^2 / g]
```

This is the only required SIDM parameter in version 1.

### 1.3 Benchmark grids

#### Dwarf case
```text
M200 = 1e10 Msun
c200 = 15.0   (default; allow override)
z = 0.3       (default; allow override)
sigma/m = [0.0, 0.2, 0.5, 1.0, 2.0]
```

#### Cluster case
```text
M200 = 1e14 Msun
c200 = 4.0    (default; allow override)
z = 0.3       (default; allow override)
sigma/m = [0.0, 0.2, 0.5, 1.0, 2.0]
```

### 1.4 Version-1 exclusions

The following are explicitly out of scope unless already easy to add:

- halo occupation distribution,
- halo mass scatter,
- satellite fraction modeling,
- miscentering,
- 2-halo term,
- intrinsic alignment,
- source redshift uncertainty propagation,
- realistic survey covariance,
- baryonic components beyond an optional placeholder.

---

## 2. External software guidance

### 2.1 Core dependency

Use:

- `parametricSIDM` as the SIDM profile engine.

### 2.2 Projection guidance

Do not use `cluster_toolkit` directly.

Implement the projection integrals in-house with `numpy` and `scipy`.

Allowed supporting packages:

- `numpy`
- `scipy`
- `matplotlib`
- `astropy`
- `pyccl` (optional reference / validation)
- `colossus` (optional reference / validation)
- `halomod` (optional future extension; not needed for v1)

### 2.3 Preferred use of other libraries

- Use CCL as a reference for conventions and, if useful, cross-checks of projected quantities for standard profiles.
- Use Colossus as a reference for halo-profile conventions and optional consistency tests against standard NFW `Sigma` and `DeltaSigma`.
- Do not make any of these packages mandatory for the main projection path.

---

## 3. Required modules

### 3.1 `config.py`

Purpose:
- centralize constants and benchmark settings.

Required contents:
- cosmological parameters (simple flat LCDM defaults),
- benchmark masses,
- benchmark concentrations,
- default redshift,
- radial grids,
- SIDM cross-section grids,
- toy error models.

### 3.2 `cosmology.py`

Purpose:
- provide basic cosmology helpers.

Required functions:

```python
def rho_crit_z(z: float, cosmo: dict) -> float:
    ...

def rdelta(M: float, z: float, cosmo: dict, definition: str = "200c") -> float:
    ...
```

Notes:
- choose one mass definition and use it consistently everywhere;
- default to `M200c` unless there is a compelling reason not to.

### 3.3 `profiles.py`

Purpose:
- generate CDM and SIDM 3D profiles.

Required functions:

```python
def nfw_profile_from_m_c(
    r: np.ndarray,
    M200: float,
    c200: float,
    z: float,
    cosmo: dict,
) -> np.ndarray:
    ...


def sidm_profile_from_parametric_model(
    r: np.ndarray,
    M200: float,
    c200: float,
    z: float,
    sigma_over_m: float,
    model_options: dict | None = None,
) -> dict:
    ...
```

Expected return dictionary keys for SIDM function:

```python
{
    "r": r,
    "rho": rho,
    "m_enclosed": m_enclosed,
    "vcirc": vcirc,
    "metadata": {...}
}
```

Requirements:
- wrap `parametricSIDM` in a stable local interface;
- keep the wrapper thin and transparent;
- gracefully fail with an informative message if the upstream code changes.

### 3.4 `projection.py`

Purpose:
- compute lensing profiles from 3D density profiles.

Required functions:

```python
def sigma_of_R(
    R: np.ndarray,
    r: np.ndarray,
    rho: np.ndarray,
    r_max: float | None = None,
) -> np.ndarray:
    ...


def sigma_bar_of_R(
    R: np.ndarray,
    sigma_R: np.ndarray,
) -> np.ndarray:
    ...


def delta_sigma_of_R(
    R: np.ndarray,
    r: np.ndarray,
    rho: np.ndarray,
) -> dict:
    ...
```

Definitions:

```text
Sigma(R)      = 2 * integral_R^infinity dr [ rho(r) * r / sqrt(r^2 - R^2) ]
Sigma_bar(R)  = (2 / R^2) * integral_0^R dR' [ R' * Sigma(R') ]
DeltaSigma(R) = Sigma_bar(R) - Sigma(R)
```

Implementation requirements:
- use robust interpolation of `rho(r)` on a logarithmic radial grid;
- avoid numerical instability near `r = R`;
- expose integration controls;
- document whether outputs are in physical `Msun / kpc^2` or another unit;
- include at least one unit test against an analytic NFW or a trusted numerical reference.

### 3.5 `forecast.py`

Purpose:
- compute distinguishability metrics.

Required functions:

```python
def fractional_error_model(
    R: np.ndarray,
    regime: str,
    scenario: str = "baseline",
) -> np.ndarray:
    ...


def delta_chi2(
    model: np.ndarray,
    reference: np.ndarray,
    sigma: np.ndarray,
) -> float:
    ...


def required_fractional_precision(
    model: np.ndarray,
    reference: np.ndarray,
    target_snr: float,
) -> float:
    ...
```

Interpretation:
- `delta_chi2` can assume diagonal covariance in v1;
- `required_fractional_precision` should return a simple scalar summary for slide-ready communication.

### 3.6 `plotting.py`

Purpose:
- generate publication-quality plots.

Required plots:

1. `rho(r)` for CDM and SIDM benchmarks, dwarf and cluster
2. `DeltaSigma(R)` for CDM and SIDM benchmarks, dwarf and cluster
3. ratio plots:
   - `rho_SIDM / rho_CDM`
   - `DeltaSigma_SIDM / DeltaSigma_CDM`
4. optional residual plots in units of assumed observational error

Style requirements:
- readable on workshop slides,
- log-log axes where appropriate,
- clear legends with halo mass and `sigma/m`,
- no hidden smoothing.

### 3.7 `io.py`

Purpose:
- save outputs.

Required outputs:
- CSV or ECSV tables for profiles,
- PNG and PDF figures,
- a compact JSON summary of benchmark distinguishability metrics.

---

## 4. Radial grids and units

### 4.1 Dwarf grids

3D grid:

```text
r = logspace(-1, 2.5) kpc   # 0.1 to ~316 kpc
```

Projected grid:

```text
R = logspace(0.5, 2.5) kpc  # ~3 to ~316 kpc
```

### 4.2 Cluster grids

3D grid:

```text
r = logspace(0.7, 3.7) kpc  # ~5 to ~5000 kpc
```

Projected grid:

```text
R = logspace(1.5, 3.5) kpc  # ~30 to ~3160 kpc
```

### 4.3 Units

Use:

- `r`, `R` in physical kpc
- `rho` in `Msun / kpc^3`
- `Sigma`, `DeltaSigma` in `Msun / kpc^2`
- `M200` in `Msun`

All unit conversions must be explicit.

---

## 5. Toy observational model

### 5.1 Dwarf case

Implement at least two scenarios:

#### baseline
```text
fractional error per bin = 0.15
```

#### conservative
```text
fractional error per bin = 0.30
```

Allow a third optional optimistic scenario:

```text
fractional error per bin = 0.10
```

### 5.2 Cluster case

Implement at least two scenarios:

#### baseline
```text
fractional error per bin = 0.05
```

#### conservative
```text
fractional error per bin = 0.10
```

Optional optimistic scenario:

```text
fractional error per bin = 0.03
```

These are placeholders for a science-case forecast, not final survey requirements.

---

## 6. Required outputs

### 6.1 Figures

Produce at least the following:

#### Figure 1: 4-panel core comparison
Panels:
1. dwarf `rho(r)`
2. dwarf `DeltaSigma(R)`
3. cluster `rho(r)`
4. cluster `DeltaSigma(R)`

#### Figure 2: ratio figure
Panels:
1. dwarf `rho_SIDM / rho_CDM`
2. dwarf `DeltaSigma_SIDM / DeltaSigma_CDM`
3. cluster `rho_SIDM / rho_CDM`
4. cluster `DeltaSigma_SIDM / DeltaSigma_CDM`

#### Figure 3: detectability summary
For each regime, show either:
- `sqrt(Delta chi^2)` vs `sigma/m`, or
- required fractional precision vs `sigma/m`.

### 6.2 Tables

#### Table A: benchmark setup
Columns:
- regime
- `M200`
- `c200`
- `z`
- `sigma_over_m`
- radial range used

#### Table B: detectability summary
Columns:
- regime
- `sigma_over_m`
- max fractional change in `DeltaSigma`
- radial range of strongest leverage
- `Delta chi^2` for each error scenario
- required fractional precision for `2 sigma` and `3 sigma`

---

## 7. Numerical validation requirements

At minimum, perform the following checks:

1. CDM sanity check
   - verify that the NFW profile integrates to the intended halo mass definition within tolerance.

2. Projection sanity check
   - compare the numerical `Sigma(R)` / `DeltaSigma(R)` for an NFW halo against an analytic or trusted reference calculation.

3. Monotonicity / positivity checks
   - no negative `rho`, `Sigma`, or `DeltaSigma` in the benchmark outputs.

4. Resolution test
   - verify that doubling radial resolution changes `DeltaSigma` negligibly over the quoted radial range.

---

## 8. Coding style and reproducibility

Requirements:

- Python 3.11+
- type hints where practical
- concise docstrings
- deterministic benchmark outputs
- all benchmark choices configurable from a single file
- avoid notebook-only logic; core computation must live in importable modules

Recommended:

- `ruff` or `black` formatting
- small unit tests for projection routines

---

## 9. Nice-to-have extensions (not required for v1)

Only attempt after the core benchmark works:

1. optional Hernquist stellar component
2. optional 2-halo term using CCL / halomod-inspired machinery
3. redshift dependence (`z = 0.1, 0.3, 0.6`)
4. concentration variation sensitivity
5. mapping between `(sigma/m)_eff` and a velocity-dependent particle model
6. optional use of the dedicated `SIDM_Lensing_Model` repository as a cross-check

---

## 10. Acceptance criteria

The project is successful when:

1. the code runs end-to-end for the dwarf and cluster benchmarks;
2. it generates `rho(r)` and `DeltaSigma(R)` for CDM and SIDM cases;
3. it outputs the required figures and summary tables;
4. it provides at least one slide-ready statement about the measurement precision needed to distinguish CDM from SIDM in each regime.

---

## 11. References to consult during implementation

- `parametricSIDM` GitHub repository for the core SIDM profile workflow.
- `SIDM_Lensing_Model` GitHub repository for lensing-specific SIDM calculations and possible cross-checks.
- CCL documentation for halo-profile projection conventions.
- Colossus documentation for standard profile and `DeltaSigma` reference behavior.
- CLMM documentation / paper for weak-lensing definitions and unit conventions.
