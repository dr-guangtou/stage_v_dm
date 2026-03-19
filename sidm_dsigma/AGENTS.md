# AGENTS.md

## Purpose

This repository builds a minimal, first-order forecast for distinguishing CDM from SIDM using projected lensing signals derived from parametric SIDM halo density profiles.

The scope is intentionally narrow:
- benchmark halos only (`1e10` and `1e14 Msun`),
- benchmark effective SIDM cross sections,
- direct `rho(r) -> Sigma(R) -> DeltaSigma(R)` projection,
- toy error model and `Delta chi^2` distinguishability.

## Core principles

1. Keep version 1 simple, transparent, and reproducible.
2. Prefer explicit units and documented conventions.
3. Avoid hidden dependencies and magic defaults.
4. Do not add survey realism unless it is explicitly requested.
5. Optimize for interpretability, not sophistication.

## Agent roles

### 1. Architecture / integration agent
Owns project structure, interfaces, configuration, and CLI/notebook wiring.

### 2. Physics / profiles agent
Owns CDM baseline profiles, SIDM wrapper logic, mass-definition consistency, and profile sanity checks.

### 3. Projection / numerics agent
Owns `Sigma(R)` and `DeltaSigma(R)` calculations, interpolation, numerical integration stability, and validation.

### 4. Forecast / plotting agent
Owns toy error models, `Delta chi^2`, summary tables, figure generation, and final notebook outputs.

## Guardrails

- Use `parametricSIDM` as the SIDM profile engine.
- Do not use `cluster_toolkit` as a required dependency.
- Projection should be implemented locally with `numpy`/`scipy`.
- `pyccl` and `colossus` are allowed as references or validation tools, not required runtime dependencies.
- Keep all outputs in well-labeled physical units unless there is a strong reason not to.

## Minimum acceptance bar

A successful first version should:
- run from a clean environment,
- generate SIDM and CDM 3D profiles for dwarf and cluster benchmarks,
- project them to `DeltaSigma(R)`,
- compute a simple distinguishability metric,
- save at least one slide-ready summary figure and one machine-readable summary table.

## When in doubt

Choose the simpler implementation that preserves scientific clarity.
