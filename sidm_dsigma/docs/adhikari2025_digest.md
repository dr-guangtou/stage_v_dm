# susmita2025_digest

## Scope

This note summarizes two cluster-scale SIDM papers and assesses how their modeling relates to the `parametricSIDM` framework for quick weak-lensing forecasts:

- Banerjee et al. 2019 / JCAP 2020, *Signatures of Self-Interacting dark matter on cluster density profile and subhalo distributions*.
- Adhikari et al. 2024, *Constraints on Dark Matter Self-interactions from Weak Lensing of Galaxies from the Dark Energy Survey around Clusters from the Atacama Cosmology Telescope Survey*.

It also compares them to the Yang et al. `parametricSIDM` approach.

---

## Executive summary

- The two papers are **consistent in spirit** with `parametricSIDM`, but **not identical in implementation**.
- For **elastic SIDM**, `parametricSIDM` is broadly compatible with the cluster-lensing studies at the level of **spherically averaged halo profiles**.
- For **dissipative SIDM**, `parametricSIDM` is **not** a drop-in replacement, because Adhikari et al. use a different fluid-cooling framework.
- A first-order reproduction of their **stacked** \(\Delta\Sigma(R)\) trends is feasible for the **elastic** models if one generates a realistic halo ensemble and projects the 3-D profiles to lensing space.
- Exact reproduction of the published \(\Delta\Sigma\) curves requires the authors’ simulation/stacking details, not just a single fiducial halo.

---

## 1. How the two papers implement SIDM

### 1.1 Banerjee et al. 2019 / JCAP 2020

This work uses **cosmological simulations with SIDM** to study how cluster observables respond to self-interactions. The focus is on:

- stacked cluster density profiles,
- subhalo counts and distributions,
- splashback-scale signatures.

A key point is that the scattering kernel can depend on **relative velocity** and **scattering angle**. The paper therefore goes beyond the simplest constant, isotropic elastic SIDM setup.

Implication for forecasting:

- this is a **simulation-based cluster-signature study**,
- not a simple analytic or single-halo profile exercise,
- and not primarily a direct weak-lensing likelihood paper.

### 1.2 Adhikari et al. 2024

This paper directly compares SIDM predictions to **observed weak-lensing** around ACT clusters measured with DES background galaxies.

Modeling choices:

- **elastic, isotropic SIDM**: modeled with **full cosmological dark-matter-only N-body simulations**;
- **dissipative SIDM**: modeled with **semi-analytic fluid simulations** including an explicit energy-loss term.

This paper is therefore much closer to a direct \(\Delta\Sigma(R)\)-based constraint analysis.

---

## 2. Relation to `parametricSIDM`

### 2.1 Where the approaches are compatible

Yang et al.'s `parametricSIDM` model is a **fast surrogate** for SIDM gravothermal evolution. It maps a CDM halo to an SIDM halo using a calibrated analytical profile and can represent more general scattering models through an **effective constant cross section** when one only cares about the halo-scale, spherically averaged structural response.

That means:

- for **elastic, isotropic, constant-\(\sigma/m\)** SIDM, `parametricSIDM` is conceptually well aligned;
- for **velocity-/angle-dependent elastic SIDM**, `parametricSIDM` can still be useful if these are compressed into an **effective cross section** at the halo level.

### 2.2 Where the approaches differ

There are two important mismatches.

1. **Single-halo surrogate vs stacked simulation ensemble**  
   `parametricSIDM` predicts the spherical density profile of an SIDM halo given a CDM counterpart. The Banerjee and Adhikari papers predict observables from **ensembles of simulated clusters**, which implicitly include halo-to-halo scatter, selection, and stacking.

2. **Elastic gravothermal evolution vs dissipative cooling**  
   Adhikari et al.'s dissipative models include an explicit **cooling/energy-loss** prescription. This is not the same physical model as the public `parametricSIDM` implementation.

### 2.3 Practical verdict

- **Elastic SIDM in Adhikari et al. 2024:** broadly reproducible in a first-order sense with `parametricSIDM` + projection + stacking.
- **Velocity-/angle-dependent elastic SIDM in Banerjee et al. 2019:** approximately reproducible only after mapping to an effective halo-level cross section.
- **Dissipative SIDM in Adhikari et al. 2024:** not reproducible with `parametricSIDM` alone.

---

## 3. Can we reproduce their \(\Delta\Sigma(R)\) profiles?

### 3.1 Elastic SIDM: yes, approximately

A feasible workflow is:

\[
\{M, c, z, \sigma/m\}_{\rm halo\ ensemble}
\rightarrow \rho_{\rm SIDM}(r)
\rightarrow \Sigma(R),\Delta\Sigma(R)
\rightarrow \text{stacked prediction}.
\]

This should recover the **qualitative** and likely **semi-quantitative** behavior of the elastic SIDM cluster weak-lensing profiles.

However, exact agreement with the published curves would require:

- the cluster mass and redshift distribution used in the analysis,
- halo concentration scatter,
- the precise stacking and sample selection,
- and the outer-halo / splashback treatment.

### 3.2 Dissipative SIDM: no, not directly

The dissipative models in Adhikari et al. are not just another elastic cross section. They involve an explicit **energy-loss scale** and are generated with a different fluid framework. A separate dissipative module would be needed.

### 3.3 Banerjee et al. 2019 nuance

For the 2019 paper, the more realistic target is to reproduce:

- the **3-D stacked density trends**,
- and possibly corresponding projected surface-density trends,

rather than every feature of the subhalo and splashback statistics.

---

## 4. SIDM parameter spaces explored

### 4.1 Banerjee et al. 2019

The paper studies **elastic SIDM** with emphasis on **velocity and angular dependence** in the differential cross section.

From the paper abstract and associated presentation material, the explored benchmarks include:

- velocity-independent / isotropic-style comparisons around
  - \(\sigma_T/m \sim 1\ \mathrm{cm^2\,g^{-1}}\),
  - \(\sigma_T/m \sim 3\ \mathrm{cm^2\,g^{-1}}\);
- velocity-dependent / anisotropic benchmarks with parameter pairs roughly of the form
  - \((w,u)=(500,1000)\ \mathrm{km\,s^{-1}}\),
  - \((1000,2000)\ \mathrm{km\,s^{-1}}\),
  - \((1600,2000)\ \mathrm{km\,s^{-1}}\).

The main idea is that cluster profiles and subhalo distributions can discriminate among these elastic SIDM scenarios.

### 4.2 Adhikari et al. 2024

The benchmark models are explicitly listed.

**Elastic, isotropic SIDM:**

- \(\sigma/m = 0.2\ \mathrm{cm^2\,g^{-1}}\)
- \(\sigma/m = 0.5\ \mathrm{cm^2\,g^{-1}}\)
- \(\sigma/m = 1.0\ \mathrm{cm^2\,g^{-1}}\)
- \(\sigma/m = 2.0\ \mathrm{cm^2\,g^{-1}}\)

**Dissipative SIDM:**

- fixed \(\sigma/m = \sigma'/m = 1\ \mathrm{cm^2\,g^{-1}}\)
- cooling / loss-velocity benchmarks
  - \(\nu_{\rm loss}=300\ \mathrm{km\,s^{-1}}\)
  - \(\nu_{\rm loss}=600\ \mathrm{km\,s^{-1}}\)
  - \(\nu_{\rm loss}=2000\ \mathrm{km\,s^{-1}}\)

---

## 5. Current constraints from these papers and nearby context

### 5.1 Cluster weak-lensing constraint from Adhikari et al. 2024

For **elastic, isotropic SIDM**, the paper reports approximately:

- **95% CL:** \(\sigma/m < 1.05\ \mathrm{cm^2\,g^{-1}}\)
- **67% CL:** \(\sigma/m < 0.5\ \mathrm{cm^2\,g^{-1}}\)

This is a useful current benchmark for cluster-scale weak-lensing constraints.

### 5.2 Relation to older cluster constraints

The Banerjee et al. 2019 paper argued that weak-lensing constraints from cluster halo profiles should already be competitive with, but somewhat weaker than, classic Bullet-Cluster-style limits, i.e. at about the **order-unity cm\(^2\)/g** level.

A reasonable compact statement for science-case work is therefore:

> Current cluster-scale constraints on **elastic SIDM** are at roughly the \(\sigma/m \lesssim 1\ \mathrm{cm^2\,g^{-1}}\) level from stacked weak lensing, with some alternative cluster probes reaching tighter but more model-dependent bounds.

### 5.3 Dissipative SIDM status

Adhikari et al. also show that current cluster weak-lensing data have sensitivity to **dissipative SIDM** benchmark models, but the results are best treated as **benchmark-model comparisons**, not yet as a clean universal bound analogous to the elastic \(\sigma/m\) constraint.

---

## 6. Implications for a fast forecasting pipeline

For a one-day Stage-V science-case forecast:

1. Use `parametricSIDM` only for **elastic SIDM**.
2. Benchmark against the Adhikari et al. elastic grid first:
   - \(\sigma/m = 0.2, 0.5, 1.0, 2.0\ \mathrm{cm^2\,g^{-1}}\).
3. Build projected profiles in terms of
   - \(\rho(r) \rightarrow \Sigma(R) \rightarrow \Delta\Sigma(R)\).
4. Prefer a **halo ensemble / stacked forecast** over a single-halo curve if trying to compare to the cluster-lensing literature.
5. Treat Banerjee et al.'s velocity-dependent models as **effective cross-section targets**, not exact scattering-kernel reproductions.
6. Do not attempt to reproduce the dissipative models without adding a separate fluid-cooling component.

---

## 7. Recommended validation targets

For an initial implementation, the most useful validation targets are:

### Tier 1: internal validation

- recover sensible CDM NFW-like \(\Delta\Sigma(R)\) behavior;
- verify SIDM-induced suppression or enhancement in the inner projected signal depending on the benchmark model;
- confirm monotonic trends with increasing elastic \(\sigma/m\) where expected.

### Tier 2: literature-facing validation

- compare projected **elastic SIDM / CDM ratios** against the qualitative trends in Adhikari et al. 2024;
- compare stacked **3-D density trends** against Banerjee et al. 2019;
- do not claim exact reproduction unless the simulated halo sample and stacking details are matched.

---

## References

1. Banerjee et al. 2019 / JCAP 2020, *Signatures of Self-Interacting dark matter on cluster density profile and subhalo distributions*.
2. Adhikari et al. 2024, *Constraints on Dark Matter Self-interactions from Weak Lensing of Galaxies from the Dark Energy Survey around Clusters from the Atacama Cosmology Telescope Survey*.
3. Yang et al. 2023/2024, *A parametric model for self-interacting dark matter halos*.
4. Yang et al. 2024, *Exploring Self-Interacting Dark Matter Halos with Diverse Galaxy Populations*. 
