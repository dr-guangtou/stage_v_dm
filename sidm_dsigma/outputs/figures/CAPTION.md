# Figure Captions

This document records unambiguous captions for each figure in this folder.

Timestamp convention used below:
- `Created timestamp` is the file-system timestamp reported by `stat` on this machine (timezone shown explicitly).
- If a figure is regenerated, its timestamp and caption entry should be updated/appended accordingly.

Common context for the benchmark forecast figures:
- Mass definition: `M200c`.
- Redshift: `z = 0.3`.
- SIDM grid used in the latest benchmark run: `sigma_over_m = [0.0, 0.2, 0.5, 1.0, 2.0] cm^2/g`.
- `sigma_over_m = 0.0` corresponds to CDM reference.

## `dwarf_rho_profiles.png`
- File: `outputs/figures/dwarf_rho_profiles.png`
- Created timestamp: `2026-03-19 02:33:05 MST`
- Caption:
  This figure compares 3D density profiles for the dwarf benchmark halo (`M200 = 1e10 Msun`, `c200 = 15`, `z = 0.3`). The left panel shows `rho(r)` in physical units (`Msun/kpc^3`) for CDM and each SIDM benchmark cross section, plotted versus physical radius `r` in `kpc` on log-log axes; the right panel shows the ratio `rho_SIDM / rho_CDM` versus `r`, where values below 1 indicate central density suppression relative to CDM and values near 1 indicate minimal modification.

## `cluster_rho_profiles.png`
- File: `outputs/figures/cluster_rho_profiles.png`
- Created timestamp: `2026-03-19 02:33:06 MST`
- Caption:
  This figure compares 3D density profiles for the cluster benchmark halo (`M200 = 1e14 Msun`, `c200 = 4`, `z = 0.3`). The left panel shows `rho(r)` (`Msun/kpc^3`) for CDM and SIDM models across physical radius in `kpc` (log-log); the right panel shows `rho_SIDM / rho_CDM` as a function of radius, isolating where SIDM-driven deviations from the CDM baseline are strongest.

## `dwarf_delta_sigma_profiles.png`
- File: `outputs/figures/dwarf_delta_sigma_profiles.png`
- Created timestamp: `2026-03-19 02:33:06 MST`
- Caption:
  This figure shows projected weak-lensing predictions for the dwarf benchmark (`M200 = 1e10 Msun`). The left panel plots `DeltaSigma(R)` (`Msun/kpc^2`) versus projected radius `R` in `kpc` for CDM and SIDM models, with a shaded gray band representing the toy baseline `1 sigma` uncertainty model applied to the CDM curve (`sigma_i = f_i * |DeltaSigma_CDM|`); the right panel shows `DeltaSigma_SIDM / DeltaSigma_CDM` versus `R`, highlighting scale-dependent fractional lensing differences used in the distinguishability forecast.

## `cluster_delta_sigma_profiles.png`
- File: `outputs/figures/cluster_delta_sigma_profiles.png`
- Created timestamp: `2026-03-19 02:33:06 MST`
- Caption:
  This figure shows projected weak-lensing predictions for the cluster benchmark (`M200 = 1e14 Msun`). The left panel presents `DeltaSigma(R)` in `Msun/kpc^2` for CDM and SIDM profiles with the toy baseline CDM `1 sigma` band overlaid; the right panel presents `DeltaSigma_SIDM / DeltaSigma_CDM`, which directly shows the radial fractional signal shifts used to compute forecasted separation metrics.

## `nfw_reference_crosscheck.png`
- File: `outputs/figures/nfw_reference_crosscheck.png`
- Created timestamp: `2026-03-19 02:33:09 MST`
- Caption:
  This figure validates the local numerical projection implementation against reference NFW behavior. The left panel compares local numeric `DeltaSigma(R)` to analytic NFW `DeltaSigma(R)` for a standard NFW halo (`M200 = 1e14 Msun`, `c200 = 4`, `z = 0.3`); if optional external backends are available, their curves are also shown. The right panel shows fractional residuals in percent (e.g., `100 * (local/reference - 1)`), providing a direct normalization and shape-consistency diagnostic for the projection kernel.

## `dwarf_precision_sweep.png`
- File: `outputs/figures/dwarf_precision_sweep.png`
- Created timestamp: `2026-03-19 02:33:11 MST`
- Caption:
  This figure summarizes required measurement precision for separating SIDM from CDM in the dwarf benchmark as a function of `sigma/m`. Each panel corresponds to a radial window (`full`, `inner`, `outer`), and each line corresponds to a target detection significance (`2 sigma`, `3 sigma`, `5 sigma`). The y-axis is the required uniform per-bin fractional uncertainty on `DeltaSigma` (in percent), defined so that `sigma_i = f * |DeltaSigma_CDM(R_i)|` for all bins in the window and `sum[(DeltaSigma_SIDM - DeltaSigma_CDM)^2 / sigma_i^2] = target_sigma^2`.

## `cluster_precision_sweep.png`
- File: `outputs/figures/cluster_precision_sweep.png`
- Created timestamp: `2026-03-19 02:33:11 MST`
- Caption:
  This figure summarizes required measurement precision for SIDM-vs-CDM separation in the cluster benchmark. As in the dwarf figure, panels are radial windows (`full`, `inner`, `outer`) and curves are significance targets (`2 sigma`, `3 sigma`, `5 sigma`), with y-axis values giving the required uniform per-bin fractional `DeltaSigma` uncertainty percentage needed to achieve the target separation.

## `ensemble_stacked_delta_sigma.png`
- File: `outputs/figures/ensemble_stacked_delta_sigma.png`
- Created timestamp: `2026-03-19 03:19:03 MST`
- Caption:
  This figure shows Tier-1 ensemble-stacked lensing predictions for a sampled cluster-like halo population (`N_halos = 100`, log-normal mass distribution with mean `3e14 Msun`, scatter `0.2 dex`, `z = 0.4`, equal weights). The left panel plots stacked `DeltaSigma(R)` in `Msun/kpc^2` for CDM and SIDM (`sigma/m = 0.2, 0.5, 1.0, 2.0 cm^2/g`) as a function of projected physical radius `R` (`kpc`); each SIDM curve is produced by halo-by-halo profile generation with `parametricSIDM`, projection, interpolation to a common `R` grid, and weighted stacking. The right panel shows `DeltaSigma_SIDM / DeltaSigma_CDM`, isolating scale-dependent fractional departures from CDM after ensemble averaging.

## `ensemble_tier1_summary.png`
- File: `outputs/figures/ensemble_tier1_summary.png`
- Created timestamp: `2026-03-19 03:19:03 MST`
- Caption:
  This is the Tier-1 four-panel summary figure for ensemble forecasting. Top-left: stacked 3D density `rho(r)` (`Msun/kpc^3`) for CDM and SIDM over a common physical `r` grid after interpolating and averaging individual halo profiles. Top-right: stacked `DeltaSigma(R)` (`Msun/kpc^2`) for CDM and SIDM over the common projected `R` grid. Bottom-left: stacked ratio `DeltaSigma_SIDM / DeltaSigma_CDM`, showing whether the inner-halo SIDM imprint survives ensemble mixing. Bottom-right: `Delta chi^2` versus `sigma/m` (computed under the toy optimistic `5%` fractional-error model with diagonal covariance), which quantifies distinguishability strength as SIDM cross section increases.

## `cluster_ensemble_stacked_delta_sigma.png`
- File: `outputs/figures/cluster_ensemble_stacked_delta_sigma.png`
- Created timestamp: `2026-03-19 03:49:25 MST`
- Caption:
  This figure is the YAML-configured cluster (HMF-mode) stacked lensing output generated from `docs/cluster_ensemble_config.yaml`. Halos are sampled from an HMF-driven distribution with logistic selection and concentration scatter, then projected and stacked with equal weights on the configured radial grid (`R = 50-2000 kpc`, `N_R = 30`). The left panel shows stacked `DeltaSigma(R)` (`Msun/kpc^2`) for CDM and SIDM (`sigma/m = 0.2, 0.5, 1.0, 2.0`), and the right panel shows `DeltaSigma_SIDM / DeltaSigma_CDM` as a function of `R`.

## `cluster_ensemble_tier1_summary.png`
- File: `outputs/figures/cluster_ensemble_tier1_summary.png`
- Created timestamp: `2026-03-19 03:49:26 MST`
- Caption:
  This 4-panel summary corresponds to the cluster YAML run (HMF mode). Top-left: stacked `rho(r)` for CDM and SIDM. Top-right: stacked `DeltaSigma(R)` on the cluster projection grid. Bottom-left: stacked `DeltaSigma` ratio relative to CDM. Bottom-right: `Delta chi^2` versus `sigma/m` using the toy optimistic 5% fractional-error model. This figure summarizes whether SIDM signatures survive stacking across the selected cluster-like halo population.

## `dwarf_ensemble_stacked_delta_sigma.png`
- File: `outputs/figures/dwarf_ensemble_stacked_delta_sigma.png`
- Created timestamp: `2026-03-19 03:50:47 MST`
- Caption:
  This figure is the YAML-configured dwarf (SHMR-mode) stacked lensing output generated from `docs/dwarf_ensemble_config.yaml`. Galaxies are sampled in stellar mass, mapped to halos through a Tier-1 SHMR model with halo scatter, and stacked with equal weights on the configured projected radial grid (`R = 1-300 kpc`, `N_R = 30`). The left panel shows stacked `DeltaSigma(R)` (`Msun/kpc^2`) for CDM and SIDM (`sigma/m = 0.5, 1.0, 3.0, 10.0`), and the right panel shows `DeltaSigma_SIDM / DeltaSigma_CDM`.

## `dwarf_ensemble_tier1_summary.png`
- File: `outputs/figures/dwarf_ensemble_tier1_summary.png`
- Created timestamp: `2026-03-19 03:50:48 MST`
- Caption:
  This 4-panel summary corresponds to the dwarf YAML run (SHMR mode). Top-left shows stacked `rho(r)` for CDM and SIDM models inferred from the SHMR-selected ensemble. Top-right shows stacked `DeltaSigma(R)` over the dwarf radial bins. Bottom-left shows `DeltaSigma_SIDM / DeltaSigma_CDM`, and bottom-right shows `Delta chi^2` versus `sigma/m` (optimistic 5% toy precision), quantifying distinguishability in the low-mass regime.

## `cluster_dwarf_delta_chi2_comparison.png`
- File: `outputs/figures/cluster_dwarf_delta_chi2_comparison.png`
- Created timestamp: `2026-03-19 03:51:12 MST`
- Caption:
  This comparison plot overlays `Delta chi^2` (optimistic 5% scenario) versus `sigma/m` for the two YAML-configured ensemble pipelines: cluster HMF and dwarf SHMR. Both axes are logarithmic. The purpose is direct cross-regime comparison of stacked SIDM distinguishability after applying each regime’s selection-aware ensemble model and radial binning.

## Update 2026-03-19T07:56:13-07:00

### cluster_ensemble_stacked_delta_sigma.png
- Created: 2026-03-19T07:56:12.338006258-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

### cluster_ensemble_tier1_summary.png
- Created: 2026-03-19T07:56:12.816890717-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

### cluster_ensemble_tier1_vs_tier2_stacked.png
- Created: 2026-03-19T07:56:13.092409134-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

### cluster_ensemble_tier1_vs_tier2_delta_chi2.png
- Created: 2026-03-19T07:56:13.195337296-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

### cluster_single_halo_tier2_diagnostic.png
- Created: 2026-03-19T07:56:13.648265839-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

## Update 2026-03-19T07:56:55-07:00

### dwarf_ensemble_stacked_delta_sigma.png
- Created: 2026-03-19T07:56:53.958179951-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

### dwarf_ensemble_tier1_summary.png
- Created: 2026-03-19T07:56:54.607635975-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

### dwarf_ensemble_tier1_vs_tier2_stacked.png
- Created: 2026-03-19T07:56:54.954513550-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

### dwarf_ensemble_tier1_vs_tier2_delta_chi2.png
- Created: 2026-03-19T07:56:55.051311493-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

### dwarf_single_halo_tier2_diagnostic.png
- Created: 2026-03-19T07:56:55.494876385-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

## Update 2026-03-19T07:58:26-07:00 (Tier-2 Detailed Captions)

## `cluster_ensemble_tier1_vs_tier2_stacked.png`
- File: `outputs/figures/cluster_ensemble_tier1_vs_tier2_stacked.png`
- Created timestamp: `2026-03-19 07:56:13 MST`
- Caption:
  Cluster HMF ensemble stacked DeltaSigma comparison between Tier-1 and Tier-2. Left panel: physical projected radius R [kpc] versus stacked DeltaSigma [Msun/kpc^2] for CDM and SIDM (sigma/m = 0.2, 0.5, 1.0, 2.0), with Tier-1 as solid curves and Tier-2 as dashed curves. Right panel: multiplicative impact ratio DeltaSigma_Tier2 / DeltaSigma_Tier1 for the same models; values above unity indicate outskirts attachment boosts the projected signal at that radius.

## `cluster_ensemble_tier1_vs_tier2_delta_chi2.png`
- File: `outputs/figures/cluster_ensemble_tier1_vs_tier2_delta_chi2.png`
- Created timestamp: `2026-03-19 07:56:13 MST`
- Caption:
  Cluster HMF ensemble distinguishability comparison. X-axis is SIDM cross section sigma/m [cm^2/g], y-axis is Delta chi^2 under the optimistic 5% fractional-error toy model. Tier-1 (inner-only) and Tier-2 (hybrid outskirts) are plotted together to show how DK14-like outskirts treatment changes CDM-vs-SIDM separation strength.

## `cluster_single_halo_tier2_diagnostic.png`
- File: `outputs/figures/cluster_single_halo_tier2_diagnostic.png`
- Created timestamp: `2026-03-19 07:56:13 MST`
- Caption:
  Single representative cluster-halo Tier-2 diagnostic. Panel 1: rho(r) [Msun/kpc^3] for CDM inner NFW, SIDM inner profile, DK14-like reference, and stitched hybrid. Panel 2: logarithmic slope d ln rho / d ln r versus r to verify outer steepening and smooth stitching. Panel 3: projected DeltaSigma(R) [Msun/kpc^2] comparing Tier-1 inner-only and Tier-2 hybrid profiles to check projection smoothness and identify where outskirts attachment activates.

## `dwarf_ensemble_tier1_vs_tier2_stacked.png`
- File: `outputs/figures/dwarf_ensemble_tier1_vs_tier2_stacked.png`
- Created timestamp: `2026-03-19 07:56:54 MST`
- Caption:
  Dwarf SHMR ensemble stacked DeltaSigma comparison between Tier-1 and Tier-2. Left panel shows stacked DeltaSigma(R) [Msun/kpc^2] for CDM and SIDM (sigma/m = 0.2, 0.5, 1.0, 2.0) with solid Tier-1 curves and dashed Tier-2 curves. Right panel shows DeltaSigma_Tier2 / DeltaSigma_Tier1, quantifying the relative effect of outer-profile attachment in the dwarf-selected sample.

## `dwarf_ensemble_tier1_vs_tier2_delta_chi2.png`
- File: `outputs/figures/dwarf_ensemble_tier1_vs_tier2_delta_chi2.png`
- Created timestamp: `2026-03-19 07:56:55 MST`
- Caption:
  Dwarf SHMR ensemble distinguishability comparison. X-axis is sigma/m [cm^2/g], y-axis is Delta chi^2 for the optimistic 5% toy precision scenario. Tier-1 and Tier-2 curves are overlaid to assess whether inner SIDM leverage survives when DK14-like outskirts are attached in the low-mass regime.

## `dwarf_single_halo_tier2_diagnostic.png`
- File: `outputs/figures/dwarf_single_halo_tier2_diagnostic.png`
- Created timestamp: `2026-03-19 07:56:55 MST`
- Caption:
  Single representative dwarf-halo Tier-2 diagnostic with the same three-panel structure as the cluster diagnostic: density components, logarithmic slope, and Tier-1 vs Tier-2 projected DeltaSigma. This figure is used to verify that stitching remains smooth in the dwarf regime and that inner SIDM profile differences remain visible at small radii.

## Update 2026-03-19T07:58:52-07:00 (Tier-1 Baseline Captions Refreshed)

## `cluster_ensemble_stacked_delta_sigma.png`
- File: `outputs/figures/cluster_ensemble_stacked_delta_sigma.png`
- Created timestamp: `2026-03-19 07:56:12 MST`
- Caption:
  Cluster HMF stacked Tier-1 baseline lensing figure regenerated in this Tier-2 run. Left panel shows stacked DeltaSigma(R) [Msun/kpc^2] for CDM and SIDM with sigma/m = 0.2, 0.5, 1.0, 2.0. Right panel shows DeltaSigma_SIDM/DeltaSigma_CDM. This figure remains the inner-only Tier-1 reference for comparison against Tier-2 hybrid outputs.

## `cluster_ensemble_tier1_summary.png`
- File: `outputs/figures/cluster_ensemble_tier1_summary.png`
- Created timestamp: `2026-03-19 07:56:12 MST`
- Caption:
  Cluster HMF Tier-1 four-panel summary regenerated with the current sigma/m grid [0.2, 0.5, 1.0, 2.0]. Panels show stacked rho(r), stacked DeltaSigma(R), DeltaSigma ratio to CDM, and Delta chi^2 versus sigma/m (5% toy precision).

## `dwarf_ensemble_stacked_delta_sigma.png`
- File: `outputs/figures/dwarf_ensemble_stacked_delta_sigma.png`
- Created timestamp: `2026-03-19 07:56:53 MST`
- Caption:
  Dwarf SHMR stacked Tier-1 baseline lensing figure regenerated in this Tier-2 run. Left panel shows stacked DeltaSigma(R) [Msun/kpc^2] for CDM and SIDM using sigma/m = 0.2, 0.5, 1.0, 2.0 (updated grid). Right panel shows DeltaSigma_SIDM/DeltaSigma_CDM.

## `dwarf_ensemble_tier1_summary.png`
- File: `outputs/figures/dwarf_ensemble_tier1_summary.png`
- Created timestamp: `2026-03-19 07:56:54 MST`
- Caption:
  Dwarf SHMR Tier-1 four-panel summary regenerated with sigma/m grid [0.2, 0.5, 1.0, 2.0]. Panels show stacked rho(r), stacked DeltaSigma(R), DeltaSigma ratio to CDM, and Delta chi^2 versus sigma/m (5% toy precision).
