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

## Update 2026-03-19T11:47:44-07:00

### cluster_ensemble_stacked_delta_sigma.png
- Created: 2026-03-19T11:47:43.252906083-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

### cluster_ensemble_tier1_summary.png
- Created: 2026-03-19T11:47:43.821531534-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

### cluster_ensemble_tier1_vs_tier2_stacked.png
- Created: 2026-03-19T11:47:44.150659323-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

### cluster_ensemble_tier1_vs_tier2_delta_chi2.png
- Created: 2026-03-19T11:47:44.255169152-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

### cluster_single_halo_tier2_diagnostic.png
- Created: 2026-03-19T11:47:44.725646973-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

## Update 2026-03-19T11:49:50-07:00

### dwarf_ensemble_stacked_delta_sigma.png
- Created: 2026-03-19T11:49:48.590800285-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

### dwarf_ensemble_tier1_summary.png
- Created: 2026-03-19T11:49:49.237215519-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

### dwarf_ensemble_tier1_vs_tier2_stacked.png
- Created: 2026-03-19T11:49:49.600851536-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

### dwarf_ensemble_tier1_vs_tier2_delta_chi2.png
- Created: 2026-03-19T11:49:49.696293831-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

### dwarf_single_halo_tier2_diagnostic.png
- Created: 2026-03-19T11:49:50.179052830-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

## Update 2026-03-19T11:50:26-07:00

### cluster_ensemble_stacked_delta_sigma.png
- Created: 2026-03-19T11:50:24.965315342-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

### cluster_ensemble_tier1_summary.png
- Created: 2026-03-19T11:50:25.612401962-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

### cluster_ensemble_tier1_vs_tier2_stacked.png
- Created: 2026-03-19T11:50:25.956269503-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

### cluster_ensemble_tier1_vs_tier2_delta_chi2.png
- Created: 2026-03-19T11:50:26.058327198-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

### cluster_single_halo_tier2_diagnostic.png
- Created: 2026-03-19T11:50:26.563595533-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, while Tier-2 curves include DK14-like outskirts attachment with smooth log-density stitching.

## Update 2026-03-19T11:58:50-07:00

### cluster_ensemble_stacked_delta_sigma.png
- Created: 2026-03-19T11:58:45.995086432-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster_ensemble_tier1_summary.png
- Created: 2026-03-19T11:58:47.362810135-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster_ensemble_tier1_vs_tier2_stacked.png
- Created: 2026-03-19T11:58:47.796515226-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster_ensemble_tier1_vs_tier2_delta_chi2.png
- Created: 2026-03-19T11:58:47.916850090-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster_single_halo_tier2_diagnostic.png
- Created: 2026-03-19T11:58:48.524509430-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster_ensemble_tier3_summary.png
- Created: 2026-03-19T11:58:49.243898869-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster_ensemble_tier2_vs_tier3_stacked.png
- Created: 2026-03-19T11:58:49.792044401-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster_ensemble_tier2_vs_tier3_delta_chi2.png
- Created: 2026-03-19T11:58:49.915891886-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster_single_halo_tier3_diagnostic.png
- Created: 2026-03-19T11:58:50.675368547-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster_tier3_sensitivity_rt_shift.png
- Created: 2026-03-19T11:58:50.858799696-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

## Update 2026-03-19T11:59:08-07:00

### dwarf_ensemble_stacked_delta_sigma.png
- Created: 2026-03-19T11:59:04.854831696-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf_ensemble_tier1_summary.png
- Created: 2026-03-19T11:59:05.561036587-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf_ensemble_tier1_vs_tier2_stacked.png
- Created: 2026-03-19T11:59:05.948815107-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf_ensemble_tier1_vs_tier2_delta_chi2.png
- Created: 2026-03-19T11:59:06.046019554-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf_single_halo_tier2_diagnostic.png
- Created: 2026-03-19T11:59:06.590419531-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf_ensemble_tier3_summary.png
- Created: 2026-03-19T11:59:07.111136675-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf_ensemble_tier2_vs_tier3_stacked.png
- Created: 2026-03-19T11:59:07.644942760-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf_ensemble_tier2_vs_tier3_delta_chi2.png
- Created: 2026-03-19T11:59:07.744937658-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf_single_halo_tier3_diagnostic.png
- Created: 2026-03-19T11:59:08.409980774-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf_tier3_sensitivity_rt_shift.png
- Created: 2026-03-19T11:59:08.569684267-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

## Update 2026-03-19T12:14:00-07:00

### cluster1_z0p1_ensemble_stacked_delta_sigma.png
- Created: 2026-03-19T12:13:56.017127752-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p1_ensemble_tier1_summary.png
- Created: 2026-03-19T12:13:56.672052622-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p1_ensemble_tier1_vs_tier2_stacked.png
- Created: 2026-03-19T12:13:57.100394487-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p1_ensemble_tier1_vs_tier2_delta_chi2.png
- Created: 2026-03-19T12:13:57.217632532-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p1_single_halo_tier2_diagnostic.png
- Created: 2026-03-19T12:13:57.876997948-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p1_ensemble_tier3_summary.png
- Created: 2026-03-19T12:13:58.695195913-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p1_ensemble_tier2_vs_tier3_stacked.png
- Created: 2026-03-19T12:13:59.137138128-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p1_ensemble_tier2_vs_tier3_delta_chi2.png
- Created: 2026-03-19T12:13:59.228872776-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p1_single_halo_tier3_diagnostic.png
- Created: 2026-03-19T12:13:59.852906227-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p1_tier3_sensitivity_rt_shift.png
- Created: 2026-03-19T12:14:00.106198072-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

## Update 2026-03-19T12:14:43-07:00

### cluster1_z0p5_ensemble_stacked_delta_sigma.png
- Created: 2026-03-19T12:14:39.561174631-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p5_ensemble_tier1_summary.png
- Created: 2026-03-19T12:14:40.160544634-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p5_ensemble_tier1_vs_tier2_stacked.png
- Created: 2026-03-19T12:14:40.523278712-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p5_ensemble_tier1_vs_tier2_delta_chi2.png
- Created: 2026-03-19T12:14:40.623074055-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p5_single_halo_tier2_diagnostic.png
- Created: 2026-03-19T12:14:41.187230825-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p5_ensemble_tier3_summary.png
- Created: 2026-03-19T12:14:41.812653065-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p5_ensemble_tier2_vs_tier3_stacked.png
- Created: 2026-03-19T12:14:42.312062502-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p5_ensemble_tier2_vs_tier3_delta_chi2.png
- Created: 2026-03-19T12:14:42.405933380-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p5_single_halo_tier3_diagnostic.png
- Created: 2026-03-19T12:14:43.119284391-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p5_tier3_sensitivity_rt_shift.png
- Created: 2026-03-19T12:14:43.283358335-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

## Update 2026-03-19T12:15:25-07:00

### cluster2_z0p1_ensemble_stacked_delta_sigma.png
- Created: 2026-03-19T12:15:22.088594913-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p1_ensemble_tier1_summary.png
- Created: 2026-03-19T12:15:22.647788286-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p1_ensemble_tier1_vs_tier2_stacked.png
- Created: 2026-03-19T12:15:23.042280912-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p1_ensemble_tier1_vs_tier2_delta_chi2.png
- Created: 2026-03-19T12:15:23.209324360-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p1_single_halo_tier2_diagnostic.png
- Created: 2026-03-19T12:15:23.694726944-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p1_ensemble_tier3_summary.png
- Created: 2026-03-19T12:15:24.300121307-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p1_ensemble_tier2_vs_tier3_stacked.png
- Created: 2026-03-19T12:15:24.732865334-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p1_ensemble_tier2_vs_tier3_delta_chi2.png
- Created: 2026-03-19T12:15:24.834270239-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p1_single_halo_tier3_diagnostic.png
- Created: 2026-03-19T12:15:25.539442539-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p1_tier3_sensitivity_rt_shift.png
- Created: 2026-03-19T12:15:25.686833858-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

## Update 2026-03-19T12:16:06-07:00

### cluster2_z0p5_ensemble_stacked_delta_sigma.png
- Created: 2026-03-19T12:16:03.126037598-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p5_ensemble_tier1_summary.png
- Created: 2026-03-19T12:16:03.716132879-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p5_ensemble_tier1_vs_tier2_stacked.png
- Created: 2026-03-19T12:16:04.050545216-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p5_ensemble_tier1_vs_tier2_delta_chi2.png
- Created: 2026-03-19T12:16:04.147610903-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p5_single_halo_tier2_diagnostic.png
- Created: 2026-03-19T12:16:04.657881737-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p5_ensemble_tier3_summary.png
- Created: 2026-03-19T12:16:05.245115519-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p5_ensemble_tier2_vs_tier3_stacked.png
- Created: 2026-03-19T12:16:05.687441826-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p5_ensemble_tier2_vs_tier3_delta_chi2.png
- Created: 2026-03-19T12:16:05.783334494-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p5_single_halo_tier3_diagnostic.png
- Created: 2026-03-19T12:16:06.467840672-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p5_tier3_sensitivity_rt_shift.png
- Created: 2026-03-19T12:16:06.647597551-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

## Update 2026-03-19T12:17:01-07:00

### dwarf1_z0p05_ensemble_stacked_delta_sigma.png
- Created: 2026-03-19T12:16:58.355311871-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p05_ensemble_tier1_summary.png
- Created: 2026-03-19T12:16:58.951833725-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p05_ensemble_tier1_vs_tier2_stacked.png
- Created: 2026-03-19T12:16:59.278313875-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p05_ensemble_tier1_vs_tier2_delta_chi2.png
- Created: 2026-03-19T12:16:59.357883930-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p05_single_halo_tier2_diagnostic.png
- Created: 2026-03-19T12:16:59.797064304-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p05_ensemble_tier3_summary.png
- Created: 2026-03-19T12:17:00.245533943-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p05_ensemble_tier2_vs_tier3_stacked.png
- Created: 2026-03-19T12:17:00.660478115-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p05_ensemble_tier2_vs_tier3_delta_chi2.png
- Created: 2026-03-19T12:17:00.744708061-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p05_single_halo_tier3_diagnostic.png
- Created: 2026-03-19T12:17:01.305227995-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p05_tier3_sensitivity_rt_shift.png
- Created: 2026-03-19T12:17:01.445106506-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

## Update 2026-03-19T12:17:53-07:00

### dwarf1_z0p2_ensemble_stacked_delta_sigma.png
- Created: 2026-03-19T12:17:50.888899803-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p2_ensemble_tier1_summary.png
- Created: 2026-03-19T12:17:51.437226057-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p2_ensemble_tier1_vs_tier2_stacked.png
- Created: 2026-03-19T12:17:51.776608706-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p2_ensemble_tier1_vs_tier2_delta_chi2.png
- Created: 2026-03-19T12:17:51.856777430-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p2_single_halo_tier2_diagnostic.png
- Created: 2026-03-19T12:17:52.313448668-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p2_ensemble_tier3_summary.png
- Created: 2026-03-19T12:17:52.747900486-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p2_ensemble_tier2_vs_tier3_stacked.png
- Created: 2026-03-19T12:17:53.176605225-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p2_ensemble_tier2_vs_tier3_delta_chi2.png
- Created: 2026-03-19T12:17:53.262454509-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p2_single_halo_tier3_diagnostic.png
- Created: 2026-03-19T12:17:53.812981844-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p2_tier3_sensitivity_rt_shift.png
- Created: 2026-03-19T12:17:53.949257135-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

## Update 2026-03-19T12:18:46-07:00

### dwarf2_z0p05_ensemble_stacked_delta_sigma.png
- Created: 2026-03-19T12:18:42.827586651-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p05_ensemble_tier1_summary.png
- Created: 2026-03-19T12:18:43.392032862-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p05_ensemble_tier1_vs_tier2_stacked.png
- Created: 2026-03-19T12:18:43.737237453-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p05_ensemble_tier1_vs_tier2_delta_chi2.png
- Created: 2026-03-19T12:18:43.826119900-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p05_single_halo_tier2_diagnostic.png
- Created: 2026-03-19T12:18:44.277556896-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p05_ensemble_tier3_summary.png
- Created: 2026-03-19T12:18:44.813823700-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p05_ensemble_tier2_vs_tier3_stacked.png
- Created: 2026-03-19T12:18:45.253279209-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p05_ensemble_tier2_vs_tier3_delta_chi2.png
- Created: 2026-03-19T12:18:45.342742443-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p05_single_halo_tier3_diagnostic.png
- Created: 2026-03-19T12:18:45.897559643-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p05_tier3_sensitivity_rt_shift.png
- Created: 2026-03-19T12:18:46.035641432-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

## Update 2026-03-19T12:19:38-07:00

### dwarf2_z0p2_ensemble_stacked_delta_sigma.png
- Created: 2026-03-19T12:19:34.908510923-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p2_ensemble_tier1_summary.png
- Created: 2026-03-19T12:19:35.461121798-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p2_ensemble_tier1_vs_tier2_stacked.png
- Created: 2026-03-19T12:19:35.806502819-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p2_ensemble_tier1_vs_tier2_delta_chi2.png
- Created: 2026-03-19T12:19:35.895833015-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p2_single_halo_tier2_diagnostic.png
- Created: 2026-03-19T12:19:36.338086605-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p2_ensemble_tier3_summary.png
- Created: 2026-03-19T12:19:36.871542692-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p2_ensemble_tier2_vs_tier3_stacked.png
- Created: 2026-03-19T12:19:37.307467937-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p2_ensemble_tier2_vs_tier3_delta_chi2.png
- Created: 2026-03-19T12:19:37.396914721-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p2_single_halo_tier3_diagnostic.png
- Created: 2026-03-19T12:19:37.944637775-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p2_tier3_sensitivity_rt_shift.png
- Created: 2026-03-19T12:19:38.089459181-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

## Update 2026-03-19T12:33:21-07:00

### cluster1_ensemble_tier3_summary.png
- Created: 2026-03-19T12:33:19.776684046-07:00
- Caption: Tier-3 summary for cluster1 with redshift overlay. Two line styles correspond to the two fixed-redshift runs in that ensemble split; color encodes sigma/m. Top-right panel uses DeltaSigma units of h Msun / pc^2. Annotation box lists stellar-mass range (or N/A for HMF), halo-mass range definition, and concentration scatter used in the ensemble configuration.

### cluster2_ensemble_tier3_summary.png
- Created: 2026-03-19T12:33:20.332443714-07:00
- Caption: Tier-3 summary for cluster2 with redshift overlay. Two line styles correspond to the two fixed-redshift runs in that ensemble split; color encodes sigma/m. Top-right panel uses DeltaSigma units of h Msun / pc^2. Annotation box lists stellar-mass range (or N/A for HMF), halo-mass range definition, and concentration scatter used in the ensemble configuration.

### dwarf1_ensemble_tier3_summary.png
- Created: 2026-03-19T12:33:20.839016676-07:00
- Caption: Tier-3 summary for dwarf1 with redshift overlay. Two line styles correspond to the two fixed-redshift runs in that ensemble split; color encodes sigma/m. Top-right panel uses DeltaSigma units of h Msun / pc^2. Annotation box lists stellar-mass range (or N/A for HMF), halo-mass range definition, and concentration scatter used in the ensemble configuration.

### dwarf2_ensemble_tier3_summary.png
- Created: 2026-03-19T12:33:21.428270340-07:00
- Caption: Tier-3 summary for dwarf2 with redshift overlay. Two line styles correspond to the two fixed-redshift runs in that ensemble split; color encodes sigma/m. Top-right panel uses DeltaSigma units of h Msun / pc^2. Annotation box lists stellar-mass range (or N/A for HMF), halo-mass range definition, and concentration scatter used in the ensemble configuration.

## Update 2026-03-19T12:41:14-07:00

### cdm_engine_convergence_parametric_vs_colossus_dk14.png
- Created: 2026-03-19T12:41:14.364957571-07:00
- Caption: Convergence diagnostic comparing CDM profiles from parametricSIDM at tiny cross section (sigma/m=1e-4 cm^2/g) versus colossus DK14 CDM. Panels show density and DeltaSigma curves and their ratios for dwarf and cluster benchmarks.

## Update 2026-03-19T12:42:40-07:00

### cdm_engine_convergence_parametric_vs_colossus_dk14.png
- Created: 2026-03-19T12:42:40.920594454-07:00
- Caption: Convergence diagnostic comparing CDM profiles from parametricSIDM at tiny cross section (sigma/m=1e-4 cm^2/g) versus colossus DK14 CDM. Panels show density and DeltaSigma curves and their ratios for dwarf and cluster benchmarks.

## Update 2026-03-19T12:43:24-07:00

### cdm_engine_convergence_parametric_vs_colossus_dk14.png
- Created: 2026-03-19T12:43:24.135110378-07:00
- Caption: Convergence diagnostic comparing CDM profiles from parametricSIDM at tiny cross section (sigma/m=1e-4 cm^2/g) versus colossus DK14 CDM. Panels show density and DeltaSigma curves and their ratios for dwarf and cluster benchmarks. Ratio panels use the convergence window 0.05*R200c to 1.0*R200c.

## Update 2026-03-19T14:22:25-07:00

### cluster1_z0p1_ensemble_stacked_delta_sigma.png
- Created: 2026-03-19T14:22:21.830804348-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p1_ensemble_tier1_summary.png
- Created: 2026-03-19T14:22:22.357378483-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p1_ensemble_tier1_vs_tier2_stacked.png
- Created: 2026-03-19T14:22:22.715981245-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p1_ensemble_tier1_vs_tier2_delta_chi2.png
- Created: 2026-03-19T14:22:22.828362226-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p1_single_halo_tier2_diagnostic.png
- Created: 2026-03-19T14:22:23.287818909-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p1_ensemble_tier3_summary.png
- Created: 2026-03-19T14:22:23.814328671-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p1_ensemble_tier2_vs_tier3_stacked.png
- Created: 2026-03-19T14:22:24.825147152-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p1_ensemble_tier2_vs_tier3_delta_chi2.png
- Created: 2026-03-19T14:22:24.932975054-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p1_single_halo_tier3_diagnostic.png
- Created: 2026-03-19T14:22:25.541112423-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p1_tier3_sensitivity_rt_shift.png
- Created: 2026-03-19T14:22:25.685726404-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

## Update 2026-03-19T14:22:59-07:00

### cluster1_z0p5_ensemble_stacked_delta_sigma.png
- Created: 2026-03-19T14:22:56.069368601-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p5_ensemble_tier1_summary.png
- Created: 2026-03-19T14:22:56.607275486-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p5_ensemble_tier1_vs_tier2_stacked.png
- Created: 2026-03-19T14:22:56.871472359-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p5_ensemble_tier1_vs_tier2_delta_chi2.png
- Created: 2026-03-19T14:22:56.977162838-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p5_single_halo_tier2_diagnostic.png
- Created: 2026-03-19T14:22:57.497240543-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p5_ensemble_tier3_summary.png
- Created: 2026-03-19T14:22:58.023528337-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p5_ensemble_tier2_vs_tier3_stacked.png
- Created: 2026-03-19T14:22:58.374357939-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p5_ensemble_tier2_vs_tier3_delta_chi2.png
- Created: 2026-03-19T14:22:58.487547874-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p5_single_halo_tier3_diagnostic.png
- Created: 2026-03-19T14:22:58.998214960-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster1_z0p5_tier3_sensitivity_rt_shift.png
- Created: 2026-03-19T14:22:59.160977840-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

## Update 2026-03-19T14:23:31-07:00

### cluster2_z0p1_ensemble_stacked_delta_sigma.png
- Created: 2026-03-19T14:23:28.267508030-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p1_ensemble_tier1_summary.png
- Created: 2026-03-19T14:23:28.812824249-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p1_ensemble_tier1_vs_tier2_stacked.png
- Created: 2026-03-19T14:23:29.113665581-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p1_ensemble_tier1_vs_tier2_delta_chi2.png
- Created: 2026-03-19T14:23:29.206602573-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p1_single_halo_tier2_diagnostic.png
- Created: 2026-03-19T14:23:29.611661196-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p1_ensemble_tier3_summary.png
- Created: 2026-03-19T14:23:30.368910074-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p1_ensemble_tier2_vs_tier3_stacked.png
- Created: 2026-03-19T14:23:30.800246239-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p1_ensemble_tier2_vs_tier3_delta_chi2.png
- Created: 2026-03-19T14:23:30.888028383-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p1_single_halo_tier3_diagnostic.png
- Created: 2026-03-19T14:23:31.434416294-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p1_tier3_sensitivity_rt_shift.png
- Created: 2026-03-19T14:23:31.579033375-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

## Update 2026-03-19T14:24:03-07:00

### cluster2_z0p5_ensemble_stacked_delta_sigma.png
- Created: 2026-03-19T14:24:00.679798603-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p5_ensemble_tier1_summary.png
- Created: 2026-03-19T14:24:01.200346947-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p5_ensemble_tier1_vs_tier2_stacked.png
- Created: 2026-03-19T14:24:01.486865997-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p5_ensemble_tier1_vs_tier2_delta_chi2.png
- Created: 2026-03-19T14:24:01.572557449-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p5_single_halo_tier2_diagnostic.png
- Created: 2026-03-19T14:24:01.974555492-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p5_ensemble_tier3_summary.png
- Created: 2026-03-19T14:24:02.496343851-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p5_ensemble_tier2_vs_tier3_stacked.png
- Created: 2026-03-19T14:24:02.843079329-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p5_ensemble_tier2_vs_tier3_delta_chi2.png
- Created: 2026-03-19T14:24:02.964863300-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p5_single_halo_tier3_diagnostic.png
- Created: 2026-03-19T14:24:03.489757061-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### cluster2_z0p5_tier3_sensitivity_rt_shift.png
- Created: 2026-03-19T14:24:03.656665802-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

## Update 2026-03-19T14:24:55-07:00

### dwarf1_z0p05_ensemble_stacked_delta_sigma.png
- Created: 2026-03-19T14:24:51.212369680-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p05_ensemble_tier1_summary.png
- Created: 2026-03-19T14:24:52.163012981-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p05_ensemble_tier1_vs_tier2_stacked.png
- Created: 2026-03-19T14:24:52.645520210-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p05_ensemble_tier1_vs_tier2_delta_chi2.png
- Created: 2026-03-19T14:24:52.739600420-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p05_single_halo_tier2_diagnostic.png
- Created: 2026-03-19T14:24:53.294876337-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p05_ensemble_tier3_summary.png
- Created: 2026-03-19T14:24:53.728591919-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p05_ensemble_tier2_vs_tier3_stacked.png
- Created: 2026-03-19T14:24:54.198842525-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p05_ensemble_tier2_vs_tier3_delta_chi2.png
- Created: 2026-03-19T14:24:54.283718348-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p05_single_halo_tier3_diagnostic.png
- Created: 2026-03-19T14:24:54.854867220-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p05_tier3_sensitivity_rt_shift.png
- Created: 2026-03-19T14:24:54.994208336-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

## Update 2026-03-19T14:25:45-07:00

### dwarf1_z0p2_ensemble_stacked_delta_sigma.png
- Created: 2026-03-19T14:25:42.323761225-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p2_ensemble_tier1_summary.png
- Created: 2026-03-19T14:25:42.904618740-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p2_ensemble_tier1_vs_tier2_stacked.png
- Created: 2026-03-19T14:25:43.257781982-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p2_ensemble_tier1_vs_tier2_delta_chi2.png
- Created: 2026-03-19T14:25:43.352015972-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p2_single_halo_tier2_diagnostic.png
- Created: 2026-03-19T14:25:43.852229357-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p2_ensemble_tier3_summary.png
- Created: 2026-03-19T14:25:44.299171686-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p2_ensemble_tier2_vs_tier3_stacked.png
- Created: 2026-03-19T14:25:44.778087854-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p2_ensemble_tier2_vs_tier3_delta_chi2.png
- Created: 2026-03-19T14:25:44.861454725-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p2_single_halo_tier3_diagnostic.png
- Created: 2026-03-19T14:25:45.445489883-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf1_z0p2_tier3_sensitivity_rt_shift.png
- Created: 2026-03-19T14:25:45.586649179-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

## Update 2026-03-19T14:26:36-07:00

### dwarf2_z0p05_ensemble_stacked_delta_sigma.png
- Created: 2026-03-19T14:26:33.422376633-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p05_ensemble_tier1_summary.png
- Created: 2026-03-19T14:26:33.991962194-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p05_ensemble_tier1_vs_tier2_stacked.png
- Created: 2026-03-19T14:26:34.309398651-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p05_ensemble_tier1_vs_tier2_delta_chi2.png
- Created: 2026-03-19T14:26:34.391035080-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p05_single_halo_tier2_diagnostic.png
- Created: 2026-03-19T14:26:34.799377918-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p05_ensemble_tier3_summary.png
- Created: 2026-03-19T14:26:35.339258671-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p05_ensemble_tier2_vs_tier3_stacked.png
- Created: 2026-03-19T14:26:35.754954338-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p05_ensemble_tier2_vs_tier3_delta_chi2.png
- Created: 2026-03-19T14:26:35.842357159-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p05_single_halo_tier3_diagnostic.png
- Created: 2026-03-19T14:26:36.444833040-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p05_tier3_sensitivity_rt_shift.png
- Created: 2026-03-19T14:26:36.586230993-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

## Update 2026-03-19T14:27:28-07:00

### dwarf2_z0p2_ensemble_stacked_delta_sigma.png
- Created: 2026-03-19T14:27:25.338500977-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p2_ensemble_tier1_summary.png
- Created: 2026-03-19T14:27:25.982053280-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p2_ensemble_tier1_vs_tier2_stacked.png
- Created: 2026-03-19T14:27:26.372807264-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p2_ensemble_tier1_vs_tier2_delta_chi2.png
- Created: 2026-03-19T14:27:26.481457949-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p2_single_halo_tier2_diagnostic.png
- Created: 2026-03-19T14:27:26.911588192-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p2_ensemble_tier3_summary.png
- Created: 2026-03-19T14:27:27.563514471-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p2_ensemble_tier2_vs_tier3_stacked.png
- Created: 2026-03-19T14:27:28.027154446-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p2_ensemble_tier2_vs_tier3_delta_chi2.png
- Created: 2026-03-19T14:27:28.112488747-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p2_single_halo_tier3_diagnostic.png
- Created: 2026-03-19T14:27:28.677552462-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

### dwarf2_z0p2_tier3_sensitivity_rt_shift.png
- Created: 2026-03-19T14:27:28.817254305-07:00
- Caption: Tier forecast figure. Axes are in physical kpc and Msun-based density units; Tier-1 curves use inner-only profiles, Tier-2 curves include DK14-like outskirts attachment, and Tier-3 curves include empirical SIDM outer corrections.

## Update 2026-03-19T14:28:04-07:00

### cluster1_ensemble_tier3_summary.png
- Created: 2026-03-19T14:28:02.649057627-07:00
- Caption: Tier-3 summary for cluster1 with redshift overlay. Two line styles correspond to the two fixed-redshift runs in that ensemble split; color encodes sigma/m. Top-right panel uses DeltaSigma units of h Msun / pc^2. Annotation box lists stellar-mass range (or N/A for HMF), halo-mass range definition, and concentration scatter used in the ensemble configuration.

### cluster2_ensemble_tier3_summary.png
- Created: 2026-03-19T14:28:03.228564739-07:00
- Caption: Tier-3 summary for cluster2 with redshift overlay. Two line styles correspond to the two fixed-redshift runs in that ensemble split; color encodes sigma/m. Top-right panel uses DeltaSigma units of h Msun / pc^2. Annotation box lists stellar-mass range (or N/A for HMF), halo-mass range definition, and concentration scatter used in the ensemble configuration.

### dwarf1_ensemble_tier3_summary.png
- Created: 2026-03-19T14:28:03.763018131-07:00
- Caption: Tier-3 summary for dwarf1 with redshift overlay. Two line styles correspond to the two fixed-redshift runs in that ensemble split; color encodes sigma/m. Top-right panel uses DeltaSigma units of h Msun / pc^2. Annotation box lists stellar-mass range (or N/A for HMF), halo-mass range definition, and concentration scatter used in the ensemble configuration.

### dwarf2_ensemble_tier3_summary.png
- Created: 2026-03-19T14:28:04.376275063-07:00
- Caption: Tier-3 summary for dwarf2 with redshift overlay. Two line styles correspond to the two fixed-redshift runs in that ensemble split; color encodes sigma/m. Top-right panel uses DeltaSigma units of h Msun / pc^2. Annotation box lists stellar-mass range (or N/A for HMF), halo-mass range definition, and concentration scatter used in the ensemble configuration.

## Update 2026-03-19T14:33:11-07:00

### cluster1_ensemble_tier3_summary.png
- Created: 2026-03-19T14:33:09.188103437-07:00
- Caption: Tier-3 summary for cluster1 with redshift overlay. Two line styles correspond to the two fixed-redshift runs in that ensemble split; color encodes sigma/m. Top-right panel uses DeltaSigma units of h Msun / pc^2. Annotation box lists stellar-mass range (or N/A for HMF), halo-mass range definition, and concentration scatter used in the ensemble configuration.

### cluster2_ensemble_tier3_summary.png
- Created: 2026-03-19T14:33:09.884085417-07:00
- Caption: Tier-3 summary for cluster2 with redshift overlay. Two line styles correspond to the two fixed-redshift runs in that ensemble split; color encodes sigma/m. Top-right panel uses DeltaSigma units of h Msun / pc^2. Annotation box lists stellar-mass range (or N/A for HMF), halo-mass range definition, and concentration scatter used in the ensemble configuration.

### dwarf1_ensemble_tier3_summary.png
- Created: 2026-03-19T14:33:10.445164204-07:00
- Caption: Tier-3 summary for dwarf1 with redshift overlay. Two line styles correspond to the two fixed-redshift runs in that ensemble split; color encodes sigma/m. Top-right panel uses DeltaSigma units of h Msun / pc^2. Annotation box lists stellar-mass range (or N/A for HMF), halo-mass range definition, and concentration scatter used in the ensemble configuration.

### dwarf2_ensemble_tier3_summary.png
- Created: 2026-03-19T14:33:11.074339628-07:00
- Caption: Tier-3 summary for dwarf2 with redshift overlay. Two line styles correspond to the two fixed-redshift runs in that ensemble split; color encodes sigma/m. Top-right panel uses DeltaSigma units of h Msun / pc^2. Annotation box lists stellar-mass range (or N/A for HMF), halo-mass range definition, and concentration scatter used in the ensemble configuration.

## Update 2026-03-19T14:33:40-07:00

### cluster1_ensemble_tier3_summary.png
- Created: 2026-03-19T14:33:38.865943909-07:00
- Caption: Tier-3 summary for cluster1 with redshift overlay. Two line styles correspond to the two fixed-redshift runs in that ensemble split; color encodes sigma/m. Top-right panel uses DeltaSigma units of h Msun / pc^2. Annotation box lists stellar-mass range (or N/A for HMF), halo-mass range definition, and concentration scatter used in the ensemble configuration.

### cluster2_ensemble_tier3_summary.png
- Created: 2026-03-19T14:33:39.590850115-07:00
- Caption: Tier-3 summary for cluster2 with redshift overlay. Two line styles correspond to the two fixed-redshift runs in that ensemble split; color encodes sigma/m. Top-right panel uses DeltaSigma units of h Msun / pc^2. Annotation box lists stellar-mass range (or N/A for HMF), halo-mass range definition, and concentration scatter used in the ensemble configuration.

### dwarf1_ensemble_tier3_summary.png
- Created: 2026-03-19T14:33:40.224849939-07:00
- Caption: Tier-3 summary for dwarf1 with redshift overlay. Two line styles correspond to the two fixed-redshift runs in that ensemble split; color encodes sigma/m. Top-right panel uses DeltaSigma units of h Msun / pc^2. Annotation box lists stellar-mass range (or N/A for HMF), halo-mass range definition, and concentration scatter used in the ensemble configuration.

### dwarf2_ensemble_tier3_summary.png
- Created: 2026-03-19T14:33:40.873161077-07:00
- Caption: Tier-3 summary for dwarf2 with redshift overlay. Two line styles correspond to the two fixed-redshift runs in that ensemble split; color encodes sigma/m. Top-right panel uses DeltaSigma units of h Msun / pc^2. Annotation box lists stellar-mass range (or N/A for HMF), halo-mass range definition, and concentration scatter used in the ensemble configuration.

## Update 2026-03-19T14:35:34-07:00

### cluster1_ensemble_tier3_summary.png
- Created: 2026-03-19T14:35:32.963280439-07:00
- Caption: Tier-3 summary for cluster1 with redshift overlay. Two line styles correspond to the two fixed-redshift runs in that ensemble split; color encodes sigma/m. Top-right panel uses DeltaSigma units of h Msun / pc^2. Annotation box lists stellar-mass range (or N/A for HMF), halo-mass range definition, and concentration scatter used in the ensemble configuration.

### cluster2_ensemble_tier3_summary.png
- Created: 2026-03-19T14:35:33.573145866-07:00
- Caption: Tier-3 summary for cluster2 with redshift overlay. Two line styles correspond to the two fixed-redshift runs in that ensemble split; color encodes sigma/m. Top-right panel uses DeltaSigma units of h Msun / pc^2. Annotation box lists stellar-mass range (or N/A for HMF), halo-mass range definition, and concentration scatter used in the ensemble configuration.

### dwarf1_ensemble_tier3_summary.png
- Created: 2026-03-19T14:35:34.139209986-07:00
- Caption: Tier-3 summary for dwarf1 with redshift overlay. Two line styles correspond to the two fixed-redshift runs in that ensemble split; color encodes sigma/m. Top-right panel uses DeltaSigma units of h Msun / pc^2. Annotation box lists stellar-mass range (or N/A for HMF), halo-mass range definition, and concentration scatter used in the ensemble configuration.

### dwarf2_ensemble_tier3_summary.png
- Created: 2026-03-19T14:35:34.788028955-07:00
- Caption: Tier-3 summary for dwarf2 with redshift overlay. Two line styles correspond to the two fixed-redshift runs in that ensemble split; color encodes sigma/m. Top-right panel uses DeltaSigma units of h Msun / pc^2. Annotation box lists stellar-mass range (or N/A for HMF), halo-mass range definition, and concentration scatter used in the ensemble configuration.
