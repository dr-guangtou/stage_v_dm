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
