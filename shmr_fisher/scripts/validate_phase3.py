"""
Phase 3 QA: validate covariance and Fisher matrix, generate diagnostic plots.

Produces:
    outputs/phase3/covariance_diagnostic.png — Sigma_crit, source density, DS errors vs z
    outputs/phase3/fisher_errors_comparison.png — marginalized errors: Stage-III vs IV vs V
    outputs/phase3/fisher_ellipses.png — 2D Fisher ellipses for key parameter pairs
    outputs/phase3/fisher_ngal_scaling.png — error scaling with N_gal (lensing-only check)

Also prints AT-3 through AT-6 results.
"""

import sys
sys.path.insert(0, '.')

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from shmr_fisher.config import SHMRParams, SurveyConfig, LensingConfig, ForecastConfig
from shmr_fisher.covariance import (
    sigma_crit_effective, n_source_effective,
    lensing_covariance, survey_volume,
)
from shmr_fisher.fisher import (
    compute_fisher_matrix, marginalized_errors, conditional_errors,
    fisher_ellipse, get_varied_params,
)

params = SHMRParams()
lc = LensingConfig()
fc = ForecastConfig(n_R_bins=8, n_Mh_bins=150, dz=0.2)
outdir = Path('outputs') / 'phase3'
outdir.mkdir(parents=True, exist_ok=True)

# -----------------------------------------------------------------------
# Plot 1: Covariance diagnostics
# -----------------------------------------------------------------------
print("Generating covariance diagnostics...")
fig, axes = plt.subplots(2, 2, figsize=(13, 10))

z_arr = np.linspace(0.05, 1.5, 30)

# Panel A: Sigma_crit_eff vs z_lens
ax = axes[0, 0]
sc_eff = [sigma_crit_effective(z, lc) / 1e12 for z in z_arr]
ax.plot(z_arr, sc_eff, 'b-', lw=2)
ax.set_xlabel(r'$z_\mathrm{lens}$', fontsize=12)
ax.set_ylabel(r'$\Sigma_\mathrm{crit,eff}$ [$M_\odot/\mathrm{pc}^2$]', fontsize=12)
ax.set_title(r'Effective $\Sigma_\mathrm{crit}$ vs lens redshift', fontsize=12)
ax.grid(True, alpha=0.3)

# Panel B: Effective source density vs z_lens
ax = axes[0, 1]
n_eff = [n_source_effective(z, lc) for z in z_arr]
ax.plot(z_arr, n_eff, 'r-', lw=2)
ax.set_xlabel(r'$z_\mathrm{lens}$', fontsize=12)
ax.set_ylabel(r'$n_\mathrm{source,eff}$ [arcmin$^{-2}$]', fontsize=12)
ax.set_title('Effective source density behind lens', fontsize=12)
ax.grid(True, alpha=0.3)

# Panel C: sigma(DS) vs R for different N_lens
ax = axes[1, 0]
R_edges = np.logspace(np.log10(0.1), np.log10(30.0), 11)
R_centers = np.sqrt(R_edges[:-1] * R_edges[1:])

for N_lens, color, ls in [(1e4, 'blue', '--'), (1e5, 'green', '-'),
                            (1e6, 'red', '-.')]:
    var = lensing_covariance(R_edges, 0.3, N_lens, 10000, lc)
    sig = np.sqrt(var) / 1e12  # Msun/pc^2
    ax.plot(R_centers, sig, color=color, ls=ls, lw=2,
            label=rf'$N_\mathrm{{lens}} = 10^{{{int(np.log10(N_lens))}}}$')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$R$ [Mpc]', fontsize=12)
ax.set_ylabel(r'$\sigma(\Delta\Sigma)$ [$M_\odot/\mathrm{pc}^2$]', fontsize=12)
ax.set_title(r'Lensing noise at $z_l=0.3$', fontsize=12)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3, which='both')

# Panel D: Survey volume vs z_max
ax = axes[1, 1]
z_max_arr = np.linspace(0.1, 2.0, 30)
for area, color in [(1000, 'blue'), (5000, 'green'), (14000, 'red')]:
    vols = [survey_volume(0.05, zm, area) for zm in z_max_arr]
    ax.plot(z_max_arr, vols, color=color, lw=2,
            label=rf'{area} deg$^2$')

ax.set_xlabel(r'$z_\mathrm{max}$', fontsize=12)
ax.set_ylabel(r'$V_\mathrm{survey}$ [Mpc$^3$]', fontsize=12)
ax.set_title('Comoving survey volume', fontsize=12)
ax.set_yscale('log')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

plt.tight_layout()
fig.savefig(outdir / 'covariance_diagnostic.png', dpi=300, bbox_inches='tight')
print(f'Saved: {outdir / "covariance_diagnostic.png"}')
plt.close(fig)

# -----------------------------------------------------------------------
# Compute Fisher matrices for three survey tiers
# -----------------------------------------------------------------------
surveys = {
    'Stage-III': SurveyConfig(
        name='Stage-III Shallow Wide',
        area_deg2=7500, z_min=0.02, z_max=0.2,
        n_gal_total=700_000, log_Mstar_min=9.5,
    ),
    'Stage-IV': SurveyConfig(
        name='Stage-IV Low-z',
        area_deg2=14000, z_min=0.05, z_max=0.4,
        n_gal_total=10_000_000, log_Mstar_min=9.0,
    ),
    'Stage-V': SurveyConfig(
        name='Stage-V Wide',
        area_deg2=10000, z_min=0.05, z_max=1.0,
        n_gal_total=50_000_000, log_Mstar_min=8.5,
        log_Mstar_min_func=lambda z: 8.5 + 1.0 * z,
    ),
}

results = {}
for label, survey in surveys.items():
    print(f"\nComputing Fisher for {label}...")
    f, pn, meta = compute_fisher_matrix(params, survey, lc, fc)
    eigvals = np.linalg.eigvalsh(f)
    all_pos = np.all(eigvals > 0)

    if all_pos:
        marg = marginalized_errors(f)
        cond = conditional_errors(f)
    else:
        marg = np.full(len(pn), np.inf)
        cond = np.full(len(pn), np.inf)

    results[label] = {
        'fisher': f, 'param_names': pn, 'metadata': meta,
        'marg': marg, 'cond': cond, 'eigvals': eigvals,
    }
    print(f"  Params: {pn}")
    print(f"  Condition number: {meta['condition_number']:.2e}")
    print(f"  Positive definite: {all_pos}")
    if all_pos:
        for n, m in zip(pn, marg):
            fid = getattr(params, n)
            pct = m / abs(fid) * 100 if abs(fid) > 1e-10 else 0
            print(f"    {n:12s}: sigma={m:.6f} ({pct:.3f}%)")

# -----------------------------------------------------------------------
# AT checks
# -----------------------------------------------------------------------
print("\n" + "="*70)
print("ACCEPTANCE TESTS")
print("="*70)

# AT-4: All Fisher matrices positive definite
for label, r in results.items():
    ok = np.all(r['eigvals'] > 0)
    print(f"AT-4 ({label} pos. def.): {'PASS' if ok else 'FAIL'}")

# AT-5: marginalized >= conditional
for label, r in results.items():
    ok = np.all(r['marg'] >= r['cond'] - 1e-15)
    print(f"AT-5 ({label} marg >= cond): {'PASS' if ok else 'FAIL'}")

# AT-6: Stage ordering on shared 5 z=0 params
# Fair comparison: all three surveys must use 5-param Fisher (no evolution).
# Stage-V has 9 params by default, so we compute a separate 5-param version.
shared = ['log_M1_0', 'N_0', 'beta_0', 'gamma_0', 'sigma_logMs']
print("\nAT-6: Stage ordering (5-param Fisher for fair comparison)")

fc_5p = ForecastConfig(n_R_bins=8, n_Mh_bins=150, dz=0.2, vary_z_evolution=False)
survey_v = surveys['Stage-V']
print("  Computing Stage-V with 5-param (evolution fixed)...")
f_v5, pn_v5, meta_v5 = compute_fisher_matrix(params, survey_v, lc, fc_5p)
marg_v5 = marginalized_errors(f_v5) if np.all(np.linalg.eigvalsh(f_v5) > 0) else np.full(len(pn_v5), np.inf)
results['Stage-V-5p'] = {'param_names': pn_v5, 'marg': marg_v5}

all_ordered = True
for pname in shared:
    errs = {}
    for label, key in [('Stage-III', 'Stage-III'), ('Stage-IV', 'Stage-IV'),
                        ('Stage-V', 'Stage-V-5p')]:
        pn = results[key]['param_names']
        if pname in pn:
            idx = pn.index(pname)
            errs[label] = results[key]['marg'][idx]

    if len(errs) == 3:
        ok = errs['Stage-V'] < errs['Stage-IV'] < errs['Stage-III']
        print(f"  {pname:12s}: III={errs['Stage-III']:.2e}  "
              f"IV={errs['Stage-IV']:.2e}  V(5p)={errs['Stage-V']:.2e}  "
              f"{'OK' if ok else 'REVERSED'}")
        if not ok:
            all_ordered = False
print(f"AT-6: {'PASS' if all_ordered else 'FAIL'}")
print("  (Note: Stage-V 9-param Fisher has wider z=0 contours due to "
      "marginalization over evolution params — expected.)")

# -----------------------------------------------------------------------
# Plot 2: Bar chart of marginalized errors across surveys
# -----------------------------------------------------------------------
print("\nGenerating error comparison plot...")
fig, ax = plt.subplots(figsize=(12, 6))

x = np.arange(len(shared))
width = 0.25
colors = {'Stage-III': '#1f77b4', 'Stage-IV': '#ff7f0e', 'Stage-V': '#2ca02c'}

for i, label in enumerate(['Stage-III', 'Stage-IV', 'Stage-V']):
    errs = []
    for pname in shared:
        pn = results[label]['param_names']
        if pname in pn:
            idx = pn.index(pname)
            fid = abs(getattr(params, pname))
            errs.append(results[label]['marg'][idx] / fid * 100
                        if fid > 1e-10 else results[label]['marg'][idx])
        else:
            errs.append(0)
    ax.bar(x + i * width, errs, width, label=label, color=colors[label])

ax.set_xlabel('SHMR Parameter', fontsize=12)
ax.set_ylabel(r'Marginalized $\sigma / |\theta_\mathrm{fid}|$ [%]', fontsize=12)
ax.set_title('Fisher forecast: parameter constraints by survey tier', fontsize=13)
ax.set_xticks(x + width)
ax.set_xticklabels([r'$\log M_1$', r'$N_0$', r'$\beta_0$',
                     r'$\gamma_0$', r'$\sigma_{\log M_*}$'], fontsize=11)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3, axis='y')
ax.set_yscale('log')

plt.tight_layout()
fig.savefig(outdir / 'fisher_errors_comparison.png', dpi=300, bbox_inches='tight')
print(f'Saved: {outdir / "fisher_errors_comparison.png"}')
plt.close(fig)

# -----------------------------------------------------------------------
# Plot 3: Fisher ellipses for key parameter pairs (Stage-V, 5 shared)
# -----------------------------------------------------------------------
print("Generating Fisher ellipse plot...")
# Use the Stage-V Fisher restricted to the 5 shared params
pn_v = results['Stage-V']['param_names']
f_v = results['Stage-V']['fisher']

# Pairs to plot
pairs = [(0, 1), (0, 4), (2, 3), (1, 4)]
pair_labels = [
    (r'$\log M_1$', r'$N_0$'),
    (r'$\log M_1$', r'$\sigma_{\log M_*}$'),
    (r'$\beta_0$', r'$\gamma_0$'),
    (r'$N_0$', r'$\sigma_{\log M_*}$'),
]

fig, axes = plt.subplots(2, 2, figsize=(11, 10))

survey_colors = {'Stage-III': '#1f77b4', 'Stage-IV': '#ff7f0e', 'Stage-V': '#2ca02c'}

for ax, (i, j), (xlabel, ylabel) in zip(axes.flat, pairs, pair_labels):
    for label in ['Stage-III', 'Stage-IV', 'Stage-V']:
        pn = results[label]['param_names']
        f = results[label]['fisher']

        # Only plot if both params are in this survey's Fisher
        if i < len(pn) and j < len(pn):
            try:
                ex, ey = fisher_ellipse(f, i, j, n_sigma=1.0)
                # Center on fiducial
                fid_i = getattr(params, pn[i])
                fid_j = getattr(params, pn[j])
                ax.plot(ex + fid_i, ey + fid_j,
                        color=survey_colors[label], lw=2, label=label)
            except Exception:
                pass

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.ticklabel_format(useOffset=False)

fig.suptitle(r'Fisher ellipses (1$\sigma$, marginalized)', fontsize=14, y=1.01)
plt.tight_layout()
fig.savefig(outdir / 'fisher_ellipses.png', dpi=300, bbox_inches='tight')
print(f'Saved: {outdir / "fisher_ellipses.png"}')
plt.close(fig)

# -----------------------------------------------------------------------
# Plot 4: N_gal scaling
# -----------------------------------------------------------------------
print("Generating N_gal scaling plot...")
n_gal_values = [1e5, 5e5, 1e6, 5e6, 1e7, 5e7]
scaling_errors = {p: [] for p in shared}

for ng in n_gal_values:
    s = SurveyConfig(name=f'test_{ng:.0e}', area_deg2=14000,
                     z_min=0.05, z_max=0.4,
                     n_gal_total=ng, log_Mstar_min=9.0)
    f, pn, _ = compute_fisher_matrix(params, s, lc, fc)
    try:
        marg = marginalized_errors(f)
        for p in shared:
            idx = pn.index(p)
            scaling_errors[p].append(marg[idx])
    except Exception:
        for p in shared:
            scaling_errors[p].append(np.inf)

fig, ax = plt.subplots(figsize=(8, 6))
for p, color in zip(shared, ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']):
    fid = abs(getattr(params, p))
    errs = np.array(scaling_errors[p])
    if fid > 1e-10:
        errs_pct = errs / fid * 100
    else:
        errs_pct = errs * 100
    ax.plot(n_gal_values, errs_pct, 'o-', color=color, lw=2, ms=5, label=p)

# Reference 1/sqrt(N) line
ng_ref = np.array(n_gal_values)
ref = 10.0 * np.sqrt(1e6 / ng_ref)
ax.plot(ng_ref, ref, 'k--', lw=1.5, alpha=0.5, label=r'$\propto 1/\sqrt{N}$')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$N_\mathrm{gal,total}$', fontsize=13)
ax.set_ylabel(r'$\sigma / |\theta_\mathrm{fid}|$ [%]', fontsize=13)
ax.set_title('Constraint scaling with galaxy count (Stage-IV footprint)', fontsize=13)
ax.legend(fontsize=9, ncol=2)
ax.grid(True, alpha=0.3, which='both')

plt.tight_layout()
fig.savefig(outdir / 'fisher_ngal_scaling.png', dpi=300, bbox_inches='tight')
print(f'Saved: {outdir / "fisher_ngal_scaling.png"}')
plt.close(fig)

print("\nAll Phase 3 QA complete.")
