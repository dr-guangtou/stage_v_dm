"""
Generate all Phase 1 figures: SHMR model setup and validation.

Produces:
    outputs/phase1/shmr_validation.png      — M*/Mh vs Mh at z=0, 0.5, 1.0, 2.0
    outputs/phase1/scatter_vs_log_Mh.png    — Mass-dependent scatter model
    outputs/phase1/shmr_with_scatter.png    — SHMR + scatter at three redshifts
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
matplotlib.use("Agg")

from shmr_fisher.config import SHMRParams
from shmr_fisher.plot_results import plot_shmr_validation
from scripts.plot_scatter_model import plot_scatter_vs_log_Mh
from scripts.plot_shmr_with_scatter import plot_shmr_with_scatter

outdir = Path("outputs") / "phase1"
outdir.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("Phase 1 figures")
print("=" * 60)

# 1. SHMR validation (M*/Mh vs Mh)
print("\n--- shmr_validation.png ---")
params = SHMRParams(use_mass_dependent_scatter=True)
plot_shmr_validation(params, save_path=outdir / "shmr_validation.png")

# 2. Scatter model
print("\n--- scatter_vs_log_Mh.png ---")
plot_scatter_vs_log_Mh()

# 3. SHMR with scatter envelopes
print("\n--- shmr_with_scatter.png ---")
plot_shmr_with_scatter()

print("\nAll Phase 1 figures generated.")
