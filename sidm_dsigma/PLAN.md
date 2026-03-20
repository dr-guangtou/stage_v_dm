# PLAN.md

## Current Development State

Core implementation milestones are complete:

1. Baseline profile and projection pipeline.
2. Tier-1 ensemble stacking (HMF + SHMR modes).
3. Tier-2 hybrid inner+outer profile path.
4. Tier-3 empirical outer-correction path.
5. Split ensemble redshift-overlay Tier-3 summaries.

## Phase-1 Organization (Current Task)

### Goal

Bring repository documentation and figure outputs into alignment with the current codebase.

### Tasks

1. Align docs with current implementation:
   - `README.md`
   - `SPEC.md`
   - `PLAN.md`
   - `docs/todo.md`
   - `docs/lessons.md`
2. Prune `outputs/figures` to retain only the latest Tier-3 overlay summaries with mass/concentration inset distributions.
2. Prune `outputs/figures` to retain only default summary figures for cluster and dwarf ensemble runs.
3. Rewrite `outputs/figures/CAPTION.md` to match retained figures only.
4. Verify clean, merge-ready state in this worktree.

### Acceptance Criteria

- Documentation accurately reflects Tier-1/2/3 architecture and active scripts/config pathways.
- `outputs/figures` contains only:
  - `cluster_ensemble_summary.png`
  - `dwarf_ensemble_summary.png`
  - `CAPTION.md`
- `CAPTION.md` contains unambiguous captions and timestamps for exactly those retained figures.

## Near-Term Next Phase (After Organization)

1. Add a `summary-only` forecast mode to avoid regenerating non-summary figure products.
2. Separate heavy debug figures from canonical deliverables via dedicated output profiles.
3. Add a small CI check to prevent stale caption entries and orphaned figure files.

## Non-Goals for This Organization Phase

- No new physics-model changes.
- No rerun of full heavy forecast production unless required for integrity checks.
- No expansion of Tier-3 model family in this phase.
