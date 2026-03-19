# lessons

## 2026-03-19
- Keep mass-definition handling explicit in every profile/projection function signature and metadata (`M200c`) to prevent silent convention drift.
- Build SIDM integration behind a strict local adapter boundary so upstream `parametricSIDM` API changes fail fast and are isolated.
- Validate the projection kernel directly against NFW analytic behavior before trusting forecast-level metrics.
- For optional external references (`colossus`), catch runtime evaluation failures (including cache-write permission issues) and fall back to deterministic local approximations to preserve pipeline stability.
- When adding Tier-2 outskirts stitching, blend in log-density with configurable match radius and width; linear blending is more prone to visible projection artifacts.
- Auto-updating caption/inventory files should still be followed by explicit, figure-specific descriptive text to keep outputs interpretable without reading code.
- SIDM benchmark grids should be regime-specific; cluster-focused lensing runs should stay in the sub-unity regime while dwarf-focused sensitivity tests can target high effective cross sections, and this split must be reflected consistently in config, code defaults, and docs.
