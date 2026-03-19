# lessons

## 2026-03-19
- Keep mass-definition handling explicit in every profile/projection function signature and metadata (`M200c`) to prevent silent convention drift.
- Build SIDM integration behind a strict local adapter boundary so upstream `parametricSIDM` API changes fail fast and are isolated.
- Validate the projection kernel directly against NFW analytic behavior before trusting forecast-level metrics.
