"""
Validation suite for the SHMR forecast pipeline.

Run after each development phase to verify correctness. Each validation
function prints diagnostic output and returns True if all checks pass.

Phase 1 validations:
- validate_shmr: SHMR shape, round-trip, peak efficiency at z=0 and z=1.
"""

from __future__ import annotations

import numpy as np

from .config import SHMRParams
from .shmr_model import mean_log_Mh, mean_log_Mstar, shmr_params_at_z


def validate_shmr(params: SHMRParams | None = None) -> bool:
    """
    Validate the SHMR model at z=0 and z=1.

    Checks
    ------
    1. AT-1 round-trip: mean_log_Mh(mean_log_Mstar(12.0)) == 12.0 within 0.01 dex.
    2. Peak M*/Mh at z=0 is between 0.02 and 0.06.
    3. Peak occurs at log Mh between 11.5 and 12.5 at z=0.
    4. mean_log_Mstar(12.0, z=0) is between 10.4 and 10.8.

    Parameters
    ----------
    params : SHMRParams or None
        If None, uses fiducial defaults.

    Returns
    -------
    passed : bool
    """
    if params is None:
        params = SHMRParams()

    all_passed = True

    # -- AT-1: Round-trip test --
    log_Mh_input = 12.0
    log_Mstar_fwd = mean_log_Mstar(log_Mh_input, params, z=0)
    log_Mh_roundtrip = mean_log_Mh(log_Mstar_fwd, params, z=0)
    roundtrip_error = abs(log_Mh_roundtrip - log_Mh_input)

    print(f"AT-1 Round-trip test:")
    print(f"  Input:  log Mh = {log_Mh_input:.2f}")
    print(f"  Forward: log M* = {log_Mstar_fwd:.4f}")
    print(f"  Inverse: log Mh = {log_Mh_roundtrip:.6f}")
    print(f"  Error:   {roundtrip_error:.6f} dex (threshold: 0.01)")

    if roundtrip_error > 0.01:
        print("  FAIL: round-trip error exceeds 0.01 dex")
        all_passed = False
    else:
        print("  PASS")

    # -- Check mean_log_Mstar(12.0, z=0) in expected range --
    print(f"\nmean_log_Mstar(12.0, z=0) = {log_Mstar_fwd:.4f}")
    if 10.4 <= log_Mstar_fwd <= 10.8:
        print("  PASS (expected 10.4-10.8)")
    else:
        print(f"  FAIL: expected 10.4-10.8, got {log_Mstar_fwd:.4f}")
        all_passed = False

    # -- Peak M*/Mh at z=0 --
    log_Mh_grid = np.linspace(10.0, 15.0, 500)
    log_Mstar_grid = mean_log_Mstar(log_Mh_grid, params, z=0)

    # M*/Mh = 10^(log_Mstar - log_Mh)
    ratio = 10.0 ** (log_Mstar_grid - log_Mh_grid)
    peak_idx = np.argmax(ratio)
    peak_ratio = ratio[peak_idx]
    peak_log_Mh = log_Mh_grid[peak_idx]

    print(f"\nPeak M*/Mh at z=0:")
    print(f"  Peak ratio = {peak_ratio:.4f} (expected 0.02-0.06)")
    print(f"  at log Mh  = {peak_log_Mh:.2f} (expected 11.5-12.5)")

    if 0.02 <= peak_ratio <= 0.06:
        print("  Peak ratio: PASS")
    else:
        print(f"  Peak ratio: FAIL")
        all_passed = False

    if 11.5 <= peak_log_Mh <= 12.5:
        print("  Peak location: PASS")
    else:
        print(f"  Peak location: FAIL")
        all_passed = False

    # -- Check redshift evolution: peak should shift to higher Mh at z=1 --
    log_Mstar_z1 = mean_log_Mstar(log_Mh_grid, params, z=1.0)
    ratio_z1 = 10.0 ** (log_Mstar_z1 - log_Mh_grid)
    peak_idx_z1 = np.argmax(ratio_z1)
    peak_log_Mh_z1 = log_Mh_grid[peak_idx_z1]

    print(f"\nRedshift evolution check:")
    print(f"  Peak log Mh at z=0: {peak_log_Mh:.2f}")
    print(f"  Peak log Mh at z=1: {peak_log_Mh_z1:.2f}")
    if peak_log_Mh_z1 >= peak_log_Mh:
        print("  PASS (peak shifts to higher Mh at z=1)")
    else:
        print("  INFO: peak shifted to lower Mh at z=1 (may be OK depending on params)")

    print(f"\n{'='*50}")
    print(f"SHMR validation: {'ALL PASSED' if all_passed else 'SOME CHECKS FAILED'}")
    print(f"{'='*50}")

    return all_passed
