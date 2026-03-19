"""Tests for target-significance precision sweep utilities."""

from sidm_stagev_forecast.config import CLUSTER_BENCHMARK, DWARF_BENCHMARK
from sidm_stagev_forecast.precision_sweep import default_radial_windows, run_precision_sweep


def test_default_radial_windows_are_ordered_and_cover_range() -> None:
    windows = default_radial_windows(DWARF_BENCHMARK)
    assert [window.name for window in windows] == ["full", "inner", "outer"]

    full, inner, outer = windows
    assert full.r_min_kpc == DWARF_BENCHMARK.r_lensing_min_kpc
    assert full.r_max_kpc == DWARF_BENCHMARK.r_lensing_max_kpc
    assert inner.r_min_kpc == full.r_min_kpc
    assert outer.r_max_kpc == full.r_max_kpc
    assert inner.r_max_kpc == outer.r_min_kpc


def test_run_precision_sweep_output_schema() -> None:
    table = run_precision_sweep(sidm_backend="surrogate")

    required_columns = {
        "benchmark",
        "m200_msun",
        "sigma_over_m_cm2_g",
        "window",
        "window_r_min_kpc",
        "window_r_max_kpc",
        "n_bins",
        "target_sigma",
        "required_uniform_fraction",
        "required_uniform_percent",
        "max_abs_ratio_shift_percent",
    }
    assert required_columns.issubset(set(table.columns))
    assert len(table) > 0

    assert set(table["benchmark"].unique()) == {
        DWARF_BENCHMARK.label,
        CLUSTER_BENCHMARK.label,
    }
    assert set(table["target_sigma"].unique()) == {2.0, 3.0, 5.0}
    assert (table["required_uniform_percent"] >= 0.0).all()
