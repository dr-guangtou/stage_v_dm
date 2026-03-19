"""I/O helpers for benchmark outputs."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def ensure_output_directories(base_directory: Path) -> dict[str, Path]:
    """Create output directory structure and return key paths."""
    figures_directory = base_directory / "figures"
    tables_directory = base_directory / "tables"
    intermediate_directory = base_directory / "intermediate"

    figures_directory.mkdir(parents=True, exist_ok=True)
    tables_directory.mkdir(parents=True, exist_ok=True)
    intermediate_directory.mkdir(parents=True, exist_ok=True)

    return {
        "figures": figures_directory,
        "tables": tables_directory,
        "intermediate": intermediate_directory,
    }


def save_table(data_frame: pd.DataFrame, output_path: Path) -> None:
    """Write a machine-readable CSV table with deterministic ordering."""
    data_frame.to_csv(output_path, index=False)
