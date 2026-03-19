"""I/O helpers for benchmark outputs."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

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


def append_figure_caption_entries(
    caption_path: Path,
    entries: list[dict[str, str]],
    timezone_name: str = "America/Phoenix",
) -> None:
    """Append figure captions with creation timestamp metadata."""
    timestamp = datetime.now(ZoneInfo(timezone_name)).isoformat(timespec="seconds")
    with caption_path.open("a", encoding="utf-8") as handle:
        handle.write(f"\n## Update {timestamp}\n")
        for entry in entries:
            handle.write(f"\n### {entry['filename']}\n")
            handle.write(f"- Created: {entry['created_at']}\n")
            handle.write(f"- Caption: {entry['caption']}\n")


def append_inventory_entries(
    inventory_path: Path,
    entries: list[dict[str, str]],
    timezone_name: str = "America/Phoenix",
) -> None:
    """Append concise inventory descriptions for outputs tables/intermediate files."""
    timestamp = datetime.now(ZoneInfo(timezone_name)).isoformat(timespec="seconds")
    with inventory_path.open("a", encoding="utf-8") as handle:
        handle.write(f"\n## Update {timestamp}\n")
        for entry in entries:
            handle.write(f"- {entry['path']}: {entry['description']}\n")
