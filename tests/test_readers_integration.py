"""End-to-end: every packaged TOML reads its corresponding test_data_download file.

Skipped when the test_data_download cache (gitignored, regenerable) is absent.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from anndata_proteomics.readers.dispatch import read_table
from anndata_proteomics.rules.loader import load_rule
from anndata_proteomics.rules.registry import iter_packaged_rules

PROJECT_ROOT = Path(__file__).resolve().parent.parent
TEST_DATA_DIR = PROJECT_ROOT / "test_data_download"
DOWNLOADED_DB = TEST_DATA_DIR / "raw_file_db_downloaded.csv"


def _find_test_data(software_name: str) -> Path | None:
    if not DOWNLOADED_DB.exists():
        return None
    with open(DOWNLOADED_DB) as f:
        for row in csv.DictReader(f):
            if row["software_name"] == software_name and row.get("status") == "ok":
                return TEST_DATA_DIR / "json_dir" / row["input_file_path"]
    return None


@pytest.mark.parametrize(
    "toml_path",
    list(iter_packaged_rules()),
    ids=lambda p: f"{p.parent.name}/{p.name}",
)
def test_reader_loads_test_data_for_packaged_rule(toml_path: Path) -> None:
    rule = load_rule(toml_path)
    data_file = _find_test_data(rule.software_name)
    if data_file is None or not data_file.exists():
        pytest.skip(
            f"no downloaded test data for {rule.software_name!r}; "
            f"regenerate via test_data_download/Makefile"
        )
    df = read_table(data_file)
    expected_min_cols = (
        len(rule.columns.var) + len(rule.columns.obs) + len(rule.layers)
    )
    assert not df.empty, f"{data_file} produced an empty DataFrame"
    assert len(df.columns) >= expected_min_cols, (
        f"{data_file}: got {len(df.columns)} columns, "
        f"rule expects at least {expected_min_cols}"
    )
