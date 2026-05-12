"""Tests for converters/recognize.py."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from anndata_proteomics.converters.recognize import matches, recognize
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
def test_recognize_picks_correct_rule_for_each_vendor(toml_path: Path) -> None:
    rule = load_rule(toml_path)
    data_file = _find_test_data(rule.software_name)
    if data_file is None or not data_file.exists():
        pytest.skip(f"no test data for {rule.software_name!r}")
    headers = list(read_table(data_file).columns)
    recognised = recognize(headers)
    assert recognised is not None, f"no rule matched headers for {rule.software_name}"
    assert recognised.software_name == rule.software_name


def test_matches_long_rule_with_extra_headers_still_matches() -> None:
    # Long rules tolerate extra unrelated columns in the source file.
    diann_rule = load_rule(
        next(p for p in iter_packaged_rules() if p.parent.name == "diann")
    )
    headers = (
        list(diann_rule.columns.obs.select.values())
        + list(diann_rule.columns.var.select.values())
        + [layer.source_column for layer in diann_rule.layers]
        + ["UnrelatedExtraColumn"]
    )
    assert matches(headers, diann_rule) is True


def test_matches_long_rule_returns_false_when_required_column_missing() -> None:
    diann_rule = load_rule(
        next(p for p in iter_packaged_rules() if p.parent.name == "diann")
    )
    headers = (
        list(diann_rule.columns.obs.select.values())
        + list(diann_rule.columns.var.select.values())
        # deliberately drop the layers' source columns
    )
    assert matches(headers, diann_rule) is False


def test_recognize_returns_none_for_empty_headers() -> None:
    assert recognize([]) is None


def test_recognize_returns_none_for_random_headers() -> None:
    assert recognize(["foo", "bar", "baz", "quux"]) is None
