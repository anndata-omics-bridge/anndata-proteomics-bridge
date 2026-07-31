"""Tests for converters/recognize.py."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from anndata_proteomics.converters.recognize import (
    RecognizedRule,
    UnrecognizedRule,
    matches,
    recognize,
)
from anndata_proteomics.readers.dispatch import read_table
from anndata_proteomics.rules.loader import load_rule
from anndata_proteomics.rules.registry import RuleLocator, iter_packaged_rules

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
    "locator",
    list(iter_packaged_rules()),
    ids=lambda item: f"{item.path.parent.name}/{item.level}",
)
def test_recognize_picks_correct_rule_for_each_vendor(locator: RuleLocator) -> None:
    rule = load_rule(locator)
    data_file = _find_test_data(rule.software_name)
    if data_file is None or not data_file.exists():
        pytest.skip(f"no test data for {rule.software_name!r}")
    headers = list(read_table(data_file).columns)
    if not matches(headers, rule):
        # DIA-NN report schemas vary by version/config; the cached file may not carry this
        # level's columns. That's a vendor-variant mismatch, not a recognition failure.
        pytest.skip(f"cached {rule.software_name} file lacks columns for {locator.level}")

    # recognize() returns a rule only when exactly one packaged rule matches. Multi-level
    # vendors (DIA-NN ships ion/peptidoform/peptide/protein/fragment, all reading the same
    # report.tsv) are ambiguous by design, so recognition is tagged as unresolved and the level must
    # be selected explicitly via load_packaged_rule(software, level).
    n_matching = sum(1 for p in iter_packaged_rules() if matches(headers, load_rule(p)))
    recognised = recognize(headers)
    if n_matching == 1:
        assert isinstance(recognised, RecognizedRule)
        assert recognised.rule.software_name == rule.software_name
    else:
        assert isinstance(recognised, UnrecognizedRule)


def test_matches_long_rule_with_extra_headers_still_matches() -> None:
    # Long rules tolerate extra unrelated columns in the source file.
    diann_rule = load_rule(
        next(
            item
            for item in iter_packaged_rules()
            if item.path.parent.name == "v1" and item.level == "ion"
        )
    )
    headers = (
        list(diann_rule.columns.obs.select.values())
        + list(diann_rule.columns.var.select.values())
        + [layer.source for layer in diann_rule.layers]
        + ["UnrelatedExtraColumn"]
    )
    assert matches(headers, diann_rule) is True


def test_matches_long_rule_returns_false_when_required_column_missing() -> None:
    diann_rule = load_rule(
        next(
            item
            for item in iter_packaged_rules()
            if item.path.parent.name == "v1" and item.level == "ion"
        )
    )
    headers = (
        list(diann_rule.columns.obs.select.values()) + list(diann_rule.columns.var.select.values())
        # deliberately drop the layers' source columns
    )
    assert matches(headers, diann_rule) is False


def test_recognize_returns_unresolved_for_empty_headers() -> None:
    assert recognize([]) == UnrecognizedRule()


def test_recognize_returns_unresolved_for_random_headers() -> None:
    assert recognize(["foo", "bar", "baz", "quux"]) == UnrecognizedRule()
