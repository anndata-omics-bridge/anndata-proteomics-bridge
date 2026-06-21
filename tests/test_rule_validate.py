"""Tests for rules/validate.py (library functions only — CLI tests in test_cli*.py)."""

from __future__ import annotations

from pathlib import Path

from anndata_proteomics.rules.registry import find_rule
from anndata_proteomics.rules.validate import (
    validate_all_packaged,
    validate_file,
)


def test_validate_file_happy() -> None:
    r = validate_file(find_rule("diann", "ion"))
    assert r.ok is True
    assert r.error is None
    assert r.rule is not None
    assert r.rule.software_name == "DIA-NN"


def test_validate_file_bad_returns_result_with_error(tmp_path: Path) -> None:
    bad = tmp_path / "bad.toml"
    bad.write_text("not = valid [[[")
    r = validate_file(bad)
    assert r.ok is False
    assert r.rule is None
    assert r.error is not None and r.error != ""


def test_validate_all_packaged_all_ok() -> None:
    results = validate_all_packaged()
    assert len(results) == 10
    failed = [r for r in results if not r.ok]
    assert not failed, "\n".join(f"{r.path}: {r.error}" for r in failed)
