"""Tests for rules/validate.py."""

from __future__ import annotations

from pathlib import Path

import pytest

from anndata_proteomics.rules import validate as validate_mod
from anndata_proteomics.rules.registry import find_rule
from anndata_proteomics.rules.validate import (
    main,
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
    assert len(results) == 6
    failed = [r for r in results if not r.ok]
    assert not failed, "\n".join(f"{r.path}: {r.error}" for r in failed)


def test_main_returns_zero_when_all_pass(capsys: pytest.CaptureFixture[str]) -> None:
    rc = main([])
    captured = capsys.readouterr()
    assert rc == 0
    assert "0 failed" in captured.out
    assert "PASS" in captured.out


def test_main_returns_one_on_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    bad = tmp_path / "bad.toml"
    bad.write_text("[[")
    monkeypatch.setattr(validate_mod, "iter_packaged_rules", lambda: iter([bad]))
    rc = main([])
    captured = capsys.readouterr()
    assert rc == 1
    assert "1 failed" in captured.out
    assert "FAIL" in captured.out
