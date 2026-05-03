"""Tests for the anndata-proteomics CLI subcommands.

We exercise the subcommand functions directly (calling them as Python)
rather than going through cyclopts' argv parsing — that's a unit-test
shortcut. The dispatch layer is exercised separately by the manual
smoke commands in TODO/PLAN_20260502_jsonschema-and-cli.md §Verification.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from anndata_proteomics.rules.registry import find_rule, packaged_rules_root
from anndata_proteomics.scripts.cli import (
    convert,
    export_schema_cmd,
    list_rules,
    validate,
)


def test_validate_no_args_walks_packaged(capsys: pytest.CaptureFixture[str]) -> None:
    rc = validate()
    out = capsys.readouterr().out
    assert rc == 0
    assert "0 failed" in out
    assert "PASS" in out


def test_validate_single_path_happy(capsys: pytest.CaptureFixture[str]) -> None:
    path = find_rule("diann", "ion")
    rc = validate(path)
    out = capsys.readouterr().out
    assert rc == 0
    assert "PASS" in out
    assert "0 failed" in out


def test_validate_single_path_bad(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    bad = tmp_path / "bad.toml"
    bad.write_text("not = valid [[[")
    rc = validate(bad)
    out = capsys.readouterr().out
    assert rc == 1
    assert "FAIL" in out
    assert "1 failed" in out


def test_validate_multiple_paths(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    good = find_rule("wombat", "peptidoform")
    bad = tmp_path / "bad.toml"
    bad.write_text("[[")
    rc = validate(good, bad)
    out = capsys.readouterr().out
    assert rc == 1  # any failure → 1
    assert "PASS" in out
    assert "FAIL" in out
    assert "2 rule(s) checked, 1 failed" in out


def test_list_shows_six_rules(capsys: pytest.CaptureFixture[str]) -> None:
    rc = list_rules()
    out = capsys.readouterr().out
    assert rc == 0
    lines = [line for line in out.splitlines() if line.strip()]
    assert len(lines) == 6
    assert "diann" in out
    assert "wombat" in out
    assert "peptidoform" in out


def test_export_schema_writes_file() -> None:
    rc = export_schema_cmd()
    assert rc == 0
    schema_path = packaged_rules_root() / "_schema" / "parse_rule.schema.json"
    assert schema_path.exists()
    assert schema_path.stat().st_size > 100


def test_convert_with_explicit_rule_toml_writes_h5ad(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # Synthesise a tiny long DataFrame matching a stripped-down rule.
    import pandas as pd

    data_path = tmp_path / "tiny.tsv"
    pd.DataFrame(
        {
            "Run": ["S1", "S2"],
            "Sequence": ["P1", "P1"],
            "Charge": [2, 2],
            "Intensity": [10.0, 20.0],
        }
    ).to_csv(data_path, sep="\t", index=False)

    rule_path = tmp_path / "rule.toml"
    rule_path.write_text(
        """
schema_version = "0.1"
file_version = "1"
software_name = "Tiny"
input_shape = "long"
quantification_level = "ion"

[axis]
obs_keys = ["Run"]
var_keys = ["Sequence", "Charge"]
x_layer = "Intensity"

[columns.obs]
Run = "Run"

[columns.var]
Sequence = "Sequence"
Charge = "Charge"

[[layers]]
name = "Intensity"
source_column = "Intensity"

[duplicates]
mode = "error"
"""
    )

    output = tmp_path / "out.h5ad"
    rc = convert(data_path, rule_toml=rule_path, output=output)
    out = capsys.readouterr().out
    assert rc == 0
    assert output.exists()
    assert "wrote" in out


def test_convert_returns_one_when_recognition_fails(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    data_path = tmp_path / "unknown.csv"
    data_path.write_text("foo,bar,baz\n1,2,3\n")
    rc = convert(data_path)
    captured = capsys.readouterr()
    assert rc == 1
    assert "auto-recognize" in (captured.out + captured.err).lower()
