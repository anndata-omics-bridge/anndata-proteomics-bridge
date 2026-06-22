"""Tests for the anndata-proteomics CLI subcommands.

We exercise the subcommand functions directly (calling them as Python)
rather than going through cyclopts' argv parsing — that's a unit-test
shortcut. The dispatch layer is exercised separately by the manual
smoke commands in TODO/PLAN_20260502_jsonschema-and-cli.md §Verification.

CLI output goes through loguru → stderr; tests read `capsys.readouterr().err`.
The `_loguru_to_pytest_capsys` fixture in conftest.py wires loguru into
pytest's stderr capture.
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
    err = capsys.readouterr().err
    assert rc == 0
    assert "0 failed" in err
    assert "PASS" in err


def test_validate_single_path_happy(capsys: pytest.CaptureFixture[str]) -> None:
    path = find_rule("diann", "ion")
    rc = validate(path)
    err = capsys.readouterr().err
    assert rc == 0
    assert "PASS" in err
    assert "0 failed" in err


def test_validate_single_path_bad(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    bad = tmp_path / "bad.toml"
    bad.write_text("not = valid [[[")
    rc = validate(bad)
    err = capsys.readouterr().err
    assert rc == 1
    assert "FAIL" in err
    assert "1 failed" in err


def test_validate_multiple_paths(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    good = find_rule("wombat", "peptidoform")
    bad = tmp_path / "bad.toml"
    bad.write_text("[[")
    rc = validate(good, bad)
    err = capsys.readouterr().err
    assert rc == 1  # any failure → 1
    assert "PASS" in err
    assert "FAIL" in err
    assert "2 rule(s) checked, 1 failed" in err


def test_list_shows_eleven_rules(capsys: pytest.CaptureFixture[str]) -> None:
    rc = list_rules()
    err = capsys.readouterr().err
    assert rc == 0
    lines = [line for line in err.splitlines() if line.strip()]
    assert len(lines) == 11
    assert "diann" in err
    assert "wombat" in err
    assert "peptidoform" in err


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

[axis.duplicates]
mode = "error"

[columns.obs.select]
Run = "Run"

[columns.var.select]
Sequence = "Sequence"
Charge = "Charge"

[[layers]]
name = "Intensity"
source_column = "Intensity"
"""
    )

    output = tmp_path / "out.h5ad"
    rc = convert(data_path, rule_toml=rule_path, output=output)
    err = capsys.readouterr().err
    assert rc == 0
    assert output.exists()
    assert "wrote" in err


def test_convert_returns_one_when_recognition_fails(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    data_path = tmp_path / "unknown.csv"
    data_path.write_text("foo,bar,baz\n1,2,3\n")
    rc = convert(data_path)
    captured = capsys.readouterr()
    assert rc == 1
    assert "auto-recognize" in (captured.out + captured.err).lower()
