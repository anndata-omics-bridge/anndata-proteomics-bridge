"""Tests for the anndata-proteomics CLI subcommands.

We exercise the subcommand functions directly (calling them as Python)
rather than going through cyclopts' argv parsing — that's a unit-test
shortcut. The dispatch layer is exercised separately by the manual
smoke commands in ../TODO/Archive/PLAN_20260502_jsonschema-and-cli.md §Verification.

CLI output goes through loguru → stderr; tests read `capsys.readouterr().err`.
The `_loguru_to_pytest_capsys` fixture in conftest.py wires loguru into
pytest's stderr capture.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from anndata_proteomics.rules.registry import find_rule_for_version, packaged_rules_root
from anndata_proteomics.scripts.cli import (
    ConvertCliOptions,
    _write_atomically,
    convert,
    export_schema_cmd,
    list_rules,
    validate,
)


def test_atomic_write_preserves_existing_output_after_failure(tmp_path: Path) -> None:
    output = tmp_path / "result.h5mu"
    output.write_text("complete")

    def fail_after_partial_write(path: Path) -> None:
        path.write_text("partial")
        raise RuntimeError("write failed")

    with pytest.raises(RuntimeError, match="write failed"):
        _write_atomically(output, fail_after_partial_write)

    assert output.read_text() == "complete"
    assert list(tmp_path.iterdir()) == [output]


def test_atomic_write_replaces_output_after_success(tmp_path: Path) -> None:
    output = tmp_path / "result.h5ad"
    output.write_text("old")

    def write_new(path: Path) -> None:
        path.write_text("new")

    _write_atomically(output, write_new)

    assert output.read_text() == "new"
    assert list(tmp_path.iterdir()) == [output]


def _tiny_rule_document(*, include_protein: bool = False) -> dict[str, object]:
    """Return a minimal self-contained long-format document for CLI tests."""
    levels: dict[str, object] = {
        "ion": {
            "axis": {
                "var_keys": ["Sequence", "Charge"],
                "x_layer": "Intensity",
            },
            "columns": {
                "var": {
                    "select": {
                        "Sequence": "Sequence",
                        "Charge": "Charge",
                    }
                }
            },
            "layers": [{"name": "Intensity", "source": "Intensity"}],
        }
    }
    if include_protein:
        levels["protein"] = {
            "axis": {"var_keys": ["Sequence"], "x_layer": "Intensity"},
            "columns": {"var": {"select": {"Sequence": "Sequence"}}},
            "layers": [{"name": "Intensity", "source": "Intensity"}],
        }
    return {
        "schema_version": "0.1",
        "file_version": "1",
        "software_name": "Tiny",
        "software_version": "1.0",
        "base": {
            "input_shape": "long",
            "axis": {"obs_keys": ["Run"], "duplicates": {"mode": "error"}},
            "columns": {"obs": {"select": {"Run": "Run"}}},
        },
        "levels": levels,
    }


def test_validate_no_args_walks_packaged(capsys: pytest.CaptureFixture[str]) -> None:
    rc = validate()
    err = capsys.readouterr().err
    assert rc == 0
    assert "0 failed" in err
    assert "PASS" in err


def test_validate_single_path_happy(capsys: pytest.CaptureFixture[str]) -> None:
    path = find_rule_for_version("diann", "ion", "2.0.0").path
    rc = validate(path)
    err = capsys.readouterr().err
    assert rc == 0
    assert "PASS" in err
    assert "0 failed" in err


def test_validate_single_path_bad(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    bad = tmp_path / "bad.json"
    bad.write_text('{"not": valid}')
    rc = validate(bad)
    err = capsys.readouterr().err
    assert rc == 1
    assert "FAIL" in err
    assert "1 failed" in err


def test_validate_multiple_paths(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    good = find_rule_for_version("wombat", "ion", "0.9.11").path
    bad = tmp_path / "bad.json"
    bad.write_text("[[")
    rc = validate(good, bad)
    err = capsys.readouterr().err
    assert rc == 1  # any failure → 1
    assert "PASS" in err
    assert "FAIL" in err
    assert "2 document(s) checked, 1 failed" in err


def test_list_shows_nineteen_document_levels(capsys: pytest.CaptureFixture[str]) -> None:
    rc = list_rules()
    err = capsys.readouterr().err
    assert rc == 0
    lines = [line for line in err.splitlines() if line.strip()]
    assert len(lines) == 19
    assert "diann" in err
    assert "wombat" in err
    assert "ion" in err
    assert "peptidoform" in err


def test_export_schema_writes_file() -> None:
    rc = export_schema_cmd()
    assert rc == 0
    schema_path = packaged_rules_root() / "_schema" / "parse_rule.schema.json"
    assert schema_path.exists()
    assert schema_path.stat().st_size > 100
    document_schema = packaged_rules_root() / "_schema" / "parse_rule_document.schema.json"
    assert document_schema.exists()


def test_convert_with_explicit_rule_config_writes_h5ad(tmp_path: Path) -> None:
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

    rule_path = tmp_path / "rule.json"
    rule_path.write_text(json.dumps(_tiny_rule_document()))

    output_base = tmp_path / "out"
    output = output_base.with_suffix(".h5ad")
    stale = output_base.with_suffix(".h5mu")
    stale.write_text("stale")
    params_path = tmp_path / "params.txt"
    params_path.write_text("unregistered tool parameters")
    with pytest.raises(KeyError, match="no parameter parser"):
        convert(
            data_path,
            "ion",
            ConvertCliOptions(
                params=params_path,
                rule_config=rule_path,
                output=output_base,
            ),
        )
    assert not output.exists()
    assert stale.exists()

    without_params = tmp_path / "without_params"
    assert (
        convert(
            data_path,
            "ion",
            ConvertCliOptions(rule_config=rule_path, output=without_params),
        )
        == 0
    )
    import anndata as ad

    no_params_metadata = ad.read_h5ad(without_params.with_suffix(".h5ad")).uns["anndata_proteomics"]
    assert no_params_metadata["rule_selection_method"] == "rule_config"
    assert "search_parameters_version_status" not in no_params_metadata


def test_convert_with_multilevel_rule_config_writes_h5mu(tmp_path: Path) -> None:
    """Without LEVEL, an external document converts every matching level."""
    import mudata
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
    document = _tiny_rule_document(include_protein=True)
    rule_path = tmp_path / "rules.json"
    rule_path.write_text(json.dumps(document))

    output_base = tmp_path / "out"
    rc = convert(
        data_path,
        options=ConvertCliOptions(rule_config=rule_path, output=output_base),
    )

    assert rc == 0
    output = output_base.with_suffix(".h5mu")
    assert output.exists()
    result = mudata.read_h5mu(output)
    assert list(result.mod) == ["ion", "protein"]
    assert result.uns["anndata_proteomics"]["rule_selection_method"] == "rule_config"
    assert all(
        modality.uns["anndata_proteomics"]["rule_selection_method"] == "rule_config"
        for modality in result.mod.values()
    )


def test_convert_with_single_level_rule_config_writes_one_modality_h5mu(
    tmp_path: Path,
) -> None:
    import mudata
    import pandas as pd

    data_path = tmp_path / "tiny.tsv"
    pd.DataFrame(
        {
            "Run": ["S1"],
            "Sequence": ["P1"],
            "Charge": [2],
            "Intensity": [10.0],
        }
    ).to_csv(data_path, sep="\t", index=False)
    rule_path = tmp_path / "rules.json"
    rule_path.write_text(json.dumps(_tiny_rule_document()))

    rc = convert(
        data_path,
        options=ConvertCliOptions(rule_config=rule_path, output=tmp_path / "out"),
    )

    assert rc == 0
    assert list(mudata.read_h5mu(tmp_path / "out.h5mu").mod) == ["ion"]


def test_convert_rejects_output_extension(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    data_path = tmp_path / "tiny.tsv"
    output = tmp_path / "out.h5ad"

    rc = convert(
        data_path,
        options=ConvertCliOptions(rule_config=tmp_path / "rule.json", output=output),
    )

    assert rc == 2
    assert "extensionless basename" in capsys.readouterr().err


def test_convert_returns_one_when_vendor_not_detected(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    data_path = tmp_path / "unknown.csv"
    data_path.write_text("foo,bar,baz\n1,2,3\n")
    rc = convert(data_path)
    captured = capsys.readouterr()
    assert rc == 1
    assert "auto-detect" in (captured.out + captured.err).lower()


def test_convert_requires_params_without_rule_config(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # --software bypasses vendor detection, so the missing-params branch is reached.
    data_path = tmp_path / "x.tsv"
    data_path.write_text("Run\tSequence\nS1\tP1\n")
    rc = convert(data_path, options=ConvertCliOptions(software="diann"))
    captured = capsys.readouterr()
    assert rc == 1
    assert "--params" in (captured.out + captured.err)
