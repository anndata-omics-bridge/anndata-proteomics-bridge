"""Focused edge-path tests for small APB support modules."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pandas as pd
import pytest

from anndata_proteomics import _logging, test_data
from anndata_proteomics.adapters.anndata.result import load_converted_result
from anndata_proteomics.annotation.loader import load_annotation
from anndata_proteomics.converters._axis import join_keys
from anndata_proteomics.converters._fragments import (
    _fragment_positions,
    _split_packed,
    explode_fragments,
)
from anndata_proteomics.rules import _discovery
from anndata_proteomics.rules.schema import ColumnLabeledFragments


def test_default_logging_sink_is_replaced(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[tuple[str, tuple[Any, ...], dict[str, Any]]] = []
    fake_logger = SimpleNamespace(
        remove=lambda: calls.append(("remove", (), {})),
        add=lambda *args, **kwargs: calls.append(("add", args, kwargs)),
    )
    monkeypatch.setattr(_logging, "logger", fake_logger)

    _logging.configure_default_sink("DEBUG")

    assert [call[0] for call in calls] == ["remove", "add"]
    assert calls[1][2] == {"format": _logging.DEFAULT_FORMAT, "level": "DEBUG"}


def test_annotation_loader_rejects_non_files_and_empty_tables(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="not a file"):
        load_annotation(tmp_path)

    empty_toml = tmp_path / "empty.toml"
    empty_toml.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="no .*samples"):
        load_annotation(empty_toml)

    invalid_samples = tmp_path / "invalid.toml"
    invalid_samples.write_text('samples = ["not-a-table"]', encoding="utf-8")
    with pytest.raises(ValueError, match="must be tables"):
        load_annotation(invalid_samples)

    empty_csv = tmp_path / "empty.csv"
    empty_csv.write_text("raw_file\n", encoding="utf-8")
    with pytest.raises(ValueError, match="at least one"):
        load_annotation(empty_csv)


def test_axis_and_fragment_helper_guards() -> None:
    assert join_keys(pd.Series(["a", 2])) == "a_2"
    assert _split_packed(None, ";") == []
    assert _split_packed(pd.NA, ";") == []
    assert _split_packed(pd.NaT, ";") == []
    assert _split_packed(float("nan"), ";") == []
    assert _split_packed("  ", ";") == []
    assert _split_packed(";;", ";") == []
    assert _fragment_positions(["b1", "y2"]) == [0, 1]

    fragments = ColumnLabeledFragments(
        delimiter=";",
        label_strategy="column",
        label_column="Fragment.Info",
        label_output="Fragment",
        value_columns=["Intensity"],
    )
    with pytest.raises(KeyError, match="missing from the input"):
        explode_fragments(pd.DataFrame({"other": [1]}), fragments)


def test_result_loader_rejects_unknown_suffix(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="unsupported converted result"):
        load_converted_result(tmp_path / "result.txt")


def test_rule_discovery_rejects_external_path_and_missing_vendor(
    tmp_path: Path,
) -> None:
    with pytest.raises(ValueError, match="not a packaged rule"):
        _discovery.document_vendor(tmp_path / "rules.json")
    assert _discovery.document_paths_for_software("definitely-missing") == ()


def test_data_lookup_empty_and_present_paths(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(test_data, "DOWNLOADED_DB", tmp_path / "missing.csv")
    monkeypatch.setattr(test_data, "TEST_DATA_DIR", tmp_path)
    assert test_data.find_test_data("missing") == test_data.VendorDataUnavailable("missing")
    assert test_data._module_for_dataset(
        tmp_path / "dataset", test_data_dir=tmp_path
    ) == test_data.DatasetModuleUnavailable(tmp_path / "dataset")
    assert test_data.find_fasta_for_module(
        "unknown", test_data_dir=tmp_path
    ) == test_data.FastaUnavailable("unknown")
    assert test_data.find_fasta_for_module(
        "dda_qexactive", test_data_dir=tmp_path
    ) == test_data.FastaUnavailable("dda_qexactive")

    index = tmp_path / "raw_file_db_downloaded.csv"
    index.write_text(
        "software_name,status,input_file_path,module\n"
        "Tool,failed,collection/failure,dia_aif\n"
        "Other,ok,collection/success,dia_astral\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(test_data, "DOWNLOADED_DB", index)
    assert test_data.find_test_data("Tool") == test_data.VendorDataUnavailable("Tool")
    assert test_data._module_for_dataset(
        tmp_path / "unmatched", test_data_dir=tmp_path
    ) == test_data.DatasetModuleUnavailable(tmp_path / "unmatched")

    monkeypatch.setattr(test_data, "PARAM_FIXTURE_DIR", tmp_path)
    assert test_data.find_param_file("unknown") == test_data.ParameterFileUnavailable("unknown")
    assert test_data.find_param_file("DIA-NN") == test_data.ParameterFileUnavailable("DIA-NN")
