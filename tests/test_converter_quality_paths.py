"""Focused line and branch coverage for conversion helpers."""

from __future__ import annotations

import logging
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from anndata_proteomics.converters import assemble, pipeline, wide
from anndata_proteomics.converters._pieces import ConversionPieces
from anndata_proteomics.converters.long import _aggfunc_for, convert_long
from anndata_proteomics.converters.numeric import warn_if_all_missing
from anndata_proteomics.params.model import Parameters
from anndata_proteomics.rules.schema import (
    ColumnCompute,
    ParseRule,
    QuantificationLevel,
    SampleNameCleanup,
)


def _long_document() -> dict[str, Any]:
    return {
        "schema_version": "0.1",
        "file_version": "test",
        "software_name": "Synthetic",
        "software_version": "1",
        "input_shape": "long",
        "quantification_level": "ion",
        "axis": {
            "obs_keys": ["Run"],
            "var_keys": ["Feature"],
            "x_layer": "Intensity",
            "duplicates": {"mode": "keep_first"},
        },
        "columns": {
            "obs": {"select": {"Run": "Run"}},
            "var": {"select": {"Feature": "Feature"}},
        },
        "layers": [{"name": "Intensity", "source": "Intensity"}],
    }


def _long_rule() -> ParseRule:
    return ParseRule.model_validate(_long_document())


def _wide_document() -> dict[str, Any]:
    document = _long_document()
    document["input_shape"] = "wide"
    document["axis"] = {
        "obs_keys": ["sample"],
        "var_keys": ["Feature"],
        "x_layer": "Intensity",
        "duplicates": {"mode": "keep_first"},
    }
    document["columns"]["obs"] = {"select": {"sample": "<sample>"}}
    document["layers"] = [
        {
            "name": "Intensity",
            "source": r"^(?P<sample>S\d+) Intensity$",
        }
    ]
    return document


def _wide_rule() -> ParseRule:
    return ParseRule.model_validate(_wide_document())


def test_assemble_copies_extra_uns_and_helper_guards(
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    rule = _long_rule()
    pieces = ConversionPieces(
        X=np.asarray([[1.0]]),
        obs=pd.DataFrame({"Run": ["S1"]}, index=["S1"]),
        var=pd.DataFrame({"Feature": ["F1"]}, index=["F1"]),
        layers={"Intensity": np.asarray([[1.0]])},
        uns={"extra": {"value": 1}},
    )
    assert assemble.to_anndata(pieces, rule).uns["extra"]["value"] == 1

    group = rule.columns.var.model_copy(update={"select": {"Missing": "not-present"}})
    with pytest.raises(ValueError, match="cannot select"):
        assemble._materialize_column_group(pd.DataFrame({"Feature": ["F1"]}), group)

    coalesce = ColumnCompute.model_validate(
        {"name": "value", "from": ["one", "two"], "how": "coalesce"}
    )
    with pytest.raises(ValueError, match="source column"):
        assemble._compute_column(pd.DataFrame({"one": [1]}), coalesce)

    proforma = ColumnCompute.model_validate(
        {
            "name": "value",
            "from": ["source"],
            "how": "proforma_sequence",
        }
    )
    with pytest.raises(ValueError, match="APB column"):
        assemble._compute_column(pd.DataFrame(), proforma)

    ion = ColumnCompute.model_validate(
        {
            "name": "value",
            "from": ["sequence", "charge"],
            "how": "proforma_ion",
        }
    )
    with pytest.raises(ValueError, match="source column"):
        assemble._compute_column(pd.DataFrame({"sequence": ["PEP"]}), ion)

    fragment = ColumnCompute.model_validate(
        {
            "name": "value",
            "from": ["ion", "fragment"],
            "how": "proforma_fragment",
        }
    )
    with pytest.raises(ValueError, match="source column"):
        assemble._compute_column(pd.DataFrame({"ion": ["PEP/2"]}), fragment)

    unsupported = coalesce.model_copy(update={"how": "unsupported"})
    with pytest.raises(ValueError, match="unsupported column"):
        assemble._compute_column(pd.DataFrame({"one": [1], "two": [2]}), unsupported)

    for value, message in [
        (None, "missing"),
        (np.nan, "missing"),
        ("", "empty"),
        ("not-numeric", "numeric"),
    ]:
        with pytest.raises(ValueError, match=message):
            assemble._format_charge(value)

    target = ad.AnnData(np.ones((1, 1)))
    target.uns["anndata_proteomics"] = {}
    monkeypatch.setattr(assemble, "available_software", lambda: ("Synthetic",))
    monkeypatch.setattr(
        assemble,
        "parse_params",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(ValueError("bad params")),
    )
    assemble._attach_search_parameters(target, "params.txt", "Synthetic")
    assert target.uns["anndata_proteomics"]["search_parameters_path"] == "params.txt"
    assert target.uns["anndata_proteomics"]["search_parameters_error"] == ("ValueError: bad params")

    missing_version = ad.AnnData(np.ones((1, 1)))
    missing_version.uns["anndata_proteomics"] = {}
    monkeypatch.setattr(
        assemble,
        "parse_params",
        lambda *_args, **_kwargs: Parameters(software_name="Synthetic"),
    )
    assemble._attach_search_parameters(missing_version, "params.txt", "Synthetic")
    assert (
        missing_version.uns["anndata_proteomics"]["search_parameters_version_status"] == "missing"
    )
    assert "no software version" in caplog.text


def test_assemble_column_selection_and_fragment_column_inventory() -> None:
    plain_rule = _long_rule()
    frame = pd.DataFrame(
        {
            "Run": ["S1"],
            "Feature": ["F1"],
            "Intensity": [1.0],
            "unused": ["x"],
        }
    )
    assert assemble._columns_needed_for_long(frame, plain_rule) == [
        "Run",
        "Feature",
        "Intensity",
    ]

    document = _long_document()
    document["quantification_level"] = "fragment"
    document["fragments"] = {
        "label_strategy": "column",
        "label_column": "Fragment.Info",
        "label_output": "Fragment",
        "value_columns": ["Fragment.Intensity"],
        "delimiter": ";",
    }
    fragment_rule = ParseRule.model_validate(document)
    fragment_frame = frame.assign(
        **{
            "Fragment.Info": ["b1/1;"],
            "Fragment.Intensity": ["2;"],
        }
    )
    assert "Fragment.Info" in assemble._columns_needed_for_long(
        fragment_frame,
        fragment_rule,
    )


def test_pipeline_empty_targets_version_failures_and_mudata_paths(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    assert pipeline.param_version(None, "tool") is None
    monkeypatch.setattr(
        pipeline,
        "parse_params",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(ValueError("bad")),
    )
    assert pipeline.param_version(tmp_path / "bad.params", "tool") is None

    monkeypatch.setattr(pipeline, "convertible_levels", lambda *_args, **_kwargs: [])
    assert pipeline.available_targets("tool", None, []) == []
    with pytest.raises(ValueError, match="no level resolves"):
        pipeline.build_mudata(pd.DataFrame(), "tool", None)
    with pytest.raises(ValueError, match="no levels supplied"):
        pipeline.build_mudata_from_rules(pd.DataFrame(), {})

    rule = _long_rule()
    monkeypatch.setattr(
        pipeline,
        "convertible_levels",
        lambda *_args, **_kwargs: ["ion"],
    )
    monkeypatch.setattr(
        pipeline,
        "_select_rule",
        lambda *_args, **_kwargs: (rule, "software_version"),
    )
    captured: dict[str, Any] = {}

    def fake_builder(
        _df: pd.DataFrame,
        rules: dict[str, ParseRule],
        **kwargs: Any,
    ) -> str:
        captured["rules"] = rules
        captured["kwargs"] = kwargs
        return "built"

    monkeypatch.setattr(pipeline, "build_mudata_from_rules", fake_builder)
    logs: list[str] = []
    assert (
        pipeline.build_mudata(
            pd.DataFrame({"x": [1]}),
            "tool",
            "1",
            log=logs.append,
        )
        == "built"
    )
    assert set(captured["rules"]) == {"ion"}
    assert logs and "skipping levels" in logs[0]

    monkeypatch.setattr(
        pipeline,
        "convertible_levels",
        lambda *_args, **_kwargs: list(pipeline.LEVELS),
    )
    logs.clear()
    assert (
        pipeline.build_mudata(
            pd.DataFrame({"x": [1]}),
            "tool",
            "1",
            log=logs.append,
        )
        == "built"
    )
    assert not logs


def test_parameter_resolution_distinguishes_missing_from_parse_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "params.txt"
    source.write_text("params", encoding="utf-8")
    monkeypatch.setattr(
        pipeline,
        "parse_params",
        lambda *_args, **_kwargs: Parameters(software_name="PEAKS"),
    )

    missing = pipeline.resolve_parameters(source, "peaks")

    assert missing.version is None
    assert missing.version_status == "missing"
    assert missing.error is None
    monkeypatch.setattr(
        pipeline,
        "parse_params",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(ValueError("wrong tool")),
    )

    failed = pipeline.resolve_parameters(source, "peaks")

    assert failed.parameters is None
    assert failed.version_status == "parse_error"
    assert failed.error == "ValueError: wrong tool"


def test_build_mudata_resolves_params_before_materializing_rules(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    parameters = Parameters(
        software_version="2.6.0",
        acquisition_method="DDA",
    )
    resolution = pipeline.ParameterResolution(
        source_path=tmp_path / "params.txt",
        parameters=parameters,
        version="2.6.0",
        version_status="present",
    )
    captured: dict[str, object] = {}
    monkeypatch.setattr(pipeline, "resolve_parameters", lambda *_args: resolution)

    def convertible_levels(
        _slug: str,
        version: str | None,
        _headers: object,
        **kwargs: object,
    ) -> list[QuantificationLevel]:
        captured["version"] = version
        captured["search_parameters"] = kwargs["search_parameters"]
        return ["ion"]

    monkeypatch.setattr(pipeline, "convertible_levels", convertible_levels)
    monkeypatch.setattr(
        pipeline,
        "_select_rule",
        lambda *_args, **_kwargs: (_long_rule(), "software_version"),
    )
    monkeypatch.setattr(
        pipeline,
        "build_mudata_from_rules",
        lambda *_args, **_kwargs: "built",
    )

    result = pipeline.build_mudata(
        pd.DataFrame({"Intensity": [1.0]}),
        "diann",
        None,
        params_path=resolution.source_path,
    )

    assert result == "built"
    assert captured == {
        "version": "2.6.0",
        "search_parameters": parameters,
    }


def test_parse_error_resolution_attaches_error_without_selection_method(
    tmp_path: Path,
) -> None:
    target = ad.AnnData(np.ones((1, 1)))
    resolution = pipeline.ParameterResolution(
        source_path=tmp_path / "params.txt",
        parameters=None,
        version=None,
        version_status="parse_error",
        error="ValueError: wrong tool",
    )

    pipeline.attach_parameter_resolution(
        target,
        resolution,
        selection_method=None,
        warn_missing=False,
    )

    metadata = target.uns["anndata_proteomics"]
    assert "rule_selection_method" not in metadata
    assert metadata["search_parameters_version_status"] == "parse_error"
    assert metadata["search_parameters_error"] == "ValueError: wrong tool"

    no_error_detail = ad.AnnData(np.ones((1, 1)))
    pipeline.attach_parameter_resolution(
        no_error_detail,
        pipeline.ParameterResolution(
            source_path=tmp_path / "unknown.params",
            parameters=None,
            version=None,
            version_status="parse_error",
        ),
        selection_method="version_unavailable",
        warn_missing=False,
    )
    assert "search_parameters_error" not in no_error_detail.uns["anndata_proteomics"]


def test_build_mudata_attaches_parameters_and_prefixes_features(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from anndata_proteomics.converters import assemble as assemble_module

    rule = _long_rule()
    converted = ad.AnnData(
        np.ones((1, 1)),
        obs=pd.DataFrame(index=["S1"]),
        var=pd.DataFrame(index=["F1"]),
    )
    monkeypatch.setattr(
        assemble_module,
        "convert",
        lambda *_args, **_kwargs: converted.copy(),
    )
    parameters = SimpleNamespace(software_name="Synthetic")
    monkeypatch.setattr(pipeline, "parse_params", lambda *_args, **_kwargs: parameters)
    writes: list[tuple[Any, Any, str]] = []
    monkeypatch.setattr(
        pipeline,
        "write_search_parameters",
        lambda target, params, *, source_path: writes.append((target, params, source_path)),
    )
    monkeypatch.setattr(pipeline, "store_quantification_summary", lambda _target: None)

    result = pipeline.build_mudata_from_rules(
        pd.DataFrame({"x": [1]}),
        {"ion": rule},
        params_path="params.txt",
        software="synthetic",
    )

    assert result.mod["ion"].var_names.tolist() == ["ion:F1"]
    assert writes[0][1:] == (parameters, "params.txt")


def test_build_mudata_records_parameter_and_selection_provenance(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    from anndata_proteomics.converters import assemble as assemble_module

    converted = ad.AnnData(
        np.ones((1, 1)),
        obs=pd.DataFrame(index=["S1"]),
        var=pd.DataFrame(index=["F1"]),
    )
    monkeypatch.setattr(
        assemble_module,
        "convert",
        lambda *_args, **_kwargs: converted.copy(),
    )
    resolution = pipeline.ParameterResolution(
        source_path=tmp_path / "params.txt",
        parameters=Parameters(software_name="PEAKS"),
        version=None,
        version_status="missing",
    )

    result = pipeline.build_mudata_from_rules(
        pd.DataFrame({"x": [1]}),
        {"ion": _long_rule()},
        parameter_resolution=resolution,
        rule_selection_method="columns",
    )

    for target in (result, result.mod["ion"]):
        metadata = target.uns["anndata_proteomics"]
        assert metadata["rule_selection_method"] == "columns"
        assert metadata["search_parameters_version_status"] == "missing"
        assert metadata["search_parameters_path"] == str(resolution.source_path)
    assert caplog.text.count("selected rule by columns") == 1


def test_matching_rules_excludes_nonmatches(monkeypatch: pytest.MonkeyPatch) -> None:
    rules: dict[QuantificationLevel, ParseRule] = {"ion": _long_rule()}
    monkeypatch.setattr(pipeline, "matches", lambda *_args: False)
    assert pipeline.matching_rules(rules, ["unrelated"]) == {}


def test_long_converter_wrong_shape_factor_and_raw_duplicate_mode() -> None:
    wide_rule = _wide_rule()
    with pytest.raises(ValueError, match="convert_long called"):
        convert_long(pd.DataFrame(), wide_rule)

    keep_all = _long_rule().model_copy(
        update={
            "axis": _long_rule().axis.model_copy(
                update={
                    "duplicates": _long_rule().axis.duplicates.model_copy(
                        update={"mode": "keep_all_as_raw_table"}
                    )
                }
            )
        }
    )
    with pytest.raises(NotImplementedError, match="keep_all"):
        _aggfunc_for(keep_all)

    document = _long_document()
    document["layers"][0].update(
        {
            "encoding_mode": "factor",
            "categories": {"missing": 0, "found": 1},
        }
    )
    frame = pd.DataFrame(
        {
            "Run": ["S1"],
            "Feature": ["F1"],
            "Intensity": ["found"],
        }
    )
    assert convert_long(frame, ParseRule.model_validate(document)).X[0, 0] == 1


def test_wide_helpers_and_converter_guards() -> None:
    assert wide._matching_columns(["S1 Intensity"], r"^S1 Intensity$") == [
        ("S1 Intensity", "S1 Intensity")
    ]

    rule = _wide_rule()
    layer = rule.layers[0].model_copy(
        update={
            "source": r"^(?P<sample>S1) Intensity (?:A|B)$",
        }
    )
    duplicate_columns = pd.DataFrame(
        {
            "Feature": ["F1"],
            "S1 Intensity A": [1.0],
            "S1 Intensity B": [2.0],
        }
    )
    with pytest.raises(ValueError, match="multiple columns"):
        wide._gather_layer_matrix(
            duplicate_columns,
            layer,
            ["S1"],
            pd.Index(["F1"]),
            ["Feature"],
            "error",
        )
    empty_for_sample = wide._gather_layer_matrix(
        pd.DataFrame({"Feature": ["F1"]}),
        rule.layers[0],
        ["S1"],
        pd.Index(["F1"]),
        ["Feature"],
        "keep_first",
    )
    assert np.isnan(empty_for_sample).all()

    grouped_cleanup = rule.model_copy(
        update={"sample_name_cleanup": SampleNameCleanup(pattern=r"run_(.)")}
    )
    assert wide._apply_sample_cleanup(["run_A", "other"], grouped_cleanup) == [
        "A",
        "other",
    ]
    whole_cleanup = rule.model_copy(
        update={"sample_name_cleanup": SampleNameCleanup(pattern=r"run_[A-Z]")}
    )
    assert wide._apply_sample_cleanup(["run_A"], whole_cleanup) == ["run_A"]

    with pytest.raises(ValueError, match="convert_wide called"):
        wide.convert_wide(pd.DataFrame(), _long_rule())

    no_match = pd.DataFrame({"Feature": ["F1"], "unrelated": [1.0]})
    with pytest.raises(ValueError, match="no columns matched"):
        wide.convert_wide(no_match, rule)

    invalid_obs_document = _wide_document()
    invalid_obs_document["columns"]["obs"] = {"select": {"sample": "unsupported"}}
    with pytest.raises(ValueError, match="only the"):
        wide.convert_wide(
            pd.DataFrame({"Feature": ["F1"], "S1 Intensity": [1.0]}),
            ParseRule.model_validate(invalid_obs_document),
        )

    optional_document = _wide_document()
    optional_document["layers"].append(
        {
            "name": "Optional",
            "source": r"^(?P<sample>S\d+) Missing$",
        }
    )
    pieces = wide.convert_wide(
        pd.DataFrame({"Feature": ["F1"], "S1 Intensity": [1.0]}),
        ParseRule.model_validate(optional_document),
    )
    assert "Optional" not in pieces.layers


def test_warn_if_all_missing_ignores_an_empty_matrix(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """A zero-feature layer has no missingness to report and must stay quiet."""
    with caplog.at_level(logging.WARNING):
        warn_if_all_missing(np.empty((0, 0), dtype="float64"), "Empty")

    assert caplog.text == ""
