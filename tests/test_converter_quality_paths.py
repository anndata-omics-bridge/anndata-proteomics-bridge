"""Focused line and branch coverage for conversion helpers."""

from __future__ import annotations

import inspect
import logging
from collections.abc import Callable
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from anndata_proteomics.adapters.anndata import conversion as conversion_adapter
from anndata_proteomics.adapters.anndata.params import require_search_parameters
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
from anndata_proteomics.workflows import conversion as conversion_workflow
from anndata_proteomics.workflows.conversion import LevelConversion


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


def test_in_memory_boundaries_drive_conversion_workflow() -> None:
    source = {
        "table": pd.DataFrame(
            {
                "Run": ["S1"],
                "Feature": ["F1"],
                "Intensity": [7.0],
            }
        )
    }
    persisted: dict[QuantificationLevel, LevelConversion] = {}

    def read_table() -> pd.DataFrame:
        return source["table"].copy()

    def write_result(result: LevelConversion) -> None:
        persisted[result.level] = result

    rule = _long_rule()
    calculated = conversion_workflow.convert_selected_level(
        read_table(),
        rule.quantification_level,
        pipeline.RuleSelection(rule, "rule_config"),
    )
    write_result(calculated)

    assert persisted == {"ion": calculated}
    assert calculated.pieces.X.shape == (1, 1)
    assert calculated.pieces.X[0, 0] == 7.0
    assert source["table"].columns.tolist() == ["Run", "Feature", "Intensity"]


def test_assemble_helper_guards(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rule = _long_rule()
    pieces = ConversionPieces(
        X=np.asarray([[1.0]]),
        obs=pd.DataFrame({"Run": ["S1"]}, index=["S1"]),
        var=pd.DataFrame({"Feature": ["F1"]}, index=["F1"]),
        layers={"Intensity": np.asarray([[1.0]])},
    )
    conversion = LevelConversion(
        level="ion",
        selection=pipeline.RuleSelection(rule, "software_version"),
        pieces=pieces,
    )
    assert set(conversion_adapter.to_anndata(conversion).uns) == {"anndata_proteomics"}

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


def test_empty_level_conversion_and_persistence_are_rejected() -> None:
    with pytest.raises(ValueError, match="no levels supplied"):
        conversion_workflow.convert_selected_levels(pd.DataFrame(), {})
    with pytest.raises(ValueError, match="no level conversions supplied"):
        conversion_adapter.to_mudata({})


@pytest.mark.parametrize(
    "converter",
    [
        conversion_workflow.convert_selected_level,
        conversion_workflow.convert_selected_levels,
        conversion_workflow.convert_level_from_parameters,
        conversion_workflow.select_rule_from_parameters,
        conversion_workflow.select_rules_from_parameters,
    ],
)
def test_pipeline_diagnostics_do_not_use_callback_contract(
    converter: Callable[..., object],
) -> None:
    assert "log" not in inspect.signature(converter).parameters


def test_available_rules_skip_unavailable_levels_without_exception_control_flow(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rule = _long_rule()
    resolved: list[QuantificationLevel] = []

    def find(
        _slug: str,
        level: QuantificationLevel,
        _version: str,
        _headers: object,
    ) -> pipeline.RuleLookup:
        resolved.append(level)
        if level == "ion":
            return pipeline.RuleSelection(rule, "software_version")
        return pipeline.RuleUnavailable(f"{level} is unavailable")

    monkeypatch.setattr(pipeline, "find_rule_for_version", find)

    selections = pipeline.available_rules_for_version("tool", "1", [])
    assert list(selections) == ["ion"]
    assert resolved == list(pipeline.LEVELS)


def test_parameter_resolution_reports_missing_version_and_propagates_parse_failure(
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

    assert missing.version == pipeline.MissingRuleVersion()
    monkeypatch.setattr(
        pipeline,
        "parse_params",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(ValueError("wrong tool")),
    )

    with pytest.raises(ValueError, match="wrong tool"):
        pipeline.resolve_parameters(source, "peaks")


def test_parameter_resolution_persistence_is_separate_from_rule_selection(
    tmp_path: Path,
) -> None:
    target = ad.AnnData(np.ones((1, 1)))
    resolution = pipeline.ParameterResolution(
        source_path=tmp_path / "params.txt",
        parameters=Parameters(software_name="Synthetic"),
        version=pipeline.MissingRuleVersion(),
    )

    conversion_adapter.write_parameter_resolution(target, resolution)

    metadata = target.uns["anndata_proteomics"]
    assert "rule_selection_method" not in metadata
    assert metadata["search_parameters_version_status"] == "missing"
    stored = require_search_parameters(target)
    assert stored.software_name == "Synthetic"


def test_mudata_adapter_prefixes_features_and_preserves_selection_provenance() -> None:
    rule = _long_rule()
    conversion = LevelConversion(
        level="ion",
        selection=pipeline.RuleSelection(rule, "software_version"),
        pieces=ConversionPieces(
            X=np.ones((1, 1)),
            obs=pd.DataFrame(index=["S1"]),
            var=pd.DataFrame(index=["F1"]),
            layers={"Intensity": np.ones((1, 1))},
        ),
    )

    result = conversion_adapter.to_mudata({"ion": conversion})

    assert result.mod["ion"].var_names.tolist() == ["ion:F1"]
    for target in (result, result.mod["ion"]):
        assert target.uns["anndata_proteomics"]["rule_selection_method"] == "software_version"


def test_parameter_resolution_can_be_persisted_to_collection_and_modalities(
    tmp_path: Path,
) -> None:
    rule = _long_rule()
    conversion = LevelConversion(
        level="ion",
        selection=pipeline.RuleSelection(rule, "columns"),
        pieces=ConversionPieces(
            X=np.ones((1, 1)),
            obs=pd.DataFrame(index=["S1"]),
            var=pd.DataFrame(index=["F1"]),
            layers={"Intensity": np.ones((1, 1))},
        ),
    )
    result = conversion_adapter.to_mudata({"ion": conversion})
    resolution = pipeline.ParameterResolution(
        source_path=tmp_path / "params.txt",
        parameters=Parameters(software_name="PEAKS"),
        version=pipeline.MissingRuleVersion(),
    )
    for target in (result, result.mod["ion"]):
        conversion_adapter.write_parameter_resolution(target, resolution)

    for target in (result, result.mod["ion"]):
        metadata = target.uns["anndata_proteomics"]
        assert metadata["rule_selection_method"] == "columns"
        assert metadata["search_parameters_version_status"] == "missing"
        assert metadata["search_parameters_path"] == str(resolution.source_path)
        stored = require_search_parameters(target)
        assert stored.software_name == "PEAKS"


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
    assert wide._matching_columns(["S1 Intensity"], r"^(?P<sample>S1) Intensity$") == [
        ("S1 Intensity", "S1")
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
            list(duplicate_columns.columns),
            ["S1"],
            pd.Index(["F1"]),
            ["Feature"],
            "error",
        )
    empty_for_sample = wide._gather_layer_matrix(
        pd.DataFrame({"Feature": ["F1"]}),
        rule.layers[0],
        ["Feature"],
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
