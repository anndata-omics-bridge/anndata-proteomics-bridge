"""Focused branch coverage for APB command functions."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from anndata_proteomics.adapters.anndata import conversion as conversion_adapter
from anndata_proteomics.converters import pipeline as conversion_pipeline
from anndata_proteomics.params.model import Parameters
from anndata_proteomics.rules.loader import load_packaged_rule_for_version
from anndata_proteomics.scripts import cli
from anndata_proteomics.workflows import conversion as conversion_workflow


class _Container:
    def __init__(self, *, modalities: dict[str, Any] | None = None) -> None:
        self.uns: dict[str, Any] = {}
        self.n_obs = 1
        if modalities is not None:
            self.mod = modalities
        self.written: Path | None = None

    def write_h5mu(self, path: Path) -> None:
        self.written = path
        path.write_text("mudata", encoding="utf-8")

    def write_h5ad(self, path: Path) -> None:
        self.written = path
        path.write_text("anndata", encoding="utf-8")


def test_rule_config_missing_level_and_no_match(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    frame = pd.DataFrame({"x": [1]})
    document = SimpleNamespace(
        levels={"ion": object()},
        effective_rules=lambda: {},
        software_name="Tool",
    )
    monkeypatch.setattr(cli, "read_table_columns", lambda _path: list(frame.columns))
    monkeypatch.setattr(cli, "_read_table_for_selections", lambda *_args: frame)
    monkeypatch.setattr(cli, "load_rule_document", lambda _path: document)
    assert (
        cli.convert(
            tmp_path / "data.tsv",
            "protein",
            cli.ConvertCliOptions(rule_config=tmp_path / "rule.json"),
        )
        == 1
    )
    assert (
        cli.convert(
            tmp_path / "data.tsv",
            options=cli.ConvertCliOptions(rule_config=tmp_path / "rule.json"),
        )
        == 1
    )


def test_rule_config_materializes_single_and_mudata_with_parameters(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    frame = pd.DataFrame({"x": [1]})
    parameters = Parameters(
        software_name="DIA-NN",
        software_version="2.6.0",
        acquisition_method="DDA",
    )
    resolution = conversion_pipeline.ParameterResolution(
        source_path=tmp_path / "params.txt",
        parameters=parameters,
        version=conversion_pipeline.PresentRuleVersion("2.6.0"),
    )
    parameterized_rule_calls: list[tuple[str, Parameters]] = []
    parameterized_rules_calls: list[Parameters] = []
    rule = object()

    def parameterized_effective_rule(level: str, search_parameters: Parameters) -> object:
        parameterized_rule_calls.append((level, search_parameters))
        return rule

    def parameterized_effective_rules(search_parameters: Parameters) -> dict[str, object]:
        parameterized_rules_calls.append(search_parameters)
        return {"ion": rule}

    document = SimpleNamespace(
        levels={"ion": object()},
        parameterized_effective_rule=parameterized_effective_rule,
        parameterized_effective_rules=parameterized_effective_rules,
        software_name="DIA-NN",
    )
    monkeypatch.setattr(cli, "read_table_columns", lambda _path: list(frame.columns))
    monkeypatch.setattr(cli, "_read_table_for_selections", lambda *_args: frame)
    monkeypatch.setattr(cli, "load_rule_document", lambda _path: document)
    monkeypatch.setattr(
        conversion_pipeline,
        "resolve_parameters",
        lambda *_args: resolution,
    )
    converted = ad.AnnData(np.ones((1, 1), dtype=np.float32))
    level_conversion = object()
    monkeypatch.setattr(
        conversion_workflow,
        "convert_selected_level",
        lambda *_args, **_kwargs: level_conversion,
    )
    monkeypatch.setattr(
        conversion_adapter,
        "to_anndata",
        lambda conversion: converted if conversion is level_conversion else None,
    )
    monkeypatch.setattr(
        conversion_adapter,
        "write_parameter_resolution",
        lambda *_args: None,
    )

    assert (
        cli.convert(
            tmp_path / "data.tsv",
            "ion",
            cli.ConvertCliOptions(
                params=resolution.source_path,
                rule_config=tmp_path / "rule.json",
                output=tmp_path / "single",
            ),
        )
        == 0
    )
    assert parameterized_rule_calls == [("ion", parameters)]

    container = _Container(modalities={"ion": converted})
    monkeypatch.setattr(
        conversion_pipeline,
        "matching_rules",
        lambda rules, _headers: rules,
    )
    monkeypatch.setattr(
        conversion_workflow,
        "convert_selected_levels",
        lambda *_args, **_kwargs: {"ion": level_conversion},
    )
    monkeypatch.setattr(conversion_adapter, "to_mudata", lambda _conversions: container)
    assert (
        cli.convert(
            tmp_path / "data.tsv",
            options=cli.ConvertCliOptions(
                params=resolution.source_path,
                rule_config=tmp_path / "rule.json",
                output=tmp_path / "multi",
            ),
        )
        == 0
    )
    assert parameterized_rules_calls == [parameters]


def test_packaged_level_and_mudata_conversion_paths(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    frame = pd.DataFrame({"x": [1]})
    adata = ad.AnnData(np.ones((1, 1), dtype=np.float32))
    monkeypatch.setattr(cli, "read_table_columns", lambda _path: list(frame.columns))
    monkeypatch.setattr(cli, "_read_table_for_selections", lambda *_args: frame)
    monkeypatch.setattr(
        conversion_pipeline,
        "recognize_software",
        lambda _columns: conversion_pipeline.RecognizedSoftware("diann"),
    )
    resolution = conversion_pipeline.ParameterResolution(
        source_path=tmp_path / "params.txt",
        parameters=Parameters(software_name="DIA-NN", software_version="2.0"),
        version=conversion_pipeline.PresentRuleVersion("2.0"),
    )
    monkeypatch.setattr(
        conversion_pipeline,
        "resolve_parameters",
        lambda *_args: resolution,
    )
    level_conversion = object()
    selected_rule = load_packaged_rule_for_version("diann", "ion", "2.0.0")
    monkeypatch.setattr(
        conversion_workflow,
        "select_rule_from_parameters",
        lambda *_args: conversion_pipeline.RuleSelection(selected_rule, "software_version"),
    )
    monkeypatch.setattr(
        conversion_workflow,
        "convert_selected_level",
        lambda *_args, **_kwargs: level_conversion,
    )
    monkeypatch.setattr(conversion_adapter, "to_anndata", lambda _conversion: adata)
    monkeypatch.setattr(
        conversion_adapter,
        "write_parameter_resolution",
        lambda *_args: None,
    )
    params = resolution.source_path
    params.write_text("params", encoding="utf-8")
    output = tmp_path / "single"
    assert (
        cli.convert(
            tmp_path / "data.tsv",
            "ion",
            cli.ConvertCliOptions(params=params, output=output),
        )
        == 0
    )
    assert output.with_suffix(".h5ad").exists()

    container = _Container(modalities={"ion": adata})
    captured_resolutions: list[conversion_pipeline.ParameterResolution] = []

    def select_rules(
        _headers: object,
        _slug: str,
        selected_resolution: conversion_pipeline.ParameterResolution,
    ) -> dict[str, conversion_pipeline.RuleSelection]:
        captured_resolutions.append(selected_resolution)
        rule = load_packaged_rule_for_version("diann", "ion", "2.0.0")
        return {"ion": conversion_pipeline.RuleSelection(rule, "software_version")}

    monkeypatch.setattr(conversion_workflow, "select_rules_from_parameters", select_rules)
    monkeypatch.setattr(
        conversion_workflow,
        "convert_selected_levels",
        lambda *_args, **_kwargs: {"ion": level_conversion},
    )
    monkeypatch.setattr(conversion_adapter, "to_mudata", lambda _conversions: container)
    stale = (tmp_path / "multi").with_suffix(".h5ad")
    stale.write_text("stale", encoding="utf-8")
    assert (
        cli.convert(
            tmp_path / "data.tsv",
            options=cli.ConvertCliOptions(params=params, output=tmp_path / "multi"),
        )
        == 0
    )
    assert (tmp_path / "multi.h5mu").exists()
    assert not stale.exists()
    assert captured_resolutions == [resolution]

    monkeypatch.setattr(
        conversion_workflow,
        "select_rules_from_parameters",
        lambda *_args: {},
    )
    assert cli.convert(tmp_path / "data.tsv", options=cli.ConvertCliOptions(params=params)) == 1


def test_compound_conversion_separates_parameter_and_rule_software(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    frame = pd.DataFrame({"x": [1]})
    monkeypatch.setattr(cli, "read_table_columns", lambda _path: list(frame.columns))
    monkeypatch.setattr(cli, "_read_table_for_selections", lambda *_args: frame)
    monkeypatch.setattr(
        conversion_pipeline,
        "recognize_software",
        lambda _columns: conversion_pipeline.RecognizedSoftware("diann"),
    )
    parameter_path = tmp_path / "fragpipe.workflow"
    parameter_path.write_text("workflow", encoding="utf-8")
    resolution = conversion_pipeline.ParameterResolution(
        source_path=parameter_path,
        parameters=Parameters(
            software_name="FragPipe",
            software_version="24.0",
            quantification_software="DIA-NN",
            quantification_software_version="1.8.2 beta 8",
        ),
        version=conversion_pipeline.PresentRuleVersion("24.0"),
    )
    parser_calls: list[str] = []

    def resolve_parameters(
        _path: Path,
        software: str,
    ) -> conversion_pipeline.ParameterResolution:
        parser_calls.append(software)
        return resolution

    selected: dict[str, object] = {}
    level_conversion = object()

    def select_rule_from_parameters(
        _headers: object,
        slug: str,
        level: str,
        selected_resolution: conversion_pipeline.ParameterResolution,
    ) -> conversion_pipeline.RuleSelection:
        selected.update(
            slug=slug,
            level=level,
            parameter_resolution=selected_resolution,
        )
        rule = load_packaged_rule_for_version("diann", "ion", "1.8.2 beta 8")
        return conversion_pipeline.RuleSelection(rule, "software_version")

    monkeypatch.setattr(conversion_pipeline, "resolve_parameters", resolve_parameters)
    monkeypatch.setattr(
        conversion_workflow,
        "select_rule_from_parameters",
        select_rule_from_parameters,
    )
    monkeypatch.setattr(
        conversion_workflow,
        "convert_selected_level",
        lambda *_args, **_kwargs: level_conversion,
    )
    monkeypatch.setattr(
        conversion_adapter,
        "to_anndata",
        lambda _conversion: ad.AnnData(np.ones((1, 1), dtype=np.float32)),
    )
    monkeypatch.setattr(
        conversion_adapter,
        "write_parameter_resolution",
        lambda *_args: None,
    )

    assert (
        cli.convert(
            tmp_path / "report.tsv",
            "ion",
            cli.ConvertCliOptions(
                params=parameter_path,
                software="diann",
                params_software="fragpipe",
                output=tmp_path / "converted",
            ),
        )
        == 0
    )
    assert parser_calls == ["fragpipe"]
    assert selected == {
        "slug": "diann",
        "level": "ion",
        "parameter_resolution": resolution,
    }
    assert conversion_pipeline.resolve_rule_version(
        resolution,
        "diann",
    ) == conversion_pipeline.PresentRuleVersion("1.8.2 beta 8")


def test_summary_and_annotation_writers(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from anndata_proteomics.adapters.anndata import annotation as annotation_adapter
    from anndata_proteomics.adapters.anndata import description as summary
    from anndata_proteomics.adapters.anndata import result as result_adapter
    from anndata_proteomics.annotation import loader
    from anndata_proteomics.workflows import sample_annotation

    monkeypatch.setattr(summary, "describe_path", lambda *_args, **_kwargs: {"n_obs": 2})
    assert cli.summary_cmd(tmp_path / "data.h5ad", json=True) == 0
    assert json.loads(capsys.readouterr().out) == {"n_obs": 2}

    loaded = SimpleNamespace(table=object(), origin=object())
    annotation_result = object()
    monkeypatch.setattr(loader, "load_annotation", lambda _path: loaded)
    monkeypatch.setattr(
        annotation_adapter,
        "read_observation_frames",
        lambda _target: (pd.DataFrame(),),
    )
    monkeypatch.setattr(
        sample_annotation,
        "run_sample_annotation",
        lambda *_args: annotation_result,
    )
    monkeypatch.setattr(
        annotation_adapter,
        "write_sample_annotation",
        lambda _target, _result: None,
    )

    mudata = _Container(modalities={})
    monkeypatch.setattr(result_adapter, "load_converted_result", lambda _path: mudata)
    assert cli.annotate(tmp_path / "data.h5mu", tmp_path / "design.tsv") == 0
    assert mudata.written == tmp_path / "data.annotated.h5mu"

    anndata = _Container()
    monkeypatch.setattr(result_adapter, "load_converted_result", lambda _path: anndata)
    explicit = tmp_path / "explicit.h5ad"
    assert (
        cli.annotate(
            tmp_path / "data.h5ad",
            tmp_path / "design.tsv",
            explicit,
        )
        == 0
    )
    assert anndata.written == explicit


def test_fasta_rejects_no_sources_or_irrelevant_input(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from anndata_proteomics.adapters.anndata import result as result_adapter

    assert cli.fasta(tmp_path / "data.h5ad") == 1
    irrelevant = _Container()
    monkeypatch.setattr(result_adapter, "load_converted_result", lambda _path: irrelevant)
    monkeypatch.setattr(
        cli.rules_adapter,
        "require_quantification_level",
        lambda _target: "transcript",
    )
    assert (
        cli.fasta(
            tmp_path / "data.h5ad",
            tmp_path / "db.fasta",
            options=cli.FastaCliOptions(validate=False),
        )
        == 1
    )


def test_proteobench_output_guards_and_writer(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from anndata_proteomics.adapters.anndata import proteobench as proteobench_adapter
    from anndata_proteomics.adapters.anndata import result as result_adapter
    from anndata_proteomics.proteobench import config

    obj = _Container()
    monkeypatch.setattr(result_adapter, "load_converted_result", lambda _path: obj)
    monkeypatch.setattr(config, "load_module_settings", lambda _path: object())
    monkeypatch.setattr(proteobench_adapter, "resolve_targets", lambda _obj: [])
    data = tmp_path / "data.h5ad"
    data.write_text("input", encoding="utf-8")
    with pytest.raises(ValueError, match="suffix"):
        cli.proteobench(
            data,
            tmp_path / "module.toml",
            output=tmp_path / "wrong.h5mu",
        )
    with pytest.raises(ValueError, match="differ"):
        cli.proteobench(
            data,
            tmp_path / "module.toml",
            output=data,
        )
    output = tmp_path / "scored.h5ad"
    assert (
        cli.proteobench(
            data,
            tmp_path / "module.toml",
            output=output,
        )
        == 0
    )
    assert output.read_text(encoding="utf-8") == "anndata"


@pytest.mark.parametrize(("app_result", "expected"), [(None, 0), (3, 3)])
def test_main_normalizes_app_result(
    app_result: int | None,
    expected: int,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    configured: list[bool] = []
    monkeypatch.setattr(cli, "configure_default_sink", lambda: configured.append(True))
    monkeypatch.setattr(cli, "app", lambda: app_result)
    assert cli.main() == expected
    assert configured == [True]
