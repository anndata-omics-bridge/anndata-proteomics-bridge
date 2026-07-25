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

from anndata_proteomics.converters import pipeline as conversion_pipeline
from anndata_proteomics.params.model import Parameters
from anndata_proteomics.scripts import cli


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
        effective_rules=lambda _search_parameters=None: {},
        software_name="Tool",
    )
    monkeypatch.setattr(cli, "read_table", lambda _path: frame)
    monkeypatch.setattr(cli, "load_rule_document", lambda _path: document)
    assert (
        cli.convert(
            tmp_path / "data.tsv",
            "protein",
            rule_config=tmp_path / "rule.json",
        )
        == 1
    )
    assert (
        cli.convert(
            tmp_path / "data.tsv",
            rule_config=tmp_path / "rule.json",
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
        version="2.6.0",
        version_status="present",
    )
    effective_calls: list[tuple[str, Parameters | None]] = []
    effective_rules_calls: list[Parameters | None] = []
    rule = object()

    def effective_rule(level: str, search_parameters: Parameters | None = None) -> object:
        effective_calls.append((level, search_parameters))
        return rule

    def effective_rules(
        search_parameters: Parameters | None = None,
    ) -> dict[str, object]:
        effective_rules_calls.append(search_parameters)
        return {"ion": rule}

    document = SimpleNamespace(
        levels={"ion": object()},
        effective_rule=effective_rule,
        effective_rules=effective_rules,
        software_name="DIA-NN",
    )
    monkeypatch.setattr(cli, "read_table", lambda _path: frame)
    monkeypatch.setattr(cli, "load_rule_document", lambda _path: document)
    monkeypatch.setattr(
        conversion_pipeline,
        "resolve_parameters",
        lambda *_args: resolution,
    )
    converted = ad.AnnData(np.ones((1, 1), dtype=np.float32))
    monkeypatch.setattr(cli, "_run_convert", lambda *_args, **_kwargs: converted)

    assert (
        cli.convert(
            tmp_path / "data.tsv",
            "ion",
            params=resolution.source_path,
            rule_config=tmp_path / "rule.json",
            output=tmp_path / "single",
        )
        == 0
    )
    assert effective_calls == [("ion", parameters)]

    container = _Container(modalities={"ion": converted})
    monkeypatch.setattr(
        conversion_pipeline,
        "matching_rules",
        lambda rules, _headers: rules,
    )
    monkeypatch.setattr(
        conversion_pipeline,
        "build_mudata_from_rules",
        lambda *_args, **_kwargs: container,
    )
    assert (
        cli.convert(
            tmp_path / "data.tsv",
            params=resolution.source_path,
            rule_config=tmp_path / "rule.json",
            output=tmp_path / "multi",
        )
        == 0
    )
    assert effective_rules_calls == [parameters]


def test_packaged_level_and_mudata_conversion_paths(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    frame = pd.DataFrame({"x": [1]})
    adata = ad.AnnData(np.ones((1, 1), dtype=np.float32))
    monkeypatch.setattr(cli, "read_table", lambda _path: frame)
    monkeypatch.setattr(conversion_pipeline, "recognize_software", lambda _columns: "diann")
    resolution = SimpleNamespace(
        version="2.0",
        version_status="present",
        parameters=Parameters(),
    )
    monkeypatch.setattr(
        conversion_pipeline,
        "resolve_parameters",
        lambda *_args: resolution,
    )
    monkeypatch.setattr(conversion_pipeline, "convert_level", lambda *_args, **_kwargs: adata)
    params = tmp_path / "params.txt"
    params.write_text("params", encoding="utf-8")
    output = tmp_path / "single"
    assert (
        cli.convert(
            tmp_path / "data.tsv",
            "ion",
            params=params,
            output=output,
        )
        == 0
    )
    assert output.with_suffix(".h5ad").exists()

    container = _Container(modalities={"ion": adata})
    target_parameters: list[Parameters | None] = []

    def convertible_levels(
        *_args: object,
        search_parameters: Parameters | None = None,
        **_kwargs: object,
    ) -> tuple[str, ...]:
        target_parameters.append(search_parameters)
        return ("ion",)

    monkeypatch.setattr(
        conversion_pipeline,
        "convertible_levels",
        convertible_levels,
    )
    monkeypatch.setattr(conversion_pipeline, "build_mudata", lambda *_args, **_kwargs: container)
    stale = (tmp_path / "multi").with_suffix(".h5ad")
    stale.write_text("stale", encoding="utf-8")
    assert (
        cli.convert(
            tmp_path / "data.tsv",
            params=params,
            output=tmp_path / "multi",
        )
        == 0
    )
    assert (tmp_path / "multi.h5mu").exists()
    assert not stale.exists()
    assert target_parameters == [resolution.parameters]

    monkeypatch.setattr(
        conversion_pipeline,
        "convertible_levels",
        lambda *_args, **_kwargs: (),
    )
    assert cli.convert(tmp_path / "data.tsv", params=params) == 1


def test_summary_and_annotation_writers(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from anndata_proteomics.annotation import apply, loader
    from anndata_proteomics.readers import result, summary

    monkeypatch.setattr(summary, "describe_path", lambda *_args, **_kwargs: {"n_obs": 2})
    assert cli.summary_cmd(tmp_path / "data.h5ad", json=True) == 0
    assert json.loads(capsys.readouterr().out) == {"n_obs": 2}

    annotation = pd.DataFrame({"sample": ["A"]})
    monkeypatch.setattr(loader, "load_annotation", lambda _path: annotation)
    monkeypatch.setattr(apply, "annotate_obs", lambda _obj, _annotation: None)

    mudata = _Container(modalities={})
    monkeypatch.setattr(result, "load_converted_result", lambda _path: mudata)
    assert cli.annotate(tmp_path / "data.h5mu", tmp_path / "design.tsv") == 0
    assert mudata.written == tmp_path / "data.annotated.h5mu"

    anndata = _Container()
    monkeypatch.setattr(result, "load_converted_result", lambda _path: anndata)
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
    from anndata_proteomics.readers import result

    assert cli.fasta(tmp_path / "data.h5ad") == 1
    irrelevant = _Container()
    irrelevant.uns = {"anndata_proteomics": {"quantification_level": "transcript"}}
    monkeypatch.setattr(result, "load_converted_result", lambda _path: irrelevant)
    assert (
        cli.fasta(
            tmp_path / "data.h5ad",
            tmp_path / "db.fasta",
            validate=False,
        )
        == 1
    )


def test_proteobench_output_guards_and_writer(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from anndata_proteomics.proteobench import config, pipeline
    from anndata_proteomics.readers import result

    obj = _Container()
    monkeypatch.setattr(result, "load_converted_result", lambda _path: obj)
    monkeypatch.setattr(config, "load_module_settings", lambda _path: object())
    monkeypatch.setattr(pipeline, "score_quantification", lambda *_args: obj)
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
