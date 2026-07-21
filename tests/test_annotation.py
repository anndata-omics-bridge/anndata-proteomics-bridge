"""Tests for the obs sample-annotation tool (anndata_proteomics.annotation).

Warnings go through loguru → stderr; the `_loguru_to_pytest_capsys` fixture in
conftest.py wires that into pytest capture, so we read `capsys.readouterr().err`.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any

import anndata as ad
import mudata
import numpy as np
import pandas as pd
import pytest
from mudata import MuData
from pydantic import ValidationError

from anndata_proteomics.annotation.apply import annotate_obs
from anndata_proteomics.annotation.loader import load_annotation
from anndata_proteomics.annotation.schema import AnnotationSpec
from anndata_proteomics.readers.summary import describe, store_quantification_summary

RUNS = ["runA1", "runA2", "runB1", "runB2"]

_BASIC = {
    "schema_version": "0.1",
    "obs": {
        "match_on": "index",
        "key_field": "raw_file",
        "samples": [
            {"raw_file": "runA1", "sample_name": "A_rep1", "condition": "A"},
            {"raw_file": "runA2", "sample_name": "A_rep2", "condition": "A"},
            {"raw_file": "runB1", "sample_name": "B_rep1", "condition": "B"},
            {"raw_file": "runB2", "sample_name": "B_rep2", "condition": "B"},
        ],
    },
}


def _adata(var_prefix: str = "ion:", n_var: int = 3, runs: list[str] = RUNS) -> ad.AnnData:
    var_names = [f"{var_prefix}{i}" for i in range(n_var)]
    return ad.AnnData(
        X=np.arange(len(runs) * n_var, dtype="float64").reshape(len(runs), n_var),
        obs=pd.DataFrame(index=pd.Index(list(runs), name="R_FileName")),
        var=pd.DataFrame(index=pd.Index(var_names, name="ProForma_ion")),
    )


def _mudata() -> MuData:
    mods = {"ion": _adata("ion:", 3), "protein": _adata("prt:", 2)}
    with mudata.set_options(pull_on_update=False):
        return MuData(mods, axis=0)


def _spec_from(tmp_path: Path, document: dict[str, Any]) -> AnnotationSpec:
    p = tmp_path / "ann.json"
    p.write_text(json.dumps(document))
    return load_annotation(p)


# --- happy paths -------------------------------------------------------------


def test_obs_join_by_index(tmp_path: Path) -> None:
    adata = _adata()
    annotate_obs(adata, _spec_from(tmp_path, _BASIC))
    assert list(adata.obs["condition"]) == ["A", "A", "B", "B"]
    assert list(adata.obs["sample_name"]) == ["A_rep1", "A_rep2", "B_rep1", "B_rep2"]


def test_annotation_preserves_quantification_summary(tmp_path: Path) -> None:
    adata = _adata()
    adata.layers["intensity"] = adata.X.copy()
    adata.uns["anndata_proteomics"] = {
        "quantification_level": "ion",
        "software_name": "Synthetic",
    }
    store_quantification_summary(adata)
    before = describe(adata)["quantification"]

    annotate_obs(adata, _spec_from(tmp_path, _BASIC))

    assert describe(adata)["quantification"] == before


def test_join_respects_obs_order(tmp_path: Path) -> None:
    """Rows are aligned by key, not by table order."""
    adata = _adata(runs=["runB2", "runA1", "runB1", "runA2"])
    annotate_obs(adata, _spec_from(tmp_path, _BASIC))
    assert list(adata.obs["condition"]) == ["B", "A", "B", "A"]


def test_match_on_named_column(tmp_path: Path) -> None:
    adata = _adata()
    adata.obs_names = ["x0", "x1", "x2", "x3"]  # index is NOT the run name
    adata.obs["Run"] = RUNS
    document = copy.deepcopy(_BASIC)
    document["obs"]["match_on"] = "Run"
    annotate_obs(adata, _spec_from(tmp_path, document))
    assert list(adata.obs["condition"]) == ["A", "A", "B", "B"]


def test_freeform_extra_columns(tmp_path: Path) -> None:
    document = {
        "schema_version": "0.1",
        "obs": {
            "key_field": "raw_file",
            "samples": [
                {"raw_file": "runA1", "batch": 1, "genotype": "wt"},
                {"raw_file": "runA2", "batch": 2, "genotype": "ko"},
                {"raw_file": "runB1", "batch": 1, "genotype": "wt"},
                {"raw_file": "runB2", "batch": 2, "genotype": "ko"},
            ],
        },
    }
    adata = _adata()
    annotate_obs(adata, _spec_from(tmp_path, document))
    assert list(adata.obs["batch"]) == [1, 2, 1, 2]
    assert list(adata.obs["genotype"]) == ["wt", "ko", "wt", "ko"]


# --- MuData ------------------------------------------------------------------


def test_mudata_annotates_global_and_modalities(tmp_path: Path) -> None:
    md = _mudata()
    annotate_obs(md, _spec_from(tmp_path, _BASIC))
    for frame in (md.obs, md.mod["ion"].obs, md.mod["protein"].obs):
        assert list(frame["condition"]) == ["A", "A", "B", "B"]


def test_mudata_roundtrip(tmp_path: Path) -> None:
    md = _mudata()
    annotate_obs(md, _spec_from(tmp_path, _BASIC))
    out = tmp_path / "md.annotated.h5mu"
    md.write_h5mu(out)
    with mudata.set_options(pull_on_update=False):
        rt = mudata.read_h5mu(out)
    assert list(rt.obs["condition"]) == ["A", "A", "B", "B"]
    assert list(rt.mod["ion"].obs["condition"]) == ["A", "A", "B", "B"]


def test_anndata_roundtrip_records_provenance(tmp_path: Path) -> None:
    adata = _adata()
    annotate_obs(adata, _spec_from(tmp_path, _BASIC))
    out = tmp_path / "a.annotated.h5ad"
    adata.write_h5ad(out)
    rt = ad.read_h5ad(out)
    assert list(rt.obs["condition"]) == ["A", "A", "B", "B"]
    assert "obs_annotations_json" in rt.uns["anndata_proteomics"]


# --- sanitisation ------------------------------------------------------------


def test_obs_column_names_sanitised(tmp_path: Path) -> None:
    document = {
        "schema_version": "0.1",
        "obs": {
            "key_field": "raw_file",
            "samples": [
                {"raw_file": run, "Sample Name": sample}
                for run, sample in zip(
                    RUNS,
                    ["A_rep1", "A_rep2", "B_rep1", "B_rep2"],
                    strict=True,
                )
            ],
        },
    }
    adata = _adata()
    annotate_obs(adata, _spec_from(tmp_path, document))
    assert "Sample_Name" in adata.obs.columns
    assert "Sample Name" not in adata.obs.columns


def test_sanitisation_collision_raises(tmp_path: Path) -> None:
    document = {
        "schema_version": "0.1",
        "obs": {
            "key_field": "raw_file",
            "samples": [{"raw_file": run, "Sample Name": "x", "Sample-Name": "y"} for run in RUNS],
        },
    }
    adata = _adata()
    with pytest.raises(ValueError, match="collision after sanitisation"):
        annotate_obs(adata, _spec_from(tmp_path, document))


# --- mismatch / error handling ----------------------------------------------


def test_no_match_raises(tmp_path: Path) -> None:
    adata = _adata(runs=["nope1", "nope2", "nope3", "nope4"])
    with pytest.raises(ValueError, match="no obs rows matched"):
        annotate_obs(adata, _spec_from(tmp_path, _BASIC))


def test_partial_match_warns(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    # obs has an extra run with no record; the JSON has runB2 with no obs row.
    adata = _adata(runs=["runA1", "runA2", "runB1", "extra_run"])
    annotate_obs(adata, _spec_from(tmp_path, _BASIC))
    err = capsys.readouterr().err
    assert "1/4 obs rows had no matching" in err
    assert "annotation record(s) matched no obs row" in err
    # matched rows still annotated; unmatched obs row is NaN
    assert list(adata.obs["condition"])[:3] == ["A", "A", "B"]
    assert pd.isna(adata.obs["condition"].iloc[3])


def test_collision_with_existing_obs_column_raises(tmp_path: Path) -> None:
    adata = _adata()
    adata.obs["condition"] = ["pre", "pre", "pre", "pre"]
    with pytest.raises(ValueError, match="already present in obs"):
        annotate_obs(adata, _spec_from(tmp_path, _BASIC))


def test_duplicate_key_field_raises(tmp_path: Path) -> None:
    document = copy.deepcopy(_BASIC)
    document["obs"]["samples"].append(
        {"raw_file": "runA1", "sample_name": "dupe", "condition": "A"}
    )
    adata = _adata()
    with pytest.raises(ValueError, match="duplicate 'raw_file'"):
        annotate_obs(adata, _spec_from(tmp_path, document))


def test_unknown_match_on_column_raises(tmp_path: Path) -> None:
    document = copy.deepcopy(_BASIC)
    document["obs"]["match_on"] = "NoSuchColumn"
    adata = _adata()
    with pytest.raises(ValueError, match="match_on column 'NoSuchColumn' not found"):
        annotate_obs(adata, _spec_from(tmp_path, document))


# --- loader / schema ---------------------------------------------------------


def test_loader_rejects_missing_key_field(tmp_path: Path) -> None:
    document = {
        "schema_version": "0.1",
        "obs": {
            "key_field": "raw_file",
            "samples": [{"sample_name": "A_rep1", "condition": "A"}],
        },
    }
    p = tmp_path / "bad.json"
    p.write_text(json.dumps(document))
    with pytest.raises(ValidationError, match="missing key_field"):
        load_annotation(p)


def test_loader_rejects_unknown_top_level_key(tmp_path: Path) -> None:
    document = {**_BASIC, "unexpected": {"foo": "bar"}}
    p = tmp_path / "bad.json"
    p.write_text(json.dumps(document))
    with pytest.raises(ValidationError):
        load_annotation(p)


def test_loader_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        load_annotation(tmp_path / "does_not_exist.json")


# --- committed fixture -------------------------------------------------------


def test_committed_aif_fixture_applies() -> None:
    fixture = Path(__file__).parent / "data" / "annotation" / "obs_samples_AIF.json"
    spec = load_annotation(fixture)
    runs = [
        "LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_01",
        "LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_02",
        "LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_03",
        "LFQ_Orbitrap_AIF_Condition_B_Sample_Alpha_01",
        "LFQ_Orbitrap_AIF_Condition_B_Sample_Alpha_02",
        "LFQ_Orbitrap_AIF_Condition_B_Sample_Alpha_03",
    ]
    adata = _adata(runs=runs)
    annotate_obs(adata, spec)
    assert list(adata.obs["condition"]) == ["A", "A", "A", "B", "B", "B"]
