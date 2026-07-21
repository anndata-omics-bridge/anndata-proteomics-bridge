"""Tests for stage-owned APB descriptive summaries."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import mudata
import numpy as np
import pandas as pd
import pytest
from mudata import MuData

from anndata_proteomics.readers import summary as summary_module
from anndata_proteomics.readers.summary import (
    describe,
    describe_path,
    store_fasta_summary,
    store_quantification_summary,
)


def _adata(prefix: str = "") -> ad.AnnData:
    values = np.array(
        [
            [1.0, np.nan, np.nan, 0.0],
            [2.0, 5.0, np.nan, 0.0],
            [3.0, np.nan, np.nan, 0.0],
        ]
    )
    obj = ad.AnnData(
        X=values.copy(),
        obs=pd.DataFrame(index=["run1", "run2", "run3"]),
        var=pd.DataFrame(index=[f"{prefix}feature{i}" for i in range(4)]),
        layers={"intensity": values.copy()},
    )
    obj.uns["anndata_proteomics"] = {
        "quantification_level": "peptide",
        "software_name": "Synthetic",
    }
    return obj


def test_quantification_summary_has_known_missingness_and_range() -> None:
    obj = _adata()

    store_quantification_summary(obj)
    result = describe(obj)

    quantification = result["quantification"]
    assert quantification["n_runs"] == 3
    assert quantification["n_features"] == 4
    assert quantification["level"] == "peptide"
    assert quantification["software_name"] == "Synthetic"
    assert quantification["software_version"] is None
    layer = quantification["layers"]["intensity"]
    assert layer["missingness_histogram"] == [1, 1, 0, 2]
    assert layer["intensity"] == {"min": 0.0, "median": 1.0, "max": 5.0}
    stored = obj.uns["anndata_proteomics"]["descriptive_summary"]
    assert json.loads(stored) == result


def test_fasta_stage_adds_proteotypic_count_without_changing_quantification() -> None:
    obj = _adata()
    store_quantification_summary(obj)
    before = describe(obj)["quantification"]
    obj.varm["fasta_validation"] = pd.DataFrame(
        {"fasta_matching_protein_count": [1, 2, 0, 1]},
        index=obj.var_names,
    )

    store_fasta_summary(obj)
    result = describe(obj)

    assert result["quantification"] == before
    assert result["fasta"] == {"proteotypic_feature_count": 2}


def test_describe_does_not_recompute_stored_quantification(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    obj = _adata()
    store_quantification_summary(obj)

    def fail_if_called(_obj: object) -> dict:
        raise AssertionError("stored quantification should be read without matrix access")

    monkeypatch.setattr(summary_module, "_quantification_summary", fail_if_called)

    assert describe(obj)["quantification"]["n_features"] == 4


def test_describe_path_round_trips_h5ad(tmp_path: Path) -> None:
    obj = _adata()
    store_quantification_summary(obj)
    path = tmp_path / "result.h5ad"
    obj.write_h5ad(path)

    result = describe_path(path)

    assert result["quantification"]["layers"]["intensity"]["missingness_histogram"] == [1, 1, 0, 2]


def test_describe_path_targets_one_mudata_modality(tmp_path: Path) -> None:
    with mudata.set_options(pull_on_update=False):
        obj = MuData({"peptide": _adata("pep:"), "protein": _adata("prt:")}, axis=0)
    store_quantification_summary(obj)
    path = tmp_path / "result.h5mu"
    obj.write_h5mu(path)

    whole = describe_path(path)
    peptide = describe_path(path, modality="peptide")

    assert set(whole["modalities"]) == {"peptide", "protein"}
    assert peptide["container_type"] == "anndata"
    assert peptide["quantification"]["level"] == "peptide"
