"""Tests for lightweight APB container presentation views."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import anndata as ad
import h5py
import mudata
import numpy as np
import pandas as pd
import pytest
from mudata import MuData
from pydantic import ValidationError

from anndata_proteomics.params.anndata_io import write_search_parameters
from anndata_proteomics.params.model import Parameters
from anndata_proteomics.readers import summary as summary_module
from anndata_proteomics.readers.summary import describe, describe_path


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


def _rule(*, input_shape: str = "long") -> dict[str, object]:
    if input_shape == "long":
        intensity_source = "Intensity"
        quality_source = "Q.Value"
        obs_source = "Run"
    else:
        intensity_source = r"^Intensity_(?P<sample>.+)$"
        quality_source = r"^Q\.Value_(?P<sample>.+)$"
        obs_source = "<sample>"
    return {
        "schema_version": "0.1",
        "file_version": "test",
        "software_name": "Synthetic",
        "software_version": ".*",
        "input_shape": input_shape,
        "quantification_level": "ion",
        "axis": {
            "obs_keys": ["Run"],
            "var_keys": ["Feature"],
            "x_layer": "intensity",
        },
        "columns": {
            "obs": {"select": {"Run": obs_source}},
            "var": {
                "select": {"Feature": "Vendor.Feature", "Other": "Vendor.Other"},
                "compute": [
                    {
                        "name": "Combined",
                        "from": ["Feature", "Other"],
                        "how": "coalesce",
                    }
                ],
            },
        },
        "layers": [
            {"name": "intensity", "source": intensity_source},
            {"name": "quality", "source": quality_source},
            {
                "name": "not_materialized",
                "source": ("Absent" if input_shape == "long" else r"^Absent_(?P<sample>.+)$"),
            },
        ],
    }


def _attach_rule(obj: ad.AnnData, *, input_shape: str = "long") -> None:
    obj.uns["anndata_proteomics"]["rule_json"] = json.dumps(_rule(input_shape=input_shape))


def _reject_quantitative_dataset_reads(monkeypatch: pytest.MonkeyPatch) -> None:
    """Fail if the on-disk summary path indexes X or a named layer dataset."""
    original_getitem = h5py.Dataset.__getitem__

    def guarded_getitem(dataset: h5py.Dataset, key: Any) -> Any:
        parts = (dataset.name or "").strip("/").split("/")
        is_x = parts[0] == "X" or (len(parts) >= 3 and parts[0] == "mod" and parts[2] == "X")
        if is_x or "layers" in parts:
            raise AssertionError(f"quantitative dataset was read: {dataset.name}")
        return original_getitem(dataset, key)

    monkeypatch.setattr(h5py.Dataset, "__getitem__", guarded_getitem)


def test_describe_is_a_non_persisted_shape_and_metadata_view(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    obj = _adata()

    def layer_names_only(_obj: object) -> dict[str, object]:
        return {"intensity": object()}

    monkeypatch.setattr(summary_module, "named_layers", layer_names_only)
    result = describe(obj)

    assert result["quantification"] == {
        "n_runs": 3,
        "n_features": 4,
        "level": "peptide",
        "software_name": "Synthetic",
        "software_version": None,
        "layers": ["intensity"],
    }
    assert result["conversion"]["quantification_level"] == "peptide"
    assert "descriptive_summary" not in obj.uns["anndata_proteomics"]


def test_column_mapping_uses_the_stored_rule_without_matrix_statistics() -> None:
    obj = _adata()
    obj.layers["quality"] = np.ones(obj.shape)
    _attach_rule(obj)

    assert describe(obj)["column_mapping"] == {
        "X": {
            "layer": "intensity",
            "source": "Intensity",
            "source_kind": "column",
        },
        "layers": {
            "intensity": {"source": "Intensity", "source_kind": "column"},
            "quality": {"source": "Q.Value", "source_kind": "column"},
        },
        "obs": {"Run": "Run"},
        "var": {
            "Feature": "Vendor.Feature",
            "Other": "Vendor.Other",
            "Combined": "computed:coalesce",
        },
    }


def test_column_mapping_marks_wide_sources_as_patterns() -> None:
    obj = _adata()
    _attach_rule(obj, input_shape="wide")

    mapping = describe(obj)["column_mapping"]

    assert mapping["X"]["source"] == r"^Intensity_(?P<sample>.+)$"
    assert mapping["X"]["source_kind"] == "pattern"
    assert mapping["obs"] == {"Run": "<sample>"}


def test_malformed_stored_rule_is_not_silently_reclassified_as_legacy() -> None:
    obj = _adata()
    obj.uns["anndata_proteomics"]["rule_json"] = "{not-json"

    with pytest.raises(ValidationError):
        describe(obj)


def test_describe_decodes_canonical_enrichments() -> None:
    obj = _adata()
    namespace = obj.uns["anndata_proteomics"]
    namespace["obs_annotations_json"] = json.dumps(
        [{"source": "samples.tsv", "obs_columns_added": ["condition"], "n_obs_matched": 3}]
    )
    namespace["var_annotations_json"] = json.dumps(
        [{"source": "fasta_validation", "n_matched_features": 2}]
    )
    namespace["fasta_config"] = json.dumps({"decoy": {"patterns": ["^REV_"]}})
    namespace["qc"] = json.dumps({"schema_version": "1", "scope": "anndata"})
    namespace["proteobench"] = {
        "schema_version": "0.1",
        "scores": {"nr_feature": 4},
    }
    write_search_parameters(obj, Parameters(software_version="1.2.3", enzyme="Trypsin"))

    result = describe(obj)

    assert result["quantification"]["software_version"] == "1.2.3"
    assert result["search_parameters"]["enzyme"] == "Trypsin"
    assert result["annotations"]["obs"][0]["n_obs_matched"] == 3
    assert result["annotations"]["var"][0]["n_matched_features"] == 2
    assert result["annotations"]["fasta_config"]["decoy"]["patterns"] == ["^REV_"]
    assert result["qc"]["scope"] == "anndata"
    assert result["proteobench"]["scores"]["nr_feature"] == 4


def test_describe_path_round_trips_h5ad_without_quantitative_reads(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    obj = _adata()
    path = tmp_path / "result.h5ad"
    obj.write_h5ad(path)
    _reject_quantitative_dataset_reads(monkeypatch)

    result = describe_path(path)

    assert result["quantification"]["n_runs"] == 3
    assert result["quantification"]["layers"] == ["intensity"]


def test_describe_path_targets_mudata_without_quantitative_reads(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    peptide = _adata("pep:")
    protein = _adata("prt:")
    protein.uns["anndata_proteomics"]["quantification_level"] = "protein"
    _attach_rule(peptide)
    _attach_rule(protein)
    with mudata.set_options(pull_on_update=False):
        obj = MuData({"peptide": peptide, "protein": protein}, axis=0)
    path = tmp_path / "result.h5mu"
    obj.write_h5mu(path)
    _reject_quantitative_dataset_reads(monkeypatch)

    whole = describe_path(path)
    peptide_view = describe_path(path, modality="peptide")
    protein_view = describe_path(path, modality="protein")

    assert set(whole["modalities"]) == {"peptide", "protein"}
    assert peptide_view["container_type"] == "anndata"
    assert peptide_view["quantification"]["level"] == "peptide"
    assert peptide_view["column_mapping"]["X"]["source"] == "Intensity"
    assert protein_view["quantification"]["level"] == "protein"
    with pytest.raises(ValueError, match="not in MuData"):
        describe_path(path, modality="missing")


def test_describe_path_rejects_invalid_targets(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="only to MuData"):
        describe_path(tmp_path / "input.h5ad", modality="ion")
    with pytest.raises(ValueError, match="unsupported"):
        describe_path(tmp_path / "input.txt")
