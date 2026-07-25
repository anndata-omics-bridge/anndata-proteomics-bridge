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
    _numeric_summary,
    describe,
    describe_path,
    store_annotation_summary,
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
    assert layer["missingness_histogram"] == {"0": 2, "1": 0, "2": 1, "3": 1}
    assert layer["summary"] == {
        "min": 0.0,
        "first_quartile": 0.0,
        "median": 1.0,
        "mean": pytest.approx(11 / 7),
        "third_quartile": 2.5,
        "max": 5.0,
    }
    stored = obj.uns["anndata_proteomics"]["descriptive_summary"]
    assert json.loads(stored) == result


def test_quantification_summary_records_applied_column_mapping() -> None:
    obj = _adata()
    obj.layers["quality"] = np.ones(obj.shape)
    _attach_rule(obj)

    store_quantification_summary(obj)

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


def test_column_mapping_marks_wide_layer_sources_as_patterns() -> None:
    obj = _adata()
    _attach_rule(obj, input_shape="wide")

    store_quantification_summary(obj)

    mapping = describe(obj)["column_mapping"]
    assert mapping["X"] == {
        "layer": "intensity",
        "source": r"^Intensity_(?P<sample>.+)$",
        "source_kind": "pattern",
    }
    assert mapping["layers"] == {
        "intensity": {
            "source": r"^Intensity_(?P<sample>.+)$",
            "source_kind": "pattern",
        }
    }
    assert mapping["obs"] == {"Run": "<sample>"}


@pytest.mark.parametrize(
    "rule_json",
    [None, "{not-json", {}, b"\xff"],
)
def test_column_mapping_tolerates_legacy_or_malformed_rule_json(
    rule_json: object,
) -> None:
    obj = _adata()
    if rule_json is not None:
        obj.uns["anndata_proteomics"]["rule_json"] = rule_json

    store_quantification_summary(obj)

    assert "column_mapping" not in describe(obj)


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
    assert result["fasta"] == {
        "feature_count": 4,
        "matched_feature_count": 3,
        "proteotypic_feature_count": 2,
    }


def test_protein_fasta_summary_counts_annotated_features() -> None:
    obj = _adata()
    obj.uns["anndata_proteomics"]["quantification_level"] = "protein"
    obj.varm["fasta"] = pd.DataFrame(
        {
            "fasta.id": pd.Categorical(["P1", "P2", None, "P4"]),
            "gene_name": pd.Categorical(["G1", "G2", None, "G4"]),
        },
        index=obj.var_names,
    )

    store_fasta_summary(obj)

    assert describe(obj)["fasta"] == {
        "feature_count": 4,
        "annotated_feature_count": 3,
    }


def test_describe_does_not_recompute_stored_quantification(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    obj = _adata()
    store_quantification_summary(obj)

    def fail_if_called(_obj: object) -> dict[str, object]:
        raise AssertionError("stored quantification should be read without matrix access")

    monkeypatch.setattr(summary_module, "_quantification_summary", fail_if_called)

    assert describe(obj)["quantification"]["n_features"] == 4


def test_describe_path_round_trips_h5ad(tmp_path: Path) -> None:
    obj = _adata()
    store_quantification_summary(obj)
    path = tmp_path / "result.h5ad"
    obj.write_h5ad(path)

    result = describe_path(path)

    assert result["quantification"]["layers"]["intensity"]["missingness_histogram"] == {
        "0": 2,
        "1": 0,
        "2": 1,
        "3": 1,
    }


def test_describe_upgrades_v1_present_run_histogram() -> None:
    obj = _adata()
    obj.uns["anndata_proteomics"]["descriptive_summary"] = json.dumps(
        {
            "schema_version": "1",
            "container_type": "anndata",
            "quantification": {
                "n_runs": 3,
                "layers": {
                    "intensity": {
                        "missingness_histogram": [1, 1, 0, 2],
                        "intensity": {"min": 0.0, "median": 1.0, "max": 5.0},
                    }
                },
            },
        }
    )

    result = describe(obj)

    assert result["schema_version"] == "5"
    layer = result["quantification"]["layers"]["intensity"]
    assert layer["missingness_histogram"] == {
        "0": 2,
        "1": 0,
        "2": 1,
        "3": 1,
    }
    assert layer["summary"] == {
        "min": 0.0,
        "first_quartile": None,
        "median": 1.0,
        "mean": None,
        "third_quartile": None,
        "max": 5.0,
    }


def test_numeric_summary_matches_r_summary_quartiles() -> None:
    assert _numeric_summary(np.arange(1.0, 7.0)) == {
        "min": 1.0,
        "first_quartile": 2.25,
        "median": 3.5,
        "mean": 3.5,
        "third_quartile": 4.75,
        "max": 6.0,
    }


def test_describe_path_targets_one_mudata_modality(tmp_path: Path) -> None:
    peptide = _adata("pep:")
    protein = _adata("prt:")
    protein.uns["anndata_proteomics"]["quantification_level"] = "protein"
    peptide.obs["condition"] = ["A", "A", "B"]
    protein.obs["condition"] = ["A", "A", "B"]
    peptide.varm["fasta_validation"] = pd.DataFrame(
        {
            "peptide_in_fasta": [True, True, False, True],
            "fasta_matching_protein_count": [1, 2, 0, 1],
        },
        index=peptide.var_names,
    )
    protein.varm["fasta"] = pd.DataFrame(
        {"fasta.id": pd.Categorical(["P1", "P2", None, "P4"])},
        index=protein.var_names,
    )
    _attach_rule(peptide)
    _attach_rule(protein)
    with mudata.set_options(pull_on_update=False):
        obj = MuData({"peptide": peptide, "protein": protein}, axis=0)
    obj.obs["condition"] = ["A", "A", "B"]
    store_quantification_summary(obj)
    store_annotation_summary(obj, fields=["condition"], n_annotated_runs=3)
    store_fasta_summary(obj)
    path = tmp_path / "result.h5mu"
    obj.write_h5mu(path)

    whole = describe_path(path)
    peptide = describe_path(path, modality="peptide")
    protein = describe_path(path, modality="protein")

    assert set(whole["modalities"]) == {"peptide", "protein"}
    assert peptide["container_type"] == "anndata"
    assert peptide["quantification"]["level"] == "peptide"
    assert peptide["column_mapping"]["X"]["source"] == "Intensity"
    assert peptide["annotation"]["group_counts"] == {"condition": 2}
    assert peptide["fasta"]["matched_feature_count"] == 3
    assert protein["annotation"]["group_counts"] == {"condition": 2}
    assert protein["column_mapping"]["X"]["source"] == "Intensity"
    assert protein["fasta"]["annotated_feature_count"] == 3
