"""Backend-independent orchestration tests.

Every orchestrator must run on ordinary pandas/NumPy values with no AnnData or MuData
object anywhere in the test. Extraction and persistence are represented by small in-memory
functions, which is the demonstration that a second backend needs adapter work only.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from anndata_proteomics.annotation.loader import AnnotationFileOrigin, AnnotationTable
from anndata_proteomics.annotation.mulink import (
    empty_feature_mapping,
    feature_positions,
    merge_owned_feature_mapping,
)
from anndata_proteomics.converters.pipeline import NoConvertibleLevelsError, RuleSelection
from anndata_proteomics.description import (
    AnnDataDescriptionSource,
    DescriptionConversionMetadata,
    DescriptionMetadata,
    MissingProteoBenchMetadata,
    MissingQcMetadata,
    MissingQuantificationLevel,
    MissingRuleMetadata,
    MissingSearchParameters,
    MissingSoftwareName,
    MuDataDescriptionSource,
)
from anndata_proteomics.rules.schema import ParseRule
from anndata_proteomics.workflows import summary as summary_workflow
from anndata_proteomics.workflows.conversion import convert_selected_level, convert_selected_levels
from anndata_proteomics.workflows.sample_annotation import run_sample_annotation


def minimal_long_rule() -> ParseRule:
    """The smallest valid long-format rule: one obs key, one var key, one layer."""
    return ParseRule.model_validate(
        {
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
    )


def _long_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Run": ["a", "a", "b", "b"],
            "Feature": ["f1", "f2", "f1", "f2"],
            "Intensity": [1.0, 2.0, 3.0, 4.0],
        }
    )


def test_conversion_workflow_runs_without_a_storage_backend() -> None:
    rule = minimal_long_rule()
    conversion = convert_selected_level(
        _long_table(),
        rule.quantification_level,
        RuleSelection(rule, "rule_config"),
    )
    assert conversion.pieces.X.shape == (2, 2)
    assert list(conversion.pieces.obs.index) == ["a", "b"]


def test_multi_level_conversion_does_not_leak_columns_between_levels() -> None:
    """Each level converts from the caller's frame; the calculation owns any copy."""
    rule = minimal_long_rule()
    data = _long_table()
    before = list(data.columns)
    conversions = convert_selected_levels(
        data,
        {rule.quantification_level: RuleSelection(rule, "rule_config")},
    )
    assert list(data.columns) == before
    assert len(conversions) == 1


def test_conversion_workflow_rejects_an_empty_selection() -> None:
    with pytest.raises(NoConvertibleLevelsError, match="no levels supplied"):
        convert_selected_levels(_long_table(), {})


def test_sample_annotation_workflow_runs_on_plain_frames() -> None:
    observations = pd.DataFrame({"raw_file": ["a", "b"]}, index=["a", "b"])
    annotation = AnnotationTable(
        samples=pd.DataFrame({"raw_file": ["a", "b"], "condition": ["x", "y"]}),
        match_on="raw_file",
        key_field="raw_file",
    )
    result = run_sample_annotation(
        (observations,),
        annotation,
        AnnotationFileOrigin(path=Path("annotation.tsv")),
    )
    assert result.provenance.columns_added == ("condition",)
    assert list(result.observations[0].frame["condition"]) == ["x", "y"]
    assert result.diagnostics.unmatched_observation_count == 0


def test_sample_annotation_propagates_to_every_supplied_frame() -> None:
    primary = pd.DataFrame({"raw_file": ["a", "b"]}, index=["a", "b"])
    modality = pd.DataFrame({"raw_file": ["a", "b"]}, index=["a", "b"])
    annotation = AnnotationTable(
        samples=pd.DataFrame({"raw_file": ["a", "b"], "condition": ["x", "y"]}),
        match_on="raw_file",
        key_field="raw_file",
    )
    result = run_sample_annotation(
        (primary, modality),
        annotation,
        AnnotationFileOrigin(path=Path("annotation.tsv")),
    )
    assert len(result.observations) == 2
    for annotated in result.observations:
        assert list(annotated.frame["condition"]) == ["x", "y"]


@dataclass(frozen=True, slots=True)
class _InMemoryDescriptionStore:
    """A second 'backend': description sources keyed by path, no HDF5 involved."""

    anndata: dict[str, AnnDataDescriptionSource]
    mudata: dict[str, MuDataDescriptionSource]

    def readers(self) -> summary_workflow.StoredDescriptionReaders:
        return summary_workflow.StoredDescriptionReaders(
            read_anndata=lambda path: self.anndata[path.name],
            read_mudata=lambda path: self.mudata[path.name],
            read_mudata_modality=lambda path, modality: self.mudata[path.name].modalities[modality],
        )


def _empty_metadata() -> DescriptionMetadata:
    """A backend that stored no APB metadata at all still describes its shape."""
    return DescriptionMetadata(
        quantification_level=MissingQuantificationLevel(),
        software_name=MissingSoftwareName(),
        conversion=DescriptionConversionMetadata(),
        search_parameters=MissingSearchParameters(),
        rule=MissingRuleMetadata(),
        annotations={},
        qc=MissingQcMetadata(),
        proteobench=MissingProteoBenchMetadata(),
    )


def _source(n_runs: int, n_features: int) -> AnnDataDescriptionSource:
    return AnnDataDescriptionSource(
        n_runs=n_runs,
        n_features=n_features,
        layers=("intensity",),
        metadata=_empty_metadata(),
    )


def test_summary_workflow_runs_against_a_non_hdf5_backend() -> None:
    store = _InMemoryDescriptionStore(
        anndata={"one.h5ad": _source(2, 3)},
        mudata={
            "two.h5mu": MuDataDescriptionSource(
                n_runs=2,
                n_features=5,
                modalities={"ion": _source(2, 5)},
                metadata=_empty_metadata(),
            )
        },
    )
    readers = store.readers()

    described = summary_workflow.describe_path(Path("one.h5ad"), readers)
    assert described["container_type"] == "anndata"

    collection = summary_workflow.describe_path(Path("two.h5mu"), readers)
    assert collection["container_type"] == "mudata"

    modality = summary_workflow.describe_modality_path(Path("two.h5mu"), "ion", readers)
    assert modality["container_type"] == "anndata"


def test_summary_workflow_rejects_unsupported_and_mismatched_requests() -> None:
    readers = _InMemoryDescriptionStore(anndata={}, mudata={}).readers()
    with pytest.raises(ValueError, match="unsupported converted result type"):
        summary_workflow.describe_path(Path("data.parquet"), readers)
    with pytest.raises(ValueError, match="applies only to MuData"):
        summary_workflow.describe_modality_path(Path("one.h5ad"), "ion", readers)


def test_mulink_merge_preserves_another_producer_and_replaces_apb_rows() -> None:
    """The sparse merge is pure: hand-built matrices, no container."""
    positions = feature_positions(pd.Index(["p1", "p2", "prot"]))
    existing = csr_matrix(np.array([[0, 0, 1], [0, 0, 0], [0, 0, 9]], dtype=np.int8))
    owned = csr_matrix(np.array([[0, 0, 1], [0, 0, 0], [0, 0, 0]], dtype=np.int8))
    new_mapping = csr_matrix(np.array([[0, 0, 0], [0, 0, 1], [0, 0, 0]], dtype=np.int8))
    mask = np.array([True, True, False])

    merge = merge_owned_feature_mapping(existing, owned, new_mapping, mask)

    assert merge.merged[0, 2] == 0, "APB withdrew its own stale edge"
    assert merge.merged[1, 2] == 1, "APB wrote the new edge"
    assert merge.merged[2, 2] == 9, "another producer's edge survives untouched"
    assert merge.owned[1, 2] == 1
    assert positions.total == 3


def test_mulink_merge_starts_from_an_empty_stored_mapping() -> None:
    shape = (2, 2)
    empty = empty_feature_mapping(shape)
    new_mapping = csr_matrix(np.array([[0, 1], [0, 0]], dtype=np.int8))

    merge = merge_owned_feature_mapping(empty, empty, new_mapping, np.array([True, True]))

    assert merge.merged.nnz == 1
    assert merge.owned[0, 1] == 1
