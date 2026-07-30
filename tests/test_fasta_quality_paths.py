"""Focused edge-path coverage for FASTA parsing, annotation, and validation."""

from __future__ import annotations

from io import StringIO
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from anndata_proteomics.annotation import validate_fasta, var_fasta
from anndata_proteomics.annotation.validate_fasta import (
    FastaValidationResult,
    _TargetInput,
)
from anndata_proteomics.fasta import annotation
from anndata_proteomics.fasta.anndata_io import read_fasta_config
from anndata_proteomics.fasta.config import (
    FastaConfig,
    resolve_fasta_config,
)
from anndata_proteomics.fasta.parser import iter_fasta


def _adata(
    *,
    level: str,
    names: list[str],
    var: pd.DataFrame | None = None,
) -> ad.AnnData:
    result = ad.AnnData(
        np.zeros((1, len(names))),
        var=var if var is not None else pd.DataFrame(index=names),
    )
    result.var_names = names
    result.uns["anndata_proteomics"] = {"quantification_level": level}
    return result


def test_fasta_config_reader_handles_missing_and_byte_payloads() -> None:
    empty = SimpleNamespace(uns={})
    assert read_fasta_config(empty) is None
    empty.uns["anndata_proteomics"] = {"other": "metadata"}
    assert read_fasta_config(empty) is None

    resolved = resolve_fasta_config(["REV_P1"], FastaConfig())
    empty.uns["anndata_proteomics"]["fasta_config"] = resolved.model_dump_json().encode()
    assert read_fasta_config(empty) == resolved


def test_fasta_parser_path_blank_lines_and_headerless_input(tmp_path: Path) -> None:
    path = tmp_path / "database.fasta"
    path.write_text(">P1\n\nPEP\n", encoding="utf-8")
    records = list(iter_fasta(str(path)))
    assert [(record.header, record.sequence) for record in records] == [("P1", "PEP")]
    assert list(iter_fasta(StringIO("SEQUENCE\n"))) == []


def test_fasta_annotation_empty_conflict_stream_and_before_cleavage() -> None:
    with pytest.raises(ValueError, match="pass either"):
        annotation.fasta_to_dataframe_with_config(
            ">P1\nPEP\n",
            fasta_config=FastaConfig(),
            decoy_pattern="REV_",
        )

    frame, _resolved = annotation.fasta_to_dataframe_with_config(
        StringIO(""),
        include_sequence=True,
    )
    assert frame.empty
    assert "sequence" in frame.columns
    frame_without_sequence, _resolved = annotation.fasta_to_dataframe_with_config(StringIO(""))
    assert "sequence" not in frame_without_sequence.columns

    stream = StringIO(">P1\nPEP\n")
    assert annotation.materialize_sources(stream) == [stream]
    assert annotation.describe_sources(stream) == ["<stream>"]

    asp_n, _name = annotation.resolve_cleavage("Asp-N")
    assert annotation._find_cleavage_sites("DAD", asp_n) == [2]


def test_var_fasta_duplicate_accession_precedence_and_no_mismatch() -> None:
    frame = pd.DataFrame(
        {
            "fasta.id": ["tr|P1|ONE", "sp|P1|ONE", "custom"],
            "proteinname": ["P1", "P1", "custom"],
            "is_decoy": [False, False, False],
        }
    )
    indexed = var_fasta._index_by_join_key(frame)
    assert indexed.loc["P1", "fasta.id"] == "sp|P1|ONE"
    assert var_fasta._database_priority("tr|P1|ONE") == 1
    assert var_fasta._database_priority("custom") == 2
    var_fasta._warn_on_mismatch(
        pd.Index(["P1", "custom"]),
        np.asarray([True, True], dtype=np.bool_),
        indexed,
    )


def test_validation_empty_fraction_contract_and_target_guards() -> None:
    result = FastaValidationResult(
        summary=pd.DataFrame(),
        matches=pd.DataFrame(),
        n_features=0,
        n_unique_sequences=0,
        n_invalid_sequences=0,
        n_matched_features=0,
        n_unmatched_features=0,
        requested_backend="auto",
        backend="ahocorapy",
        sequence_field="ProForma_peptide",
        leading_protein_field=None,
        il_equivalent=False,
        fasta_config=resolve_fasta_config([], FastaConfig()),
        unmatched_sequences=[],
    )
    assert result.fraction_unmatched == 0

    with pytest.raises(ValueError, match="pass either"):
        validate_fasta._validate_targets(
            SimpleNamespace(uns={}),
            {"ion": object()},
            ">P1\nPEP\n",
            sequence_field="ProForma_peptide",
            backend="auto",
            fasta_config=FastaConfig(),
            decoy_pattern="REV_",
            contaminant_pattern=None,
            leading_protein_field=None,
            protein_match_on=None,
            il_equivalent=False,
            is_uniprot=True,
            store=False,
        )
    with pytest.raises(ValueError, match="no peptide-derived"):
        validate_fasta._validate_targets(
            SimpleNamespace(uns={}),
            {},
            ">P1\nPEP\n",
            sequence_field="ProForma_peptide",
            backend="auto",
            fasta_config=None,
            decoy_pattern=None,
            contaminant_pattern=None,
            leading_protein_field=None,
            protein_match_on=None,
            il_equivalent=False,
            is_uniprot=True,
            store=False,
        )


def test_feature_target_and_leading_protein_guards() -> None:
    ion = _adata(level="ion", names=["ion"])
    protein = _adata(level="protein", names=["protein"])
    container = SimpleNamespace(mod={"ion": ion})
    with pytest.raises(ValueError, match="not in MuData"):
        validate_fasta._resolve_feature_target(container, "missing")

    protein_container = SimpleNamespace(mod={"protein": protein})
    with pytest.raises(ValueError, match="not peptide-derived"):
        validate_fasta._resolve_feature_target(protein_container, "protein")
    with pytest.raises(ValueError, match="no peptide-derived"):
        validate_fasta._resolve_feature_target(protein_container, None)

    with pytest.raises(ValueError, match="leading_protein_field"):
        validate_fasta._resolve_leading_protein_field(ion, "missing")
    ion.var["Leading"] = [None]
    assert validate_fasta._resolve_leading_protein_field(ion, "Leading") == "Leading"
    leading = validate_fasta._feature_leading_proteins(
        ion,
        "Leading",
        is_uniprot=True,
    )
    assert leading.isna().all()


def _mapping_fixture(
    *,
    global_names: list[str],
    sequence: str | None = None,
) -> tuple[Any, dict[str, _TargetInput], pd.DataFrame]:
    protein = _adata(
        level="protein",
        names=["protein"],
        var=pd.DataFrame({"Protein_Group": ["P1"]}, index=["protein"]),
    )
    feature = _adata(level="ion", names=["feature"])
    item = _TargetInput(
        name="ion",
        target=feature,
        normalized_sequences=pd.Series([sequence], index=["feature"], dtype="object"),
        leading_proteins=pd.Series([None], index=["feature"], dtype="object"),
        leading_protein_field=None,
    )
    owner = SimpleNamespace(
        var_names=pd.Index(global_names),
        n_vars=len(global_names),
        mod={"protein": protein},
        varp={},
    )
    matches = pd.DataFrame(columns=["sequence", "proteinname"])
    return owner, {"ion": item}, matches


def test_mulink_topology_and_existing_shape_guards() -> None:
    owner, targets, matches = _mapping_fixture(global_names=["duplicate", "duplicate"])
    with pytest.raises(ValueError, match="globally unique"):
        validate_fasta._store_mulink_feature_mapping(
            owner,
            targets,
            matches,
            protein_match_on=None,
            is_uniprot=True,
        )

    owner, targets, matches = _mapping_fixture(global_names=["protein"])
    with pytest.raises(ValueError, match="absent from the MuData global"):
        validate_fasta._store_mulink_feature_mapping(
            owner,
            targets,
            matches,
            protein_match_on=None,
            is_uniprot=True,
        )

    owner, targets, matches = _mapping_fixture(
        global_names=["feature", "protein"],
        sequence=None,
    )
    stats = validate_fasta._store_mulink_feature_mapping(
        owner,
        targets,
        matches,
        protein_match_on=None,
        is_uniprot=True,
    )
    assert stats.n_fasta_edges == 0

    owner, targets, matches = _mapping_fixture(
        global_names=["feature", "protein"],
        sequence=None,
    )
    owner.varp["feature_mapping"] = csr_matrix((1, 1))
    with pytest.raises(ValueError, match="feature_mapping.*shape"):
        validate_fasta._store_mulink_feature_mapping(
            owner,
            targets,
            matches,
            protein_match_on=None,
            is_uniprot=True,
        )

    owner, targets, matches = _mapping_fixture(
        global_names=["feature", "protein"],
        sequence=None,
    )
    owner.varp["_apb_fasta_feature_mapping_contribution"] = csr_matrix((1, 1))
    with pytest.raises(ValueError, match="contribution.*shape"):
        validate_fasta._store_mulink_feature_mapping(
            owner,
            targets,
            matches,
            protein_match_on=None,
            is_uniprot=True,
        )


def test_single_pattern_distinguishes_none_from_disabled() -> None:
    assert validate_fasta._single_pattern(None) is None
    assert validate_fasta._single_pattern("") == ()
