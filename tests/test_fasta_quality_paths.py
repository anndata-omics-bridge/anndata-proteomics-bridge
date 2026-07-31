"""Focused edge-path coverage for pure FASTA calculations and their adapter."""

from __future__ import annotations

from io import StringIO
from pathlib import Path

import anndata as ad
import mudata
import numpy as np
import pandas as pd
import pytest

from anndata_proteomics.adapters.anndata import fasta as fasta_adapter
from anndata_proteomics.annotation import validate_fasta, var_fasta
from anndata_proteomics.fasta import annotation
from anndata_proteomics.fasta.config import (
    ExplicitPatterns,
    FastaConfig,
    InferredPatterns,
    resolve_fasta_config,
)
from anndata_proteomics.fasta.parser import iter_fasta


def _adata(*, level: str, names: list[str]) -> ad.AnnData:
    result = ad.AnnData(np.zeros((1, len(names))), var=pd.DataFrame(index=names))
    result.uns["anndata_proteomics"] = {"quantification_level": level}
    return result


def test_fasta_config_reader_handles_missing_and_byte_payloads() -> None:
    empty = ad.AnnData()
    assert not fasta_adapter.has_fasta_config(empty)
    with pytest.raises(ValueError, match="no stored FASTA"):
        fasta_adapter.require_fasta_config(empty)

    resolved = resolve_fasta_config(["REV_P1"], FastaConfig())
    empty.uns["anndata_proteomics"] = {"fasta_config": resolved.model_dump_json().encode()}
    assert fasta_adapter.require_fasta_config(empty) == resolved


def test_fasta_parser_path_blank_lines_and_headerless_input(tmp_path: Path) -> None:
    path = tmp_path / "database.fasta"
    path.write_text(">P1\n\nPEP\n", encoding="utf-8")
    records = list(iter_fasta(str(path)))
    assert [(record.header, record.sequence) for record in records] == [("P1", "PEP")]
    assert list(iter_fasta(StringIO("SEQUENCE\n"))) == []


def test_fasta_annotation_empty_stream_and_before_cleavage() -> None:
    frame, _resolved = annotation.fasta_to_dataframe_with_config(
        StringIO(""),
        annotation.FastaAnnotationConfig(include_sequence=True),
    )
    assert frame.empty
    assert "sequence" in frame.columns
    frame_without_sequence, _resolved = annotation.fasta_to_dataframe_with_config(StringIO(""))
    assert "sequence" not in frame_without_sequence.columns

    stream = StringIO(">P1\nPEP\n")
    assert annotation.materialize_sources(stream) == [stream]
    assert annotation.describe_sources(stream) == ["<stream>"]

    asp_n = annotation.resolve_cleavage_name("Asp-N")
    assert annotation._find_cleavage_sites("DAD", asp_n.rule) == [2]


def test_protein_annotation_prefers_curated_duplicate_accession() -> None:
    result = var_fasta.annotate_proteins_from_fasta(
        pd.Series(["P1"], index=["protein"]),
        ">tr|P1|ONE unreviewed\nPEPTIDE\n>sp|P1|ONE curated\nPEPTIDE\n",
    )
    assert result.frame.loc["protein", "fasta.id"] == "sp|P1|ONE"


def test_empty_peptide_matching_has_a_concrete_result() -> None:
    result = validate_fasta.match_peptides_to_fasta(
        pd.Series([], index=pd.Index([], dtype="object"), dtype="object"),
        ">P1\nPEP\n",
    )
    assert result.fraction_unmatched == 0
    assert result.matches.empty
    assert list(result.matches.columns) == [
        "sequence",
        "protein_id",
        "proteinname",
        "start",
        "end",
        "length",
        "is_decoy",
        "is_contaminant",
    ]


def test_reported_protein_validation_needs_no_container() -> None:
    matches = validate_fasta.match_peptides_to_fasta(
        pd.Series(["PEP", "PEP"], index=["feature-a", "feature-b"]),
        ">sp|P1|ONE\nMPEPK\n>sp|P2|TWO\nMOTHERK\n",
    )
    validation = validate_fasta.validate_reported_proteins(
        pd.Series(["P1", "P2"], index=["feature-a", "feature-b"]),
        matches,
    )
    assert validation.summary["leading_protein_in_fasta"].tolist() == [True, True]
    assert validation.summary["peptide_in_leading_protein"].tolist() == [True, False]


def test_mudata_target_guards_are_adapter_responsibilities() -> None:
    protein = _adata(level="protein", names=["protein"])
    with mudata.set_options(pull_on_update=False):
        container = mudata.MuData({"protein": protein})
    with pytest.raises(ValueError, match="no peptide-derived"):
        fasta_adapter.read_peptide_mudata_inputs(container)
    with pytest.raises(ValueError, match="not in MuData"):
        fasta_adapter.read_peptide_mudata_modality_input(
            container,
            "missing",
        )
    with pytest.raises(ValueError, match="not peptide-derived"):
        fasta_adapter.read_peptide_mudata_modality_input(
            container,
            "protein",
        )


def test_feature_mapping_uses_only_typed_domain_inputs() -> None:
    peptide_nodes = validate_fasta.PeptideFeatureNodes(
        nodes=(validate_fasta.PeptideFeatureNode(position=0, sequence="PEP"),),
        total_nodes=2,
    )
    protein_nodes = validate_fasta.ProteinFeatureNodes(
        positions_by_accession={"P1": (1,)},
        total_nodes=2,
    )
    matches = validate_fasta.PeptideProteinMatches(
        accessions_by_sequence={"PEP": frozenset({"P1", "P2"})},
        all_accessions=frozenset({"P1", "P2"}),
    )
    result = validate_fasta.build_feature_mapping(peptide_nodes, protein_nodes, matches)
    assert result.mapping[0, 1] == 1
    assert result.represented_accessions == frozenset({"P1"})
    assert result.unrepresented_accessions == frozenset({"P2"})

    different_axis = validate_fasta.ProteinFeatureNodes(
        positions_by_accession={"P1": (1,)},
        total_nodes=3,
    )
    with pytest.raises(ValueError, match="different global feature-axis sizes"):
        validate_fasta.build_feature_mapping(peptide_nodes, different_axis, matches)


def test_pattern_policies_distinguish_inference_from_disabled() -> None:
    inferred = FastaConfig(decoy=InferredPatterns(candidates=(r"^REV_",)))
    disabled = FastaConfig(decoy=ExplicitPatterns(patterns=()))
    assert inferred.decoy.mode == "infer"
    assert disabled.decoy.mode == "explicit"
