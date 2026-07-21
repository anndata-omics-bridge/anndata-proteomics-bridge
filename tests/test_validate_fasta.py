"""Tests for FASTA identification validation (Aho-Corasick).

Covers the test matrix in TODO/TODO_aho_cor.md: unmodified/modified sequences,
shared and repeated matches, unmatched flagging, multi-file databases, decoy
policy, normalization edge cases, backend equivalence, h5ad/h5mu round-trips,
and automatic validation through the ``apb fasta`` CLI.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import mudata
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from anndata_proteomics.annotation.validate_fasta import (
    FastaValidationResult,
    _replace_owned_feature_mapping,
    validate_peptide_modalities_against_fasta,
    validate_peptides_against_fasta,
)
from anndata_proteomics.fasta.anndata_io import read_fasta_config
from anndata_proteomics.fasta.config import FastaConfig
from prozor.ahocorasick import get_available_backends

FASTA = (
    ">sp|P12345|PROT1 first protein\n"
    "MKWVTGVFRRDTHKGVFRR\n"  # GVFRR occurs twice; contains DTHK
    ">sp|P67890|PROT2 second protein\n"
    "MRGVFRRSEQXUNIQUEPEPX\n"  # GVFRR (shared) and UNIQUEPEP
    ">REV_sp|P00000|DECOY decoy\n"
    "GVFRRDECOYONLY\n"  # decoy carrying GVFRR + a decoy-only peptide
)


def _peptide_adata(sequences, index=None, level="ion", extra=None):
    """A peptide-derived AnnData with a ``ProForma_peptide`` var column."""
    index = index or [f"feat{i}" for i in range(len(sequences))]
    var = pd.DataFrame({"ProForma_peptide": sequences}, index=index)
    if extra:
        for col, vals in extra.items():
            var[col] = vals
    adata = ad.AnnData(X=np.zeros((1, len(sequences)), dtype=float), var=var)
    adata.uns["anndata_proteomics"] = {"quantification_level": level, "schema_version": "1"}
    return adata


# --- basic presence -------------------------------------------------------


def test_returns_result_and_mutates():
    a = _peptide_adata(["GVFRR", "DTHK"])
    result = validate_peptides_against_fasta(a, FASTA)
    assert isinstance(result, FastaValidationResult)
    assert "fasta_validation" in a.varm
    assert read_fasta_config(a) is not None
    assert "fasta_matches" not in a.uns["anndata_proteomics"]


def test_every_unmodified_peptide_found():
    a = _peptide_adata(["DTHK", "UNIQUEPEP"])
    validate_peptides_against_fasta(a, FASTA)
    assert a.varm["fasta_validation"]["peptide_in_fasta"].all()


def test_modified_peptidoform_validates_via_unmodified_sequence():
    """The modified peptidoform column is irrelevant; validation uses ProForma_peptide."""
    a = _peptide_adata(
        ["DTHK"],
        extra={"ProForma_peptidoform": ["DT[Phospho]HK"]},
    )
    result = validate_peptides_against_fasta(a, FASTA)
    assert result.n_matched_features == 1


# --- shared / repeated ----------------------------------------------------


def test_shared_peptide_maps_to_all_proteins():
    a = _peptide_adata(["GVFRR"])
    validate_peptides_against_fasta(a, FASTA)
    proteins = a.varm["fasta_validation"].loc["feat0", "fasta_matching_protein_ids"].split(";")
    assert set(proteins) == {"P12345", "P67890", "REV_P00000"}
    assert a.varm["fasta_validation"].loc["feat0", "fasta_match_site_count"] == 4
    assert a.varm["fasta_validation"].loc["feat0", "fasta_matching_protein_count"] == 3


def test_repeated_occurrence_preserves_both_coordinates():
    a = _peptide_adata(["GVFRR"])
    result = validate_peptides_against_fasta(a, FASTA)
    prot1 = result.matches[result.matches["protein_id"] == "sp|P12345|PROT1"]
    assert sorted(prot1["start"].tolist()) == [5, 14]


def test_ion_and_fragment_share_one_match_result():
    """Two features with the same peptide reuse one entry in the unique-keyed table."""
    a = _peptide_adata(["DTHK", "DTHK"], index=["ion", "frag"])
    result = validate_peptides_against_fasta(a, FASTA)
    # one unique sequence -> one match row (DTHK occurs once in the DB target)
    assert (result.matches["sequence"] == "DTHK").sum() == 1
    counts = a.varm["fasta_validation"]["fasta_match_site_count"]
    assert counts["ion"] == counts["frag"] == 1


# --- unmatched ------------------------------------------------------------


def test_unmatched_peptide_flagged_and_reported():
    a = _peptide_adata(["DTHK", "ZZZZZ"])
    result = validate_peptides_against_fasta(a, FASTA)
    assert not a.varm["fasta_validation"].loc["feat1", "peptide_in_fasta"]
    assert "ZZZZZ" in result.sample_unmatched()
    assert result.n_unmatched_features == 1


# --- multi-file database --------------------------------------------------


def test_multiple_fasta_files_form_one_database(tmp_path: Path):
    f1 = tmp_path / "a.fasta"
    f2 = tmp_path / "b.fasta"
    f1.write_text(">sp|A|ONE\nMKDTHKXX\n")
    f2.write_text(">sp|B|TWO\nMRUNIQUEPEPX\n")
    a = _peptide_adata(["DTHK", "UNIQUEPEP"])
    result = validate_peptides_against_fasta(a, [f1, f2])
    assert result.n_matched_features == 2
    assert set(result.matches["proteinname"]) == {"A", "B"}


def test_duplicate_fasta_ids_do_not_hide_later_sequences() -> None:
    fasta = ">sp|P1|ONE\nFIRSTPEP\n>sp|P1|ONE\nSECONDPEP\n"
    a = _peptide_adata(["SECONDPEP"])

    result = validate_peptides_against_fasta(a, fasta)

    assert result.summary.loc["feat0", "peptide_in_fasta"]
    assert result.fasta_config.n_fasta_ids == 2
    assert result.matches["protein_id"].tolist() == ["sp|P1|ONE"]


def test_non_uniprot_proteinname_is_full_id():
    """With is_uniprot=False, proteinname is the full header id (no middle-field extraction)."""
    a = _peptide_adata(["DTHK"])
    result = validate_peptides_against_fasta(a, FASTA, is_uniprot=False, store=False)
    row = result.matches[result.matches["sequence"] == "DTHK"].iloc[0]
    assert row["protein_id"] == "sp|P12345|PROT1"
    assert row["proteinname"] == "sp|P12345|PROT1"  # not "P12345"


# --- decoy / contaminant policy ------------------------------------------


def test_decoy_only_match_is_valid_and_classified():
    a = _peptide_adata(["DECOYONLY"])
    result = validate_peptides_against_fasta(a, FASTA)
    assert a.varm["fasta_validation"].loc["feat0", "peptide_in_fasta"]
    assert result.matches.loc[result.matches["sequence"] == "DECOYONLY", "is_decoy"].all()


def test_contaminant_only_match_is_valid_and_classified():
    a = _peptide_adata(["CONTAMINANT"])
    result = validate_peptides_against_fasta(
        a,
        ">zz|C1|CONT contaminant\nXXCONTAMINANTXX\n",
    )
    assert a.varm["fasta_validation"].loc["feat0", "peptide_in_fasta"]
    assert result.matches["is_contaminant"].all()


def test_custom_candidate_configuration_is_inferred_during_validation() -> None:
    a = _peptide_adata(["CUSTOMPEP"])
    result = validate_peptides_against_fasta(
        a,
        ">CUSTOM_P1 custom decoy\nXXCUSTOMPEPXX\n",
        fasta_config=FastaConfig(decoy_candidates=(r"^CUSTOM_",)),
    )
    assert result.fasta_config.decoy.patterns == (r"^CUSTOM_",)
    assert result.matches["is_decoy"].all()


def test_invalid_backend_raises():
    a = _peptide_adata(["DTHK"])
    with pytest.raises(ValueError, match="backend"):
        validate_peptides_against_fasta(a, FASTA, backend="typo")


# --- normalization edge cases --------------------------------------------


def test_empty_input():
    a = _peptide_adata([], index=[])
    result = validate_peptides_against_fasta(a, FASTA)
    assert result.n_features == 0
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


def test_duplicate_peptides_deduplicated_in_matches():
    a = _peptide_adata(["GVFRR", "GVFRR"])
    result = validate_peptides_against_fasta(a, FASTA)
    # unique sequence scanned once; both features validate
    assert a.varm["fasta_validation"]["peptide_in_fasta"].all()
    assert result.n_unique_sequences == 1


def test_missing_sequence_flagged_invalid():
    a = _peptide_adata(["DTHK", None])
    result = validate_peptides_against_fasta(a, FASTA)
    assert result.n_invalid_sequences == 1
    assert not a.varm["fasta_validation"].loc["feat1", "peptide_in_fasta"]


def test_lowercase_sequence_normalized():
    a = _peptide_adata(["dthk"])
    validate_peptides_against_fasta(a, FASTA)
    assert a.varm["fasta_validation"].loc["feat0", "peptide_in_fasta"]


def test_il_not_equivalent_by_default():
    """Peptide with I must NOT match its L-variant in the DB when il_equivalent=False."""
    a = _peptide_adata(["PEPTIDE"])  # has I
    validate_peptides_against_fasta(a, ">P\nXPEPTLDEX\n")  # DB has L; default strict
    assert not a.varm["fasta_validation"].loc["feat0", "peptide_in_fasta"]


def test_il_equivalent_opt_in():
    a = _peptide_adata(["PEPTIDE"])  # I
    validate_peptides_against_fasta(a, ">P\nXPEPTLDEX\n", il_equivalent=True)  # L
    assert a.varm["fasta_validation"].loc["feat0", "peptide_in_fasta"]


def test_ambiguous_residue_matched_literally():
    a = _peptide_adata(["ABXZ"])
    validate_peptides_against_fasta(a, ">P\nMKABXZQ\n")
    assert a.varm["fasta_validation"].loc["feat0", "peptide_in_fasta"]


def test_ambiguous_residue_is_not_a_wildcard():
    """X is a literal residue: 'ABXZ' must NOT match 'ABYZ'."""
    a = _peptide_adata(["ABXZ"])
    validate_peptides_against_fasta(a, ">P\nMKABYZQ\n")
    assert not a.varm["fasta_validation"].loc["feat0", "peptide_in_fasta"]


def test_empty_string_sequence_flagged_invalid():
    a = _peptide_adata(["DTHK", "   "])  # whitespace-only normalizes to invalid
    result = validate_peptides_against_fasta(a, FASTA)
    assert result.n_invalid_sequences == 1
    assert not a.varm["fasta_validation"].loc["feat1", "peptide_in_fasta"]


def test_missing_sequence_field_raises():
    a = _peptide_adata(["DTHK"])
    with pytest.raises(ValueError, match="sequence_field"):
        validate_peptides_against_fasta(a, FASTA, sequence_field="NotThere")


# --- target resolution ----------------------------------------------------


def test_protein_level_rejected():
    a = _peptide_adata(["P12345"], level="protein")
    with pytest.raises(ValueError, match="peptide-derived"):
        validate_peptides_against_fasta(a, FASTA)


def test_mudata_auto_picks_single_peptide_modality():
    pep = _peptide_adata(["DTHK", "GVFRR"], level="peptidoform")
    prot = ad.AnnData(
        X=np.zeros((1, 1)), var=pd.DataFrame({"Protein_Group": ["P12345"]}, index=["P12345"])
    )
    prot.uns["anndata_proteomics"] = {"quantification_level": "protein"}
    md = mudata.MuData({"peptidoform": pep, "protein": prot})
    result = validate_peptides_against_fasta(md, FASTA)
    assert result.n_matched_features == 2
    assert "fasta_validation" in md.mod["peptidoform"].varm


def test_mudata_ambiguous_requires_modality():
    pep = _peptide_adata(["DTHK"], level="peptidoform")
    ion = _peptide_adata(["GVFRR"], level="ion")
    md = mudata.MuData({"peptidoform": pep, "ion": ion})
    with pytest.raises(ValueError, match="multiple peptide-derived"):
        validate_peptides_against_fasta(md, FASTA)
    # explicit modality resolves it
    validate_peptides_against_fasta(md, FASTA, modality="ion")
    assert "fasta_validation" in md.mod["ion"].varm


def test_leading_protein_validations_are_independent() -> None:
    a = _peptide_adata(
        ["DTHK", "DTHK", "DTHK"],
        extra={"Protein_Group": ["P12345", "P67890", "NOT_IN_FASTA"]},
    )
    validate_peptides_against_fasta(a, FASTA)
    validation = a.varm["fasta_validation"]
    assert validation["leading_protein_in_fasta"].tolist() == [True, True, False]
    assert validation["peptide_in_leading_protein"].tolist() == [True, False, False]


def test_all_modalities_validate_and_store_mulink_edges(tmp_path: Path) -> None:
    ion = _peptide_adata(
        ["DTHK", "GVFRR"],
        index=["ion:dthk", "ion:gvfrr"],
        level="ion",
        extra={"Protein_Group": ["P12345", "P67890"]},
    )
    fragment = _peptide_adata(
        ["DTHK"],
        index=["frg:dthk"],
        level="fragment",
        extra={"Protein_Group": ["P12345"]},
    )
    protein = ad.AnnData(
        X=np.zeros((1, 2)),
        var=pd.DataFrame(
            {"Protein_Group": ["P12345", "P67890"]},
            index=["prt:P12345", "prt:P67890"],
        ),
    )
    protein.uns["anndata_proteomics"] = {"quantification_level": "protein"}
    with mudata.set_options(pull_on_update=False):
        md = mudata.MuData({"ion": ion, "fragment": fragment, "protein": protein}, axis=0)

    results = validate_peptide_modalities_against_fasta(md, FASTA)

    assert set(results) == {"ion", "fragment"}
    assert "fasta_validation" in md.mod["ion"].varm
    assert "fasta_validation" in md.mod["fragment"].varm
    assert set(results["ion"].matches["sequence"]) == {"DTHK", "GVFRR"}
    assert set(results["fragment"].matches["sequence"]) == {"DTHK"}
    mapping = md.varp["feature_mapping"]
    positions = {name: i for i, name in enumerate(md.var_names)}
    assert mapping[positions["ion:dthk"], positions["prt:P12345"]] == 1
    assert mapping[positions["frg:dthk"], positions["prt:P12345"]] == 1
    assert mapping[positions["ion:gvfrr"], positions["prt:P67890"]] == 1
    assert _last_provenance(md.mod["ion"])["n_unrepresented_fasta_proteins"] == 1

    path = tmp_path / "linked.h5mu"
    md.write_h5mu(path)
    with mudata.set_options(pull_on_update=False):
        restored = mudata.read_h5mu(path)
    assert (restored.varp["feature_mapping"] != mapping).nnz == 0
    assert (
        restored.varp["_apb_fasta_feature_mapping_contribution"]
        != md.varp["_apb_fasta_feature_mapping_contribution"]
    ).nnz == 0


def test_mulink_update_preserves_existing_weights_and_is_idempotent() -> None:
    ion = _peptide_adata(
        ["DTHK"],
        index=["ion:dthk"],
        extra={"Protein_Group": ["P12345"]},
    )
    protein = ad.AnnData(
        X=np.zeros((1, 1)),
        var=pd.DataFrame(
            {"Protein_Group": ["P12345"]},
            index=["prt:P12345"],
        ),
    )
    protein.uns["anndata_proteomics"] = {"quantification_level": "protein"}
    with mudata.set_options(pull_on_update=False):
        md = mudata.MuData({"ion": ion, "protein": protein}, axis=0)
    positions = {name: i for i, name in enumerate(md.var_names)}
    existing = csr_matrix(
        (
            [7.0],
            (
                [positions["prt:P12345"]],
                [positions["ion:dthk"]],
            ),
        ),
        shape=(md.n_vars, md.n_vars),
    )
    md.varp["feature_mapping"] = existing

    validate_peptide_modalities_against_fasta(md, FASTA)
    once = md.varp["feature_mapping"].copy()
    validate_peptide_modalities_against_fasta(md, FASTA)

    assert md.varp["feature_mapping"][positions["prt:P12345"], positions["ion:dthk"]] == 7.0
    assert md.varp["feature_mapping"][positions["ion:dthk"], positions["prt:P12345"]] == 1.0
    assert (md.varp["feature_mapping"] != once).nnz == 0


def test_mulink_revalidation_replaces_only_apb_owned_edges() -> None:
    ion = _peptide_adata(["PEPTIDE"], index=["ion:peptide"])
    protein = ad.AnnData(
        X=np.zeros((1, 2)),
        var=pd.DataFrame(
            {"Protein_Group": ["P1", "P2"]},
            index=["prt:P1", "prt:P2"],
        ),
    )
    protein.uns["anndata_proteomics"] = {"quantification_level": "protein"}
    with mudata.set_options(pull_on_update=False):
        md = mudata.MuData({"ion": ion, "protein": protein}, axis=0)
    positions = {name: i for i, name in enumerate(md.var_names)}
    fasta_p1 = ">sp|P1|ONE\nMPEPTIDEK\n>sp|P2|TWO\nMQQQQQK\n"
    fasta_p2 = ">sp|P1|ONE\nMQQQQQK\n>sp|P2|TWO\nMPEPTIDEK\n"

    validate_peptide_modalities_against_fasta(md, fasta_p1)
    assert md.varp["feature_mapping"][positions["ion:peptide"], positions["prt:P1"]] == 1
    assert md.varp["feature_mapping"][positions["ion:peptide"], positions["prt:P2"]] == 0

    validate_peptide_modalities_against_fasta(md, fasta_p2)
    mapping = md.varp["feature_mapping"]
    owned = md.varp["_apb_fasta_feature_mapping_contribution"]
    assert mapping[positions["ion:peptide"], positions["prt:P1"]] == 0
    assert mapping[positions["ion:peptide"], positions["prt:P2"]] == 1
    assert owned[positions["ion:peptide"], positions["prt:P1"]] == 0
    assert owned[positions["ion:peptide"], positions["prt:P2"]] == 1
    assert (
        md.mod["ion"].varm["fasta_validation"].loc["ion:peptide", "fasta_matching_protein_ids"]
        == "P2"
    )


def test_mulink_does_not_claim_or_remove_an_existing_overlapping_edge() -> None:
    ion = _peptide_adata(["PEPTIDE"], index=["ion:peptide"])
    protein = ad.AnnData(
        X=np.zeros((1, 2)),
        var=pd.DataFrame(
            {"Protein_Group": ["P1", "P2"]},
            index=["prt:P1", "prt:P2"],
        ),
    )
    protein.uns["anndata_proteomics"] = {"quantification_level": "protein"}
    with mudata.set_options(pull_on_update=False):
        md = mudata.MuData({"ion": ion, "protein": protein}, axis=0)
    positions = {name: i for i, name in enumerate(md.var_names)}
    row = positions["ion:peptide"]
    p1 = positions["prt:P1"]
    p2 = positions["prt:P2"]
    md.varp["feature_mapping"] = csr_matrix(
        ([7.0], ([row], [p1])),
        shape=(md.n_vars, md.n_vars),
    )
    fasta_p1 = ">sp|P1|ONE\nMPEPTIDEK\n>sp|P2|TWO\nMQQQQQK\n"
    fasta_p2 = ">sp|P1|ONE\nMQQQQQK\n>sp|P2|TWO\nMPEPTIDEK\n"

    validate_peptide_modalities_against_fasta(md, fasta_p1)
    assert md.varp["feature_mapping"][row, p1] == 7.0
    assert md.varp["_apb_fasta_feature_mapping_contribution"][row, p1] == 0

    validate_peptide_modalities_against_fasta(md, fasta_p2)
    assert md.varp["feature_mapping"][row, p1] == 7.0
    assert md.varp["feature_mapping"][row, p2] == 1.0


def test_single_modality_revalidation_retains_other_apb_owned_edges() -> None:
    ion = _peptide_adata(["PEPTIDE"], index=["ion:peptide"], level="ion")
    fragment = _peptide_adata(
        ["PEPTIDE"],
        index=["frg:peptide"],
        level="fragment",
    )
    protein = ad.AnnData(
        X=np.zeros((1, 2)),
        var=pd.DataFrame(
            {"Protein_Group": ["P1", "P2"]},
            index=["prt:P1", "prt:P2"],
        ),
    )
    protein.uns["anndata_proteomics"] = {"quantification_level": "protein"}
    with mudata.set_options(pull_on_update=False):
        md = mudata.MuData(
            {"ion": ion, "fragment": fragment, "protein": protein},
            axis=0,
        )
    positions = {name: i for i, name in enumerate(md.var_names)}
    fasta_p1 = ">sp|P1|ONE\nMPEPTIDEK\n>sp|P2|TWO\nMQQQQQK\n"
    fasta_p2 = ">sp|P1|ONE\nMQQQQQK\n>sp|P2|TWO\nMPEPTIDEK\n"

    validate_peptide_modalities_against_fasta(md, fasta_p1)
    validate_peptides_against_fasta(md, fasta_p2, modality="ion")

    mapping = md.varp["feature_mapping"]
    assert mapping[positions["ion:peptide"], positions["prt:P1"]] == 0
    assert mapping[positions["ion:peptide"], positions["prt:P2"]] == 1
    assert mapping[positions["frg:peptide"], positions["prt:P1"]] == 1
    assert mapping[positions["frg:peptide"], positions["prt:P2"]] == 0


def test_mulink_preserves_uint64_dtype_and_large_weight_exactly() -> None:
    ion = _peptide_adata(["DTHK"], index=["ion:dthk"])
    protein = ad.AnnData(
        X=np.zeros((1, 1)),
        var=pd.DataFrame(
            {"Protein_Group": ["P12345"]},
            index=["prt:P12345"],
        ),
    )
    protein.uns["anndata_proteomics"] = {"quantification_level": "protein"}
    with mudata.set_options(pull_on_update=False):
        md = mudata.MuData({"ion": ion, "protein": protein}, axis=0)
    positions = {name: i for i, name in enumerate(md.var_names)}
    large_weight = 2**63 + 1
    md.varp["feature_mapping"] = csr_matrix(
        (
            np.array([large_weight], dtype=np.uint64),
            (
                [positions["prt:P12345"]],
                [positions["ion:dthk"]],
            ),
        ),
        shape=(md.n_vars, md.n_vars),
        dtype=np.uint64,
    )

    validate_peptide_modalities_against_fasta(md, FASTA)

    mapping = md.varp["feature_mapping"]
    assert mapping.dtype == np.dtype("uint64")
    assert mapping[positions["prt:P12345"], positions["ion:dthk"]] == large_weight
    assert mapping[positions["ion:dthk"], positions["prt:P12345"]] == 1


def test_owned_edge_replacement_uses_overflow_safe_sparse_coordinates() -> None:
    size = 70_000
    owned_row, owned_col = 61_356, 65_000
    unrelated_row, unrelated_col = 0, 17_704
    existing = csr_matrix(
        (
            [1, 7],
            ([owned_row, unrelated_row], [owned_col, unrelated_col]),
        ),
        shape=(size, size),
    )
    old_owned = csr_matrix(
        ([1], ([owned_row], [owned_col])),
        shape=(size, size),
        dtype="int8",
    )
    target_rows = np.zeros(size, dtype=bool)
    target_rows[owned_row] = True

    merged, owned = _replace_owned_feature_mapping(
        existing,
        old_owned,
        csr_matrix((size, size), dtype="int8"),
        target_rows,
    )

    assert merged[owned_row, owned_col] == 0
    assert merged[unrelated_row, unrelated_col] == 7
    assert owned.nnz == 0


# --- backend equivalence --------------------------------------------------


# A protein with NESTED peptides (a fully-cleaved peptide and its missed-cleavage
# extension share a start) plus a repeated peptide — the case where a
# non-overlapping Aho-Corasick backend would silently drop matches.
OVERLAP_FASTA = ">sp|PN|NEST\nMKSAMPLERPEPTIDERKGVFRRXGVFRR\n"
OVERLAP_PEPTIDES = ["SAMPLER", "SAMPLERPEPTIDER", "GVFRR", "ZZZZZ"]


def _sorted_matches(result):
    return result.matches.sort_values(["sequence", "protein_id", "start"]).reset_index(drop=True)


@pytest.mark.parametrize("backend", get_available_backends())
def test_backends_produce_identical_match_tables(backend):
    """Every available backend agrees with ahocorapy, incl. nested/overlapping peptides."""
    result = validate_peptides_against_fasta(
        _peptide_adata(OVERLAP_PEPTIDES), OVERLAP_FASTA, backend=backend, store=False
    )
    ref = validate_peptides_against_fasta(
        _peptide_adata(OVERLAP_PEPTIDES), OVERLAP_FASTA, backend="ahocorapy", store=False
    )
    pd.testing.assert_frame_equal(_sorted_matches(result), _sorted_matches(ref))
    # guard: the fixture must actually exercise the overlapping path
    assert "SAMPLERPEPTIDER" in set(result.matches["sequence"])


@pytest.mark.parametrize("backend", get_available_backends())
def test_nested_peptide_matches_on_every_backend(backend):
    """A missed-cleavage extension nested in a protein validates on all backends."""
    a = _peptide_adata(["SAMPLER", "SAMPLERPEPTIDER"], index=["short", "long"])
    validate_peptides_against_fasta(a, OVERLAP_FASTA, backend=backend)
    assert a.varm["fasta_validation"]["peptide_in_fasta"].all()


def test_match_table_absolute_coordinates():
    """Pin exact start/end (half-open)/length for a repeated peptide."""
    a = _peptide_adata(["GVFRR"])
    result = validate_peptides_against_fasta(a, FASTA, store=False)
    prot1 = (
        result.matches[result.matches["protein_id"] == "sp|P12345|PROT1"]
        .sort_values("start")
        .reset_index(drop=True)
    )
    assert prot1["start"].tolist() == [5, 14]
    assert prot1["end"].tolist() == [10, 19]
    assert prot1["length"].tolist() == [5, 5]


# --- round-trips ----------------------------------------------------------


def test_h5ad_round_trip(tmp_path: Path):
    a = _peptide_adata(["GVFRR", "DTHK", "ZZZZZ"])
    validate_peptides_against_fasta(a, FASTA)
    p = tmp_path / "r.h5ad"
    a.write_h5ad(p)
    b = ad.read_h5ad(p)
    pd.testing.assert_frame_equal(a.varm["fasta_validation"], b.varm["fasta_validation"])
    assert read_fasta_config(b) == read_fasta_config(a)


def test_h5mu_round_trip(tmp_path: Path):
    pep = _peptide_adata(["GVFRR", "DTHK"], level="peptidoform")
    md = mudata.MuData({"peptidoform": pep})
    validate_peptides_against_fasta(md, FASTA)
    p = tmp_path / "r.h5mu"
    with mudata.set_options(pull_on_update=False):
        md.write_h5mu(p)
        mb = mudata.read_h5mu(p)
    pd.testing.assert_frame_equal(
        md.mod["peptidoform"].varm["fasta_validation"],
        mb.mod["peptidoform"].varm["fasta_validation"],
    )


def test_revalidation_overwrites():
    a = _peptide_adata(["DTHK"])
    validate_peptides_against_fasta(a, FASTA)
    # a second run must not raise and must replace the prior summary
    validate_peptides_against_fasta(a, FASTA)
    assert a.varm["fasta_validation"].shape[0] == 1


# --- provenance -----------------------------------------------------------


def _last_provenance(adata):
    import json

    entries = json.loads(adata.uns["anndata_proteomics"]["var_annotations_json"])
    return entries[-1]


def test_provenance_entry_content():
    a = _peptide_adata(["DTHK", "ZZZZZ"])
    validate_peptides_against_fasta(a, FASTA, backend="ahocorapy")
    entry = _last_provenance(a)
    assert entry["source"] == "fasta_validation"
    assert entry["destination"] == "varm['fasta_validation']"
    assert entry["requested_backend"] == "ahocorapy"
    assert entry["backend"] == "ahocorapy"
    assert entry["sequence_field"] == "ProForma_peptide"
    assert entry["fasta_config"]["decoy"]["patterns"] == ["^REV_"]
    assert entry["il_equivalent"] is False
    assert entry["n_features"] == 2
    assert entry["n_matched_features"] == 1
    assert entry["n_unmatched_features"] == 1


def test_auto_backend_provenance_records_requested_and_resolved_backend():
    a = _peptide_adata(["DTHK"])
    result = validate_peptides_against_fasta(a, FASTA, backend="auto")

    entry = _last_provenance(a)
    assert result.requested_backend == "auto"
    assert entry["requested_backend"] == "auto"
    assert result.backend in get_available_backends()
    assert entry["backend"] == result.backend


def test_provenance_appends_without_clobbering():
    import json

    a = _peptide_adata(["DTHK"])
    # a pre-existing entry (e.g. from an obs/var annotation step) must survive
    a.uns["anndata_proteomics"]["var_annotations_json"] = json.dumps([{"source": "preexisting"}])
    validate_peptides_against_fasta(a, FASTA)
    entries = json.loads(a.uns["anndata_proteomics"]["var_annotations_json"])
    assert [e["source"] for e in entries] == ["preexisting", "fasta_validation"]
    # a second validation appends another entry
    validate_peptides_against_fasta(a, FASTA)
    entries = json.loads(a.uns["anndata_proteomics"]["var_annotations_json"])
    assert [e["source"] for e in entries] == [
        "preexisting",
        "fasta_validation",
        "fasta_validation",
    ]


def test_provenance_round_trips(tmp_path: Path):
    a = _peptide_adata(["DTHK"])
    validate_peptides_against_fasta(a, FASTA)
    p = tmp_path / "prov.h5ad"
    a.write_h5ad(p)
    b = ad.read_h5ad(p)
    assert _last_provenance(b)["source"] == "fasta_validation"


# --- CLI ------------------------------------------------------------------


def test_cli_fasta_validates_peptide_only_input_by_default(tmp_path: Path):
    from anndata_proteomics.scripts.cli import fasta as cli_fasta

    data = tmp_path / "in.h5ad"
    _peptide_adata(["DTHK", "ZZZZZ"]).write_h5ad(data)
    fasta = tmp_path / "db.fasta"
    fasta.write_text(">sp|A|ONE\nMKDTHKXX\n")

    rc = cli_fasta(data, fasta)
    assert rc == 0
    out = data.with_name("in.annotated.h5ad")
    assert out.exists()
    assert "fasta_validation" in ad.read_h5ad(out).varm

    result = ad.read_h5ad(out)
    assert result.n_vars == 2
    assert result.varm["fasta_validation"]["peptide_in_fasta"].tolist() == [True, False]
    from anndata_proteomics.readers.summary import describe

    assert describe(result)["fasta"]["proteotypic_feature_count"] == 1


def test_generator_sources_keep_provenance() -> None:
    a = _peptide_adata(["DTHK"])
    validate_peptides_against_fasta(a, (source for source in [FASTA]))
    assert _last_provenance(a)["fasta_sources"] == ["<inline-fasta>"]
