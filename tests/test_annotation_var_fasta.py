"""Tests for FASTA → protein-var annotation (anndata_proteomics.annotation.var_fasta).

Synthetic AnnData/MuData keep these independent of the cached test-data catalog.
Warnings go through loguru → stderr; the `_loguru_to_pytest_capsys` fixture in
conftest.py wires that into pytest capture, so we read `capsys.readouterr().err`.
"""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import mudata
import numpy as np
import pandas as pd
import pytest
from mudata import MuData

from anndata_proteomics.annotation.var_fasta import annotate_var_from_fasta
from anndata_proteomics.fasta.annotation import count_peptides
from anndata_proteomics.params.anndata_io import write_search_parameters
from anndata_proteomics.params.model import Parameters
from anndata_proteomics.scripts.cli import fasta as fasta_cmd

# A few forward UniProt records (lifted from prolfquapp's fixture) + one contaminant
# + one REV_ decoy. ≥2 carry GN=, so the gene_name column is emitted.
FASTA = """\
>sp|A0A385XJL2|YGDT_ECOLI Protein YgdT OS=Escherichia coli OX=83333 GN=ygdT PE=4 SV=1
MLSTESWDNCEKPPLLFPFTALTCDETPVFSGSVLNLVAHSVDKYGIG
>sp|P03018|UVRD_ECOLI DNA helicase II OS=Escherichia coli OX=83333 GN=uvrD PE=1 SV=1
MDVSYLLDSLNDKQREAVAAPRSNLLVLAGAGSGKTRVLVHRIAWLMSVENCSPYSIMAV
>sp|P04982|RBSD_ECOLI D-ribose pyranase OS=Escherichia coli OX=83333 GN=rbsD PE=1 SV=3
MKKGTVLNSDISSVISRLGHTDTLVVCDAGLPIPKSTTRIDMALTQGVPSFMQVLGVVTN
>sp|P04994|EX7L_ECOLI Exodeoxyribonuclease 7 OS=Escherichia coli OX=83333 GN=xseA PE=1 SV=2
MLPSQSPAIFTVSRLNQTVRLLLEHEMGQVWISGEISNFTQPASGHWYFTLKDDTAQVRC
>REV_sp|Q13515|BFSP2_HUMAN Phakinin OS=Homo sapiens OX=9606 GN=BFSP2 PE=1 SV=1
GSEERDLLAHYSAVDKQLQCKRALLHAREQQQQEAEARIERLEAELRGVVAGLNQLEMDH
"""

SEQ_P03018 = "MDVSYLLDSLNDKQREAVAAPRSNLLVLAGAGSGKTRVLVHRIAWLMSVENCSPYSIMAV"


def _protein_adata(
    var_names: list[str],
    *,
    level: str = "protein",
    with_group_column: bool = True,
) -> ad.AnnData:
    n = len(var_names)
    var = pd.DataFrame(index=pd.Index(var_names))
    if with_group_column:
        var["Protein_Group"] = list(var_names)
    adata = ad.AnnData(
        X=np.arange(2 * n, dtype="float64").reshape(2, n),
        obs=pd.DataFrame(index=pd.Index(["run1", "run2"])),
        var=var,
    )
    adata.uns["anndata_proteomics"] = {"quantification_level": level}
    return adata


def _ion_adata(var_names: list[str]) -> ad.AnnData:
    n = len(var_names)
    return ad.AnnData(
        X=np.arange(2 * n, dtype="float64").reshape(2, n),
        obs=pd.DataFrame(index=pd.Index(["run1", "run2"])),
        var=pd.DataFrame(index=pd.Index(var_names)),
    )


# --- happy paths -------------------------------------------------------------


def test_varm_fasta_has_expected_columns() -> None:
    adata = _protein_adata(["P03018", "A0A385XJL2"])
    annotate_var_from_fasta(adata, FASTA, cleavage="Trypsin")
    fa = adata.varm["fasta"]
    assert {
        "fasta.id",
        "fasta.header",
        "protein_length",
        "nr_peptides",
        "gene_name",
    } <= set(fa.columns)
    assert list(fa.index) == list(adata.var_names)  # var-aligned
    assert fa.loc["P03018", "fasta.id"] == "sp|P03018|UVRD_ECOLI"
    assert fa.loc["P03018", "gene_name"] == "uvrD"
    assert fa.loc["P03018", "protein_length"] == len(SEQ_P03018)


def test_leading_accession_join_splits_group_and_uniprot_form() -> None:
    # "P04982;Q99999" -> first token P04982; "sp|P04994|..." -> middle P04994.
    adata = _protein_adata(["P04982;Q99999", "sp|P04994|EX7L_ECOLI"])
    annotate_var_from_fasta(adata, FASTA, cleavage="Trypsin")
    fa = adata.varm["fasta"]
    assert fa.loc["P04982;Q99999", "gene_name"] == "rbsD"
    assert fa.loc["sp|P04994|EX7L_ECOLI", "gene_name"] == "xseA"


def test_match_on_index_strips_prt_prefix() -> None:
    adata = _protein_adata(["prt:P03018", "prt:A0A385XJL2"], with_group_column=False)
    annotate_var_from_fasta(adata, FASTA, match_on="index", cleavage="Trypsin")
    assert adata.varm["fasta"].loc["prt:P03018", "fasta.id"] == "sp|P03018|UVRD_ECOLI"


def test_columns_subset_restricts_stored_columns() -> None:
    adata = _protein_adata(["P03018"])
    annotate_var_from_fasta(adata, FASTA, cleavage="Trypsin", columns=["nr_peptides"])
    assert list(adata.varm["fasta"].columns) == ["nr_peptides"]


# --- cleavage / enzyme -------------------------------------------------------


def test_enzyme_read_from_search_parameters_drives_count() -> None:
    lysc = count_peptides(SEQ_P03018, cleavage="Lys-C", min_length=7, max_length=30)
    trypsin = count_peptides(SEQ_P03018, cleavage="Trypsin", min_length=7, max_length=30)
    assert lysc != trypsin  # precondition: the enzyme must matter for this protein

    adata = _protein_adata(["P03018"])
    write_search_parameters(
        adata,
        Parameters(enzyme="Lys-C", min_peptide_length=7, max_peptide_length=30),
    )
    annotate_var_from_fasta(adata, FASTA)  # no cleavage arg => read from params
    assert adata.varm["fasta"].loc["P03018", "nr_peptides"] == lysc


def test_cleavage_override_wins_over_params() -> None:
    adata = _protein_adata(["P03018"])
    write_search_parameters(adata, Parameters(enzyme="Lys-C"))
    annotate_var_from_fasta(adata, FASTA, cleavage="Trypsin/P", min_length=7, max_length=30)
    expected = count_peptides(SEQ_P03018, cleavage="Trypsin/P", min_length=7, max_length=30)
    assert adata.varm["fasta"].loc["P03018", "nr_peptides"] == expected


def test_no_params_warns_and_defaults_to_trypsin(capsys: pytest.CaptureFixture[str]) -> None:
    adata = _protein_adata(["P03018"])  # no search parameters stored
    annotate_var_from_fasta(adata, FASTA)
    err = capsys.readouterr().err
    assert "no enzyme in search parameters" in err
    expected = count_peptides(SEQ_P03018, cleavage="Trypsin", min_length=7, max_length=30)
    assert adata.varm["fasta"].loc["P03018", "nr_peptides"] == expected


def test_unknown_enzyme_override_warns_and_falls_back(
    capsys: pytest.CaptureFixture[str],
) -> None:
    adata = _protein_adata(["P03018"])
    annotate_var_from_fasta(adata, FASTA, cleavage="Pepsin", min_length=7, max_length=30)
    assert "unknown enzyme 'Pepsin'" in capsys.readouterr().err
    expected = count_peptides(SEQ_P03018, cleavage="Trypsin", min_length=7, max_length=30)
    assert adata.varm["fasta"].loc["P03018", "nr_peptides"] == expected


# --- MuData ------------------------------------------------------------------


def test_mudata_annotates_protein_modality_only() -> None:
    prot = _protein_adata(["prt:P03018", "prt:A0A385XJL2"], with_group_column=False)
    ion = _ion_adata(["ion:a", "ion:b"])
    with mudata.set_options(pull_on_update=False):
        md = MuData({"ion": ion, "protein": prot}, axis=0)
    annotate_var_from_fasta(md, FASTA, match_on="index", cleavage="Trypsin")
    assert "fasta" in md.mod["protein"].varm
    assert "fasta" not in md.mod["ion"].varm


def test_mudata_roundtrips_through_h5mu(tmp_path: Path) -> None:
    prot = _protein_adata(["prt:P03018", "prt:A0A385XJL2"], with_group_column=False)
    ion = _ion_adata(["ion:a", "ion:b"])
    with mudata.set_options(pull_on_update=False):
        md = MuData({"ion": ion, "protein": prot}, axis=0)
    annotate_var_from_fasta(md, FASTA, match_on="index", cleavage="Trypsin")

    out = tmp_path / "md.annotated.h5mu"
    md.write_h5mu(out)
    with mudata.set_options(pull_on_update=False):
        rt = mudata.read_h5mu(out)
    fa = rt.mod["protein"].varm["fasta"]
    assert fa.loc["prt:P03018", "fasta.id"] == "sp|P03018|UVRD_ECOLI"


def test_mudata_without_protein_modality_raises() -> None:
    ion = _ion_adata(["ion:a", "ion:b"])
    with mudata.set_options(pull_on_update=False):
        md = MuData({"ion": ion}, axis=0)
    with pytest.raises(ValueError, match="no 'protein' modality"):
        annotate_var_from_fasta(md, FASTA, cleavage="Trypsin")


# --- guards / mismatch -------------------------------------------------------


def test_non_protein_anndata_raises() -> None:
    adata = _protein_adata(["P03018"], level="ion")
    with pytest.raises(ValueError, match="protein layer only"):
        annotate_var_from_fasta(adata, FASTA, cleavage="Trypsin")


def test_zero_match_raises() -> None:
    adata = _protein_adata(["NOSUCH1", "NOSUCH2"])
    with pytest.raises(ValueError, match="no var rows matched"):
        annotate_var_from_fasta(adata, FASTA, cleavage="Trypsin")


def test_partial_match_warns_and_leaves_nan(capsys: pytest.CaptureFixture[str]) -> None:
    adata = _protein_adata(["P03018", "NOSUCH"])
    annotate_var_from_fasta(adata, FASTA, cleavage="Trypsin")
    assert "1/2 var rows had no matching" in capsys.readouterr().err
    fa = adata.varm["fasta"]
    assert fa.loc["P03018", "gene_name"] == "uvrD"
    assert pd.isna(fa.loc["NOSUCH", "gene_name"])


def test_rerun_varm_collision_raises() -> None:
    adata = _protein_adata(["P03018"])
    annotate_var_from_fasta(adata, FASTA, cleavage="Trypsin")
    with pytest.raises(ValueError, match="already present"):
        annotate_var_from_fasta(adata, FASTA, cleavage="Trypsin")


def test_unknown_match_on_column_raises() -> None:
    adata = _protein_adata(["P03018"], with_group_column=False)
    with pytest.raises(ValueError, match="match_on column 'Protein_Group' not found"):
        annotate_var_from_fasta(adata, FASTA, cleavage="Trypsin")


# --- provenance --------------------------------------------------------------


def test_provenance_records_enzyme_and_sources() -> None:
    adata = _protein_adata(["P03018"])
    write_search_parameters(adata, Parameters(enzyme="Lys-C"))
    annotate_var_from_fasta(adata, FASTA)
    entries = json.loads(adata.uns["anndata_proteomics"]["var_annotations_json"])
    assert len(entries) == 1
    entry = entries[0]
    assert entry["source"] == "fasta"
    assert entry["destination"] == "varm['fasta']"
    assert entry["cleavage_enzyme"] == "Lys-C"
    assert entry["columns"]
    assert entry["fasta_sources"] == ["<inline-fasta>"]


# --- CLI ---------------------------------------------------------------------


def test_cli_fasta_writes_annotated_file(tmp_path: Path) -> None:
    adata = _protein_adata(["P03018", "A0A385XJL2"])
    data_path = tmp_path / "proteins.h5ad"
    adata.write_h5ad(data_path)

    fasta_path = tmp_path / "db.fasta"
    fasta_path.write_text(FASTA)

    out = tmp_path / "proteins.annotated.h5ad"
    rc = fasta_cmd(data_path, fasta_path, output=out, cleavage="Trypsin")
    assert rc == 0
    assert out.exists()

    rt = ad.read_h5ad(out)
    assert rt.varm["fasta"].loc["P03018", "fasta.id"] == "sp|P03018|UVRD_ECOLI"


def test_cli_fasta_requires_a_fasta_file(tmp_path: Path) -> None:
    adata = _protein_adata(["P03018"])
    data_path = tmp_path / "proteins.h5ad"
    adata.write_h5ad(data_path)
    assert fasta_cmd(data_path) == 1
