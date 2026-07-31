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

from anndata_proteomics.adapters.anndata.fasta import (
    DEFAULT_ANNDATA_PROTEIN_FASTA_CONFIG,
    FASTA_ACCESSIONS_INDEX,
    AnnDataProteinFastaConfig,
    FastaAccessionsColumn,
    read_protein_annotation_input,
    require_protein_mudata_target,
    write_anndata_fasta_config,
    write_mudata_fasta_config,
    write_protein_annotation,
)
from anndata_proteomics.adapters.anndata.params import write_search_parameters
from anndata_proteomics.annotation.var_fasta import (
    ProteinFastaAnnotationResult,
    SelectedFastaColumns,
    annotate_proteins_from_fasta,
)
from anndata_proteomics.fasta.annotation import count_peptides, resolve_cleavage_name
from anndata_proteomics.fasta.parser import FastaSources
from anndata_proteomics.params.model import Parameters
from anndata_proteomics.rules.loader import load_rule
from anndata_proteomics.rules.registry import find_rule_for_version
from anndata_proteomics.rules.schema import ParseRule, QuantificationLevel
from anndata_proteomics.scripts.cli import FastaCliOptions
from anndata_proteomics.scripts.cli import fasta as fasta_cmd
from anndata_proteomics.workflows.fasta import (
    MaximumPeptideLength,
    MinimumPeptideLength,
    NamedCleavage,
    resolve_protein_annotation_input,
)

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
TRYPSIN_CONFIG = AnnDataProteinFastaConfig(cleavage=NamedCleavage("Trypsin"))


def annotate_var_from_fasta(
    target: ad.AnnData,
    fasta_sources: FastaSources,
    config: AnnDataProteinFastaConfig = DEFAULT_ANNDATA_PROTEIN_FASTA_CONFIG,
) -> ProteinFastaAnnotationResult:
    """Compose protein annotation and AnnData persistence for tests."""
    extracted = resolve_protein_annotation_input(read_protein_annotation_input(target, config))
    result = annotate_proteins_from_fasta(
        extracted.protein_groups,
        fasta_sources,
        extracted.config,
    )
    write_protein_annotation(target, result, extracted.provenance)
    write_anndata_fasta_config(target, result.fasta_config)
    return result


def annotate_protein_mudata_from_fasta(
    target: MuData,
    fasta_sources: FastaSources,
    config: AnnDataProteinFastaConfig = DEFAULT_ANNDATA_PROTEIN_FASTA_CONFIG,
) -> ProteinFastaAnnotationResult:
    """Compose protein annotation and MuData persistence for tests."""
    protein = require_protein_mudata_target(target)
    extracted = resolve_protein_annotation_input(read_protein_annotation_input(protein, config))
    result = annotate_proteins_from_fasta(
        extracted.protein_groups,
        fasta_sources,
        extracted.config,
    )
    write_protein_annotation(protein, result, extracted.provenance)
    write_mudata_fasta_config(target, result.fasta_config)
    return result


def _protein_adata(
    var_names: list[str],
    *,
    level: QuantificationLevel = "protein",
    with_group_column: bool = True,
    with_rule: bool = True,
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
    if with_group_column and with_rule:
        _attach_test_role_rule(adata, level, "Protein_Group")
    return adata


def _attach_test_role_rule(
    adata: ad.AnnData,
    level: QuantificationLevel,
    column: str,
) -> None:
    """Attach a minimal effective rule declaring both protein semantics."""
    rule = ParseRule.model_validate(
        {
            "schema_version": "0.1",
            "file_version": "test",
            "software_name": "test",
            "software_version": ".*",
            "input_shape": "long",
            "quantification_level": level,
            "axis": {
                "obs_keys": ["Run"],
                "var_keys": [column],
                "x_layer": "Intensity",
            },
            "columns": {
                "obs": {"select": {"Run": "Run"}},
                "var": {"select": {column: column}},
            },
            "column_roles": {
                "protein_assignment": column,
                "fasta_accessions": column,
            },
            "layers": [{"name": "Intensity", "source": "Intensity"}],
        }
    )
    adata.uns["anndata_proteomics"]["rule_json"] = rule.model_dump_json(by_alias=True)


def _attach_diann_rule(adata: ad.AnnData, level: QuantificationLevel) -> None:
    """Attach a packaged DIA-NN rule as conversion provenance."""
    rule = load_rule(find_rule_for_version("diann", level, "1.9.2"))
    adata.uns["anndata_proteomics"]["rule_json"] = rule.model_dump_json(by_alias=True)


def _ion_adata(var_names: list[str]) -> ad.AnnData:
    n = len(var_names)
    return ad.AnnData(
        X=np.arange(2 * n, dtype="float64").reshape(2, n),
        obs=pd.DataFrame(index=pd.Index(["run1", "run2"])),
        var=pd.DataFrame(index=pd.Index(var_names)),
    )


def _varm_frame(adata: ad.AnnData | MuData, key: str = "fasta") -> pd.DataFrame:
    """Return a DataFrame-valued varm entry after checking its runtime type."""
    value = adata.varm[key]
    if not isinstance(value, pd.DataFrame):
        raise AssertionError(f"expected varm[{key!r}] to contain a DataFrame")
    return value


# --- happy paths -------------------------------------------------------------


def test_varm_fasta_has_expected_columns() -> None:
    adata = _protein_adata(["P03018", "A0A385XJL2"])
    annotate_var_from_fasta(
        adata,
        FASTA,
        TRYPSIN_CONFIG,
    )
    fa = _varm_frame(adata)
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
    annotate_var_from_fasta(
        adata,
        FASTA,
        TRYPSIN_CONFIG,
    )
    fa = _varm_frame(adata)
    assert fa.loc["P04982;Q99999", "gene_name"] == "rbsD"
    assert fa.loc["sp|P04994|EX7L_ECOLI", "gene_name"] == "xseA"


def test_declared_protein_role_wins_over_plausible_columns() -> None:
    adata = _protein_adata(["feature"], with_group_column=False)
    adata.var["Protein_Group"] = ["NOT_IN_FASTA"]
    adata.var["Protein_Ids"] = ["P03018"]
    _attach_diann_rule(adata, "protein")

    annotate_var_from_fasta(
        adata,
        FASTA,
        TRYPSIN_CONFIG,
    )

    assert _varm_frame(adata).loc["feature", "fasta.id"] == "sp|P03018|UVRD_ECOLI"


def test_match_column_is_not_guessed_without_conversion_rule() -> None:
    adata = _protein_adata(["feature"], with_rule=False)

    with pytest.raises(ValueError, match="does not declare column_roles.fasta_accessions"):
        annotate_var_from_fasta(adata, FASTA, TRYPSIN_CONFIG)


def test_match_column_does_not_cascade_from_assignment_role() -> None:
    adata = _protein_adata(["feature"], level="protein", with_group_column=False)
    adata.var["Genes"] = ["PLSL_HUMAN"]
    rule = load_rule(find_rule_for_version("alphadia", "ion", "1.10.3"))
    adata.uns["anndata_proteomics"]["rule_json"] = rule.model_dump_json(by_alias=True)

    with pytest.raises(ValueError, match="does not declare column_roles.fasta_accessions"):
        annotate_var_from_fasta(adata, FASTA, TRYPSIN_CONFIG)


def test_match_on_index_strips_prt_prefix() -> None:
    adata = _protein_adata(["prt:P03018", "prt:A0A385XJL2"], with_group_column=False)
    annotate_var_from_fasta(
        adata,
        FASTA,
        AnnDataProteinFastaConfig(
            match_on=FASTA_ACCESSIONS_INDEX,
            cleavage=NamedCleavage("Trypsin"),
        ),
    )
    assert _varm_frame(adata).loc["prt:P03018", "fasta.id"] == "sp|P03018|UVRD_ECOLI"


def test_columns_subset_restricts_stored_columns() -> None:
    adata = _protein_adata(["P03018"])
    annotate_var_from_fasta(
        adata,
        FASTA,
        AnnDataProteinFastaConfig(
            cleavage=NamedCleavage("Trypsin"),
            columns=SelectedFastaColumns(("nr_peptides",)),
        ),
    )
    assert list(_varm_frame(adata).columns) == ["nr_peptides"]


# --- cleavage / enzyme -------------------------------------------------------


def test_enzyme_read_from_search_parameters_drives_count() -> None:
    lysc = count_peptides(
        SEQ_P03018,
        cleavage=resolve_cleavage_name("Lys-C").rule,
        min_length=7,
        max_length=30,
    )
    trypsin = count_peptides(
        SEQ_P03018,
        cleavage=resolve_cleavage_name("Trypsin").rule,
        min_length=7,
        max_length=30,
    )
    assert lysc != trypsin  # precondition: the enzyme must matter for this protein

    adata = _protein_adata(["P03018"])
    write_search_parameters(
        adata,
        Parameters(enzyme="Lys-C", min_peptide_length=7, max_peptide_length=30),
    )
    annotate_var_from_fasta(adata, FASTA)  # no cleavage arg => read from params
    assert _varm_frame(adata).loc["P03018", "nr_peptides"] == lysc


def test_cleavage_override_wins_over_params() -> None:
    adata = _protein_adata(["P03018"])
    write_search_parameters(adata, Parameters(enzyme="Lys-C"))
    annotate_var_from_fasta(
        adata,
        FASTA,
        AnnDataProteinFastaConfig(
            cleavage=NamedCleavage("Trypsin/P"),
            minimum_length=MinimumPeptideLength(7),
            maximum_length=MaximumPeptideLength(30),
        ),
    )
    expected = count_peptides(
        SEQ_P03018,
        cleavage=resolve_cleavage_name("Trypsin/P").rule,
        min_length=7,
        max_length=30,
    )
    assert _varm_frame(adata).loc["P03018", "nr_peptides"] == expected


def test_no_params_warns_and_defaults_to_trypsin(capsys: pytest.CaptureFixture[str]) -> None:
    adata = _protein_adata(["P03018"])  # no search parameters stored
    annotate_var_from_fasta(adata, FASTA)
    err = capsys.readouterr().err
    assert "no enzyme in search parameters" in err
    expected = count_peptides(
        SEQ_P03018,
        cleavage=resolve_cleavage_name("Trypsin").rule,
        min_length=7,
        max_length=30,
    )
    assert _varm_frame(adata).loc["P03018", "nr_peptides"] == expected


def test_unknown_enzyme_override_warns_and_falls_back(
    capsys: pytest.CaptureFixture[str],
) -> None:
    adata = _protein_adata(["P03018"])
    annotate_var_from_fasta(
        adata,
        FASTA,
        AnnDataProteinFastaConfig(
            cleavage=NamedCleavage("Pepsin"),
            minimum_length=MinimumPeptideLength(7),
            maximum_length=MaximumPeptideLength(30),
        ),
    )
    assert "unknown enzyme 'Pepsin'" in capsys.readouterr().err
    expected = count_peptides(
        SEQ_P03018,
        cleavage=resolve_cleavage_name("Trypsin").rule,
        min_length=7,
        max_length=30,
    )
    assert _varm_frame(adata).loc["P03018", "nr_peptides"] == expected


# --- MuData ------------------------------------------------------------------


def test_mudata_annotates_protein_modality_only() -> None:
    prot = _protein_adata(["prt:P03018", "prt:A0A385XJL2"], with_group_column=False)
    ion = _ion_adata(["ion:a", "ion:b"])
    with mudata.set_options(pull_on_update=False):
        md = MuData({"ion": ion, "protein": prot}, axis=0)
    annotate_protein_mudata_from_fasta(
        md,
        FASTA,
        AnnDataProteinFastaConfig(
            match_on=FASTA_ACCESSIONS_INDEX,
            cleavage=NamedCleavage("Trypsin"),
        ),
    )
    assert "fasta" in md.mod["protein"].varm
    assert "fasta" not in md.mod["ion"].varm


def test_mudata_roundtrips_through_h5mu(tmp_path: Path) -> None:
    prot = _protein_adata(["prt:P03018", "prt:A0A385XJL2"], with_group_column=False)
    ion = _ion_adata(["ion:a", "ion:b"])
    with mudata.set_options(pull_on_update=False):
        md = MuData({"ion": ion, "protein": prot}, axis=0)
    annotate_protein_mudata_from_fasta(
        md,
        FASTA,
        AnnDataProteinFastaConfig(
            match_on=FASTA_ACCESSIONS_INDEX,
            cleavage=NamedCleavage("Trypsin"),
        ),
    )

    out = tmp_path / "md.annotated.h5mu"
    md.write_h5mu(out)
    with mudata.set_options(pull_on_update=False):
        rt = mudata.read_h5mu(out)
    fa = _varm_frame(rt.mod["protein"])
    assert fa.loc["prt:P03018", "fasta.id"] == "sp|P03018|UVRD_ECOLI"


def test_mudata_without_protein_modality_raises() -> None:
    ion = _ion_adata(["ion:a", "ion:b"])
    with mudata.set_options(pull_on_update=False):
        md = MuData({"ion": ion}, axis=0)
    with pytest.raises(ValueError, match="no 'protein' modality"):
        annotate_protein_mudata_from_fasta(
            md,
            FASTA,
            TRYPSIN_CONFIG,
        )


# --- guards / mismatch -------------------------------------------------------


def test_non_protein_anndata_raises() -> None:
    adata = _protein_adata(["P03018"], level="ion")
    with pytest.raises(ValueError, match="protein layer only"):
        annotate_var_from_fasta(
            adata,
            FASTA,
            TRYPSIN_CONFIG,
        )


def test_zero_match_raises() -> None:
    adata = _protein_adata(["NOSUCH1", "NOSUCH2"])
    with pytest.raises(ValueError, match="no var rows matched"):
        annotate_var_from_fasta(
            adata,
            FASTA,
            TRYPSIN_CONFIG,
        )


def test_partial_match_warns_and_roundtrips_nullable_flags(
    capsys: pytest.CaptureFixture[str],
    tmp_path: Path,
) -> None:
    adata = _protein_adata(["P03018", "NOSUCH"])
    annotate_var_from_fasta(
        adata,
        FASTA,
        TRYPSIN_CONFIG,
    )
    assert "1/2 var rows had no matching" in capsys.readouterr().err
    fa = _varm_frame(adata)
    assert fa.loc["P03018", "gene_name"] == "uvrD"
    assert pd.isna(fa.loc["NOSUCH", "gene_name"])
    assert str(fa["is_decoy"].dtype) == "boolean"
    assert isinstance(fa["fasta.id"].dtype, pd.CategoricalDtype)
    assert pd.isna(fa.loc["NOSUCH", "is_decoy"])

    path = tmp_path / "partial.h5ad"
    adata.write_h5ad(path)
    restored = _varm_frame(ad.read_h5ad(path))
    assert not restored.loc["P03018", "is_decoy"]
    assert pd.isna(restored.loc["NOSUCH", "is_decoy"])


def test_rerun_varm_collision_raises() -> None:
    adata = _protein_adata(["P03018"])
    annotate_var_from_fasta(
        adata,
        FASTA,
        TRYPSIN_CONFIG,
    )
    with pytest.raises(ValueError, match="already present"):
        annotate_var_from_fasta(
            adata,
            FASTA,
            TRYPSIN_CONFIG,
        )


def test_unknown_match_on_column_raises() -> None:
    adata = _protein_adata(["P03018"], with_group_column=False)
    with pytest.raises(ValueError, match="match_on 'Protein_Group' not in var columns"):
        annotate_var_from_fasta(
            adata,
            FASTA,
            AnnDataProteinFastaConfig(
                match_on=FastaAccessionsColumn("Protein_Group"),
                cleavage=NamedCleavage("Trypsin"),
            ),
        )


def test_decoy_quantification_is_retained_and_annotated() -> None:
    adata = _protein_adata(["REV_Q13515"])
    matrix = adata.X
    assert isinstance(matrix, np.ndarray)
    before = matrix.copy()
    annotate_var_from_fasta(
        adata,
        FASTA,
        TRYPSIN_CONFIG,
    )
    assert adata.n_vars == 1
    np.testing.assert_array_equal(adata.X, before)
    fasta_frame = _varm_frame(adata)
    assert bool(fasta_frame.loc["REV_Q13515", "is_decoy"])
    assert fasta_frame.loc["REV_Q13515", "fasta.id"] == "REV_sp|Q13515|BFSP2_HUMAN"


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
    rc = fasta_cmd(
        data_path,
        fasta_path,
        options=FastaCliOptions(output=out, cleavage="Trypsin"),
    )
    assert rc == 0
    assert out.exists()

    rt = ad.read_h5ad(out)
    assert _varm_frame(rt).loc["P03018", "fasta.id"] == "sp|P03018|UVRD_ECOLI"
    annotation = json.loads(rt.uns["anndata_proteomics"]["var_annotations_json"])[-1]
    assert annotation["source"] == "fasta"
    assert annotation["destination"] == "varm['fasta']"
    assert annotation["n_var_matched"] == 2


def test_cli_fasta_requires_a_fasta_file(tmp_path: Path) -> None:
    adata = _protein_adata(["P03018"])
    data_path = tmp_path / "proteins.h5ad"
    adata.write_h5ad(data_path)
    assert fasta_cmd(data_path) == 1


def test_cli_fasta_validates_all_mudata_layers_by_default(tmp_path: Path) -> None:
    protein = _protein_adata(["prt:P03018"], with_group_column=False)
    _attach_diann_rule(protein, "protein")
    ion = ad.AnnData(
        X=np.ones((2, 1)),
        obs=pd.DataFrame(index=["run1", "run2"]),
        var=pd.DataFrame(
            {
                "ProForma_peptide": ["MDVSY"],
                "Protein_Group": ["P03018"],
            },
            index=["ion:MDVSY"],
        ),
    )
    ion.uns["anndata_proteomics"] = {"quantification_level": "ion"}
    _attach_diann_rule(ion, "ion")
    with mudata.set_options(pull_on_update=False):
        mdata = MuData({"ion": ion, "protein": protein}, axis=0)
    data_path = tmp_path / "levels.h5mu"
    output_path = tmp_path / "levels.fasta.h5mu"
    fasta_path = tmp_path / "db.fasta"
    mdata.write_h5mu(data_path)
    fasta_path.write_text(FASTA)

    assert (
        fasta_cmd(
            data_path,
            fasta_path,
            options=FastaCliOptions(
                output=output_path,
                cleavage="Trypsin",
                match_on="index",
                leading_protein_field="Protein_Group",
            ),
        )
        == 0
    )

    with mudata.set_options(pull_on_update=False):
        restored = mudata.read_h5mu(output_path)
    assert "fasta" in restored.mod["protein"].varm
    assert _varm_frame(restored.mod["ion"], "fasta_validation").loc[
        "ion:MDVSY",
        "peptide_in_leading_protein",
    ]
    assert restored.varp["feature_mapping"].nnz == 1


def test_cli_fasta_no_validate_only_annotates_proteins(tmp_path: Path) -> None:
    protein = _protein_adata(["prt:P03018"], with_group_column=False)
    _attach_diann_rule(protein, "protein")
    ion = ad.AnnData(
        X=np.ones((2, 1)),
        obs=pd.DataFrame(index=["run1", "run2"]),
        var=pd.DataFrame(
            {"ProForma_peptide": ["MDVSY"]},
            index=["ion:MDVSY"],
        ),
    )
    ion.uns["anndata_proteomics"] = {"quantification_level": "ion"}
    _attach_diann_rule(ion, "ion")
    with mudata.set_options(pull_on_update=False):
        mdata = MuData({"ion": ion, "protein": protein}, axis=0)
    data_path = tmp_path / "levels.h5mu"
    output_path = tmp_path / "levels.no-validation.h5mu"
    fasta_path = tmp_path / "db.fasta"
    mdata.write_h5mu(data_path)
    fasta_path.write_text(FASTA)

    assert (
        fasta_cmd(
            data_path,
            fasta_path,
            options=FastaCliOptions(
                output=output_path,
                cleavage="Trypsin",
                match_on="index",
                validate=False,
            ),
        )
        == 0
    )

    with mudata.set_options(pull_on_update=False):
        restored = mudata.read_h5mu(output_path)
    assert "fasta" in restored.mod["protein"].varm
    assert "fasta_validation" not in restored.mod["ion"].varm
    assert "feature_mapping" not in restored.varp
