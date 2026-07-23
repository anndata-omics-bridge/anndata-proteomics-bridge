"""ProteoBench matrix scoring and storage contracts."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import mudata
import numpy as np
import pandas as pd
import pytest
from mudata import MuData
from scipy import sparse

from anndata_proteomics.annotation.apply import annotate_obs
from anndata_proteomics.annotation.loader import AnnotationTable
from anndata_proteomics.annotation.validate_fasta import validate_peptides_against_fasta
from anndata_proteomics.proteobench.config import (
    ExpectedRatio,
    ModuleGeneral,
    ModuleSettings,
    ModificationParserSettings,
    SampleSettings,
    ToolGeneral,
    ToolSettings,
    load_module_settings,
    load_tool_settings,
)
from anndata_proteomics.proteobench.intermediate import align_runs, compute_intermediate
from anndata_proteomics.proteobench.mapping import (
    parse_proteobench_features,
    render_proteobench_features,
)
from anndata_proteomics.proteobench.metrics import build_scores, compute_roc_auc
from anndata_proteomics.proteobench.pipeline import score_quantification
from anndata_proteomics.proteobench.resolve import resolve_roles
from anndata_proteomics.readers.summary import describe
from anndata_proteomics.rules.schema import ParseRule
from anndata_proteomics.scripts.cli import proteobench as proteobench_cmd


def _module_settings() -> ModuleSettings:
    return ModuleSettings(
        species_expected_ratio={
            "YEAST": ExpectedRatio(A_vs_B=2.0),
            "ECOLI": ExpectedRatio(A_vs_B=0.25),
            "HUMAN": ExpectedRatio(A_vs_B=1.0),
        },
        species_mapper={"_YEAST": "YEAST", "_ECOLI": "ECOLI", "_HUMAN": "HUMAN"},
        general=ModuleGeneral(min_count_multispec=1, level="ion"),
        samples=[
            SampleSettings(raw_file="run_A1", sample_name="A1", condition="A"),
            SampleSettings(raw_file="run_A2", sample_name="A2", condition="A"),
            SampleSettings(raw_file="run_B1", sample_name="B1", condition="B"),
            SampleSettings(raw_file="run_B2", sample_name="B2", condition="B"),
        ],
    )


def _tool_settings() -> ToolSettings:
    return ToolSettings(
        mapper={
            "Run": "Raw file",
            "Protein.Ids": "Proteins",
            "Precursor.Normalised": "Intensity",
        },
        general=ToolGeneral(contaminant_flag="Cont_", decoy_flag=True),
    )


def _rule() -> ParseRule:
    return ParseRule.model_validate(
        {
            "schema_version": "0.1",
            "file_version": "test",
            "software_name": "DIA-NN",
            "software_version": ".*",
            "input_shape": "long",
            "quantification_level": "ion",
            "axis": {
                "obs_keys": ["Run"],
                "var_keys": ["Protein_Ids"],
                "x_layer": "Precursor_Normalised",
            },
            "columns": {
                "obs": {"select": {"Run": "Run"}},
                "var": {"select": {"Protein_Ids": "Protein.Ids"}},
            },
            "layers": [{"name": "Precursor_Normalised", "source": "Precursor.Normalised"}],
        }
    )


def _adata(*, sparse_x: bool = False) -> ad.AnnData:
    values = np.asarray(
        [
            [10, 20, 5, 10, 10, 10],
            [10, 20, 5, 10, 10, 0],
            [10, 10, 20, 10, 10, np.nan],
            [10, 10, 20, 10, 10, -1],
        ],
        dtype=np.float64,
    )
    matrix = sparse.csr_matrix(np.nan_to_num(values, nan=0.0)) if sparse_x else values
    obs = pd.DataFrame(
        {"Run": ["run_A1", "run_A2", "run_B1", "run_B2"]},
        index=["run_A1", "run_A2", "run_B1", "run_B2"],
    )
    var = pd.DataFrame(
        {
            "Protein_Ids": [
                "P1_HUMAN",
                "P2_YEAST",
                "P3_ECOLI",
                "Cont_P4_HUMAN",
                "P5_HUMAN_YEAST",
                "P6_UNKNOWN",
            ],
            "ProForma_ion": ["H/2", "Y/2", "E/2", "C/2", "M/2", "N/2"],
            "ProForma_peptide": ["H", "Y", "E", "C", "M", "N"],
        },
        index=[f"ion:{name}" for name in ["H/2", "Y/2", "E/2", "C/2", "M/2", "N/2"]],
    )
    result = ad.AnnData(X=matrix, obs=obs, var=var)
    result.uns["anndata_proteomics"] = {
        "quantification_level": "ion",
        "software_name": "DIA-NN",
        "rule_json": json.dumps(_rule().model_dump(mode="json", by_alias=True)),
    }
    return result


def _intermediate(adata: ad.AnnData):
    module = _module_settings()
    tool = _tool_settings()
    rule, roles = resolve_roles(adata, module, tool)
    design = align_runs(adata, rule, roles, module, tool)
    return compute_intermediate(adata, module, tool, roles, design)


def test_matrix_intermediate_matches_hand_computed_values() -> None:
    result = _intermediate(_adata())
    frame = result.varm

    assert frame.index.tolist() == [
        "ion:H/2",
        "ion:Y/2",
        "ion:E/2",
        "ion:C/2",
        "ion:M/2",
        "ion:N/2",
    ]
    assert frame.loc["ion:H/2", "CV_A"] == 0
    assert frame.loc["ion:Y/2", "log2_A_vs_B"] == pytest.approx(1.0)
    assert frame.loc["ion:E/2", "log2_A_vs_B"] == pytest.approx(-2.0)
    assert frame.loc["ion:E/2", "epsilon"] == pytest.approx(0.0)
    assert frame["nr_observed"].tolist() == [4, 4, 4, 4, 4, 1]
    assert frame["included"].tolist() == [True, True, True, False, False, False]
    assert result.legacy["precursor ion"].tolist() == ["E/2", "H/2", "Y/2"]
    assert result.legacy.index.tolist() == [0, 1, 3]


def test_dense_and_sparse_intermediates_are_equal() -> None:
    dense = _intermediate(_adata())
    sparse_result = _intermediate(_adata(sparse_x=True))
    pd.testing.assert_frame_equal(dense.varm, sparse_result.varm)
    pd.testing.assert_frame_equal(dense.legacy, sparse_result.legacy)


def test_legacy_feature_rendering_does_not_change_canonical_proforma() -> None:
    features = pd.Series(["AC[UNIMOD:4]DM[UNIMOD:35]/2", "PEPTIDE/3"])

    rendered = render_proteobench_features(features)

    assert rendered.tolist() == ["AC[Carbamidomethyl]DM[Oxidation]/2", "PEPTIDE/3"]
    assert features.tolist() == ["AC[UNIMOD:4]DM[UNIMOD:35]/2", "PEPTIDE/3"]


def test_legacy_feature_rendering_can_match_before_aa_false_parser() -> None:
    features = pd.Series(
        [
            "PEPTIDEC[UNIMOD:4]/2",
            "PEPTIDEC-[UNIMOD:4]/2",
            "M[UNIMOD:35]PEPTIDE/2",
        ]
    )

    rendered = render_proteobench_features(
        features,
        drop_final_residue_modifications=True,
    )

    assert rendered.tolist() == [
        "PEPTIDEC/2",
        "PEPTIDEC/2",
        "M[Oxidation]PEPTIDE/2",
    ]


def test_raw_modification_parser_preserves_proteobench_tool_semantics() -> None:
    target = _adata()
    raw = pd.Series(["n[42.0106]PEPTIDEC[57.0215]", "M[15.9949]PEPTIDE"])
    canonical = pd.Series(["[UNIMOD:1]-PEPTIDEC[UNIMOD:4]/2", "M[UNIMOD:35]PEPTIDE/3"])
    settings = ModificationParserSettings(
        parse_column="Modified Sequence",
        before_aa=False,
        pattern=r"(?<=\[).+?(?=\])",
        modification_dict={
            "42.0106": "Acetyl",
            "57.0215": "Carbamidomethyl",
            "15.9949": "Oxidation",
        },
    )

    rendered = parse_proteobench_features(
        raw,
        canonical,
        settings,
        level="ion",
        target=target,
    )

    assert rendered.tolist() == ["[Acetyl]-PEPTIDEC/2", "M[Oxidation]PEPTIDE/3"]


def test_scores_keep_proteobench_names_and_thresholds() -> None:
    result = _intermediate(_adata())
    scores = build_scores(result.legacy, result.intermediate_hash)

    assert list(scores["results"]) == ["1", "2", "3", "4", "5", "6"]
    assert scores["results"]["1"]["nr_feature"] == 3
    assert scores["results"]["1"]["roc_auc"] == 1.0
    assert scores["results"]["5"]["nr_feature"] == 0
    assert scores["results"]["5"]["CV_median"] is None
    assert scores["nr_feature"] == 3
    assert scores["intermediate_hash"] == result.intermediate_hash


def test_roc_auc_handles_ties_and_missing_class() -> None:
    tied = pd.DataFrame(
        {
            "species": ["HUMAN", "HUMAN", "YEAST", "YEAST"],
            "log2_A_vs_B": [0.0, 1.0, 0.0, 1.0],
            "log2_expectedRatio": [0.0, 0.0, 1.0, 1.0],
        }
    )
    assert compute_roc_auc(tied) == 0.5
    assert np.isnan(compute_roc_auc(tied[tied["species"] == "HUMAN"]))


def test_pipeline_stores_nested_scores_and_preserves_other_enrichments(
    tmp_path: Path,
) -> None:
    adata = _adata()
    adata.obs["condition"] = ["A", "A", "B", "B"]
    adata.varm["fasta"] = pd.DataFrame({"accession": ["x"] * adata.n_vars}, index=adata.var_names)

    returned = score_quantification(adata, _module_settings(), _tool_settings())
    assert returned is adata
    assert "proteobench" in adata.varm
    assert adata.uns["proteobench"]["column_roles"]["Proteins"] == "var:Protein_Ids"
    assert adata.uns["proteobench"]["protein_mapping"]["species_mapper"] == {
        "_YEAST": "YEAST",
        "_ECOLI": "ECOLI",
        "_HUMAN": "HUMAN",
    }
    accession_mapping = adata.uns["proteobench"]["protein_mapping"]["accession_mapper"]
    assert accession_mapping["entries"] > 30_000
    assert accession_mapping["sha256"] == (
        "032034e2f9bea3fc41290c7461417280b1d37cec41ee8ef9a44c250781a4b997"
    )
    assert adata.uns["proteobench"]["scores"]["nr_feature"] == 3
    assert "fasta" in adata.varm
    assert adata.obs["condition"].tolist() == ["A", "A", "B", "B"]

    path = tmp_path / "scored.h5ad"
    adata.write_h5ad(path)
    restored = ad.read_h5ad(path)
    assert restored.uns["proteobench"]["scores"]["results"]["1"]["nr_feature"] == 3
    assert (
        restored.uns["proteobench"]["protein_mapping"]["accession_mapper"]["sha256"]
        == accession_mapping["sha256"]
    )
    restored_proteobench = restored.varm["proteobench"]
    assert isinstance(restored_proteobench, pd.DataFrame)
    assert restored_proteobench.index.equals(restored.var_names)


def test_pipeline_refuses_collisions_and_conflicting_annotations() -> None:
    adata = _adata()
    score_quantification(adata, _module_settings(), _tool_settings())
    with pytest.raises(ValueError, match="refusing to overwrite"):
        score_quantification(adata, _module_settings(), _tool_settings())

    conflicting = _adata()
    conflicting.obs["condition"] = ["B", "A", "B", "B"]
    with pytest.raises(ValueError, match="disagrees"):
        score_quantification(conflicting, _module_settings(), _tool_settings())


def test_annotation_fasta_and_scoring_are_order_independent() -> None:
    annotation = AnnotationTable(
        samples=pd.DataFrame(
            {
                "raw_file": ["run_A1", "run_A2", "run_B1", "run_B2"],
                "sample_name": ["A1", "A2", "B1", "B2"],
                "condition": ["A", "A", "B", "B"],
            }
        )
    )
    fasta = ">sp|P1|TEST_HUMAN\nHYECMN\n"

    enrich_first = _adata()
    annotate_obs(enrich_first, annotation)
    validate_peptides_against_fasta(enrich_first, fasta, backend="ahocorapy")
    score_quantification(enrich_first, _module_settings(), _tool_settings())

    score_first = _adata()
    score_quantification(score_first, _module_settings(), _tool_settings())
    annotate_obs(score_first, annotation)
    validate_peptides_against_fasta(score_first, fasta, backend="ahocorapy")

    assert enrich_first.uns["proteobench"]["scores"] == score_first.uns["proteobench"]["scores"]
    enrich_scores = enrich_first.varm["proteobench"]
    score_scores = score_first.varm["proteobench"]
    assert isinstance(enrich_scores, pd.DataFrame)
    assert isinstance(score_scores, pd.DataFrame)
    pd.testing.assert_frame_equal(enrich_scores, score_scores)
    assert "fasta_validation" in enrich_first.varm
    assert "fasta_validation" in score_first.varm
    assert enrich_first.obs["condition"].tolist() == score_first.obs["condition"].tolist()


def test_role_resolution_reports_unretained_protein_source() -> None:
    adata = _adata()
    tool = _tool_settings().model_copy(
        update={"mapper": {"Missing.Proteins": "Proteins", "Run": "Raw file"}}
    )
    with pytest.raises(ValueError, match="does not retain"):
        resolve_roles(adata, _module_settings(), tool)


def test_config_defaults_match_golden_json_projection(tmp_path: Path) -> None:
    module_path = tmp_path / "module.toml"
    module_path.write_text(_module_toml())
    tool_path = tmp_path / "tool.toml"
    tool_path.write_text(_tool_toml())

    module = load_module_settings(module_path)
    tool = load_tool_settings(tool_path)

    assert module.general.default_cutoff_min_feature == 1
    assert module.general.max_nr_observed == 6
    assert tool.source_for("Proteins") == "Protein.Ids"


def test_cli_scores_plain_converted_h5ad_and_describe_exposes_scores(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "converted.h5ad"
    _adata().write_h5ad(input_path)
    module_path = tmp_path / "module.toml"
    module_path.write_text(_module_toml())
    tool_path = tmp_path / "tool.toml"
    tool_path.write_text(_tool_toml())
    output_path = tmp_path / "scored.h5ad"

    result = proteobench_cmd(
        input_path,
        module_path,
        tool_path,
        output=output_path,
    )

    assert result == 0
    restored = ad.read_h5ad(output_path)
    assert restored.uns["proteobench"]["scores"]["nr_feature"] == 3
    assert describe(restored)["proteobench"]["scores"]["nr_feature"] == 3
    assert "proteobench" not in ad.read_h5ad(input_path).uns


def test_cli_scores_only_the_module_modality_in_mudata(tmp_path: Path) -> None:
    ion = _adata()
    protein = _adata()
    ion.var_names = [f"ion:{name}" for name in ion.var_names]
    protein.var_names = [f"prt:{name}" for name in protein.var_names]
    protein.uns["anndata_proteomics"]["quantification_level"] = "protein"
    with mudata.set_options(pull_on_update=False):
        container = MuData({"ion": ion, "protein": protein})
    input_path = tmp_path / "converted.h5mu"
    container.write_h5mu(input_path)
    module_path = tmp_path / "module.toml"
    module_path.write_text(_module_toml())
    tool_path = tmp_path / "tool.toml"
    tool_path.write_text(_tool_toml())
    output_path = tmp_path / "scored.h5mu"

    proteobench_cmd(
        input_path,
        module_path,
        tool_path,
        output=output_path,
    )

    with mudata.set_options(pull_on_update=False):
        restored = mudata.read_h5mu(output_path)
    assert restored.mod["ion"].uns["proteobench"]["scores"]["nr_feature"] == 3
    assert "proteobench" not in restored.mod["protein"].uns


def _module_toml() -> str:
    return """
[species_expected_ratio.YEAST]
A_vs_B = 2.0
[species_expected_ratio.ECOLI]
A_vs_B = 0.25
[species_expected_ratio.HUMAN]
A_vs_B = 1.0
[species_mapper]
_YEAST = "YEAST"
_ECOLI = "ECOLI"
_HUMAN = "HUMAN"
[general]
min_count_multispec = 1
level = "ion"
[[samples]]
raw_file = "run_A1"
sample_name = "A1"
condition = "A"
[[samples]]
raw_file = "run_A2"
sample_name = "A2"
condition = "A"
[[samples]]
raw_file = "run_B1"
sample_name = "B1"
condition = "B"
[[samples]]
raw_file = "run_B2"
sample_name = "B2"
condition = "B"
"""


def _tool_toml() -> str:
    return """
[mapper]
Run = "Raw file"
"Protein.Ids" = "Proteins"
"Precursor.Normalised" = "Intensity"
[general]
contaminant_flag = "Cont_"
decoy_flag = true
"""
