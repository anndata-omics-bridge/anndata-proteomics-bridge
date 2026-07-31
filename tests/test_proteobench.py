"""ProteoBench matrix scoring and storage contracts."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast

import anndata as ad
import mudata
import numpy as np
import pandas as pd
import pytest
from mudata import MuData
from pydantic import ValidationError
from scipy import sparse

from anndata_proteomics.annotation.apply import annotate_obs
from anndata_proteomics.annotation.loader import AnnotationTable
from anndata_proteomics.annotation.validate_fasta import (
    FastaValidationConfig,
    validate_peptides_against_fasta,
)
from anndata_proteomics.proteobench import intermediate, mapping, metrics, resolve
from anndata_proteomics.proteobench.config import (
    ExpectedRatio,
    ModuleGeneral,
    ModuleSettings,
    SampleSettings,
    load_module_settings,
)
from anndata_proteomics.proteobench.intermediate import align_runs, compute_intermediate
from anndata_proteomics.proteobench.mapping import (
    render_proteobench_features,
)
from anndata_proteomics.proteobench.metrics import build_scores, compute_roc_auc
from anndata_proteomics.proteobench.pipeline import score_quantification
from anndata_proteomics.proteobench.resolve import resolve_roles
from anndata_proteomics.readers.summary import describe
from anndata_proteomics.rules.schema import ParseRule
from anndata_proteomics.scripts.cli import proteobench as proteobench_cmd

GOLDEN_LEGACY_INTERMEDIATE = (
    Path(__file__).parent / "data" / "proteobench" / "small_legacy_intermediate.txt"
)
GOLDEN_LEGACY_INTERMEDIATE_HASH = "9077847f733c12b1297a4928a0b4c509e50e4ed9"


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
            "column_roles": {
                "protein_assignment": "Protein_Ids",
                "fasta_accessions": "Protein_Ids",
            },
            "layers": [{"name": "Precursor_Normalised", "source": "Precursor.Normalised"}],
        }
    )


def _adata(*, sparse_x: bool = False, annotated: bool = True) -> ad.AnnData:
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
    obs_data = {"Run": ["run_A1", "run_A2", "run_B1", "run_B2"]}
    if annotated:
        obs_data.update(
            {
                "sample_name": ["A1", "A2", "B1", "B2"],
                "condition": ["A", "A", "B", "B"],
            }
        )
    obs = pd.DataFrame(obs_data, index=["run_A1", "run_A2", "run_B1", "run_B2"])
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
        # A real ion rule keeps ProForma_ion in axis.var_keys, so var_names *are* those values.
        index=["H/2", "Y/2", "E/2", "C/2", "M/2", "N/2"],
    )
    result = ad.AnnData(X=matrix, obs=obs, var=var)
    result.uns["anndata_proteomics"] = {
        "quantification_level": "ion",
        "software_name": "DIA-NN",
        "rule_json": json.dumps(_rule().model_dump(mode="json", by_alias=True)),
    }
    return result


def _intermediate(adata: ad.AnnData, *, level: str | None = None):
    module = _module_settings()
    rule, roles = resolve_roles(adata)
    design = align_runs(adata, module)
    return compute_intermediate(
        adata,
        module,
        roles,
        design,
        level=level or rule.quantification_level,
    )


def test_matrix_intermediate_matches_hand_computed_values() -> None:
    result = _intermediate(_adata())
    frame = result.varm

    assert frame.index.tolist() == ["H/2", "Y/2", "E/2", "C/2", "M/2", "N/2"]
    assert frame.loc["H/2", "CV_A"] == 0
    assert frame.loc["Y/2", "log2_A_vs_B"] == pytest.approx(1.0)
    assert frame.loc["E/2", "log2_A_vs_B"] == pytest.approx(-2.0)
    assert frame.loc["E/2", "epsilon"] == pytest.approx(0.0)
    assert frame["nr_observed"].tolist() == [4, 4, 4, 4, 4, 1]
    assert frame["included"].tolist() == [True, True, True, False, False, False]
    assert result.legacy["precursor ion"].tolist() == ["E/2", "H/2", "Y/2"]
    assert result.legacy.index.tolist() == [0, 1, 3]


def test_complete_legacy_intermediate_and_hash_match_golden() -> None:
    result = _intermediate(_adata())
    expected = pd.read_csv(GOLDEN_LEGACY_INTERMEDIATE, index_col=0)

    pd.testing.assert_frame_equal(result.legacy, expected, check_dtype=False)
    assert result.intermediate_hash == GOLDEN_LEGACY_INTERMEDIATE_HASH


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
    adata.varm["fasta"] = pd.DataFrame({"accession": ["x"] * adata.n_vars}, index=adata.var_names)

    returned = score_quantification(adata, _module_settings())
    assert returned is adata
    assert "proteobench" in adata.varm
    assert "proteobench" not in adata.uns
    proteobench = adata.uns["anndata_proteomics"]["proteobench"]
    assert proteobench["column_roles"]["Proteins"] == "var:Protein_Ids"
    assert proteobench["protein_mapping"]["species_mapper"] == {
        "_YEAST": "YEAST",
        "_ECOLI": "ECOLI",
        "_HUMAN": "HUMAN",
    }
    accession_mapping = proteobench["protein_mapping"]["accession_mapper"]
    assert accession_mapping["entries"] > 30_000
    assert accession_mapping["sha256"] == (
        "032034e2f9bea3fc41290c7461417280b1d37cec41ee8ef9a44c250781a4b997"
    )
    assert proteobench["scores"]["nr_feature"] == 3
    assert "fasta" in adata.varm
    assert adata.obs["condition"].tolist() == ["A", "A", "B", "B"]

    path = tmp_path / "scored.h5ad"
    adata.write_h5ad(path)
    restored = ad.read_h5ad(path)
    restored_namespace = restored.uns["anndata_proteomics"]["proteobench"]
    assert restored_namespace["scores"]["results"]["1"]["nr_feature"] == 3
    assert (
        restored_namespace["protein_mapping"]["accession_mapper"]["sha256"]
        == accession_mapping["sha256"]
    )
    restored_proteobench = restored.varm["proteobench"]
    assert isinstance(restored_proteobench, pd.DataFrame)
    assert restored_proteobench.index.equals(restored.var_names)


def test_pipeline_refuses_collisions_and_invalid_annotations() -> None:
    adata = _adata()
    score_quantification(adata, _module_settings())
    with pytest.raises(ValueError, match="refusing to overwrite"):
        score_quantification(adata, _module_settings())

    conflicting = _adata()
    conflicting.obs["condition"] = ["B", "A", "B", "B"]
    with pytest.raises(ValueError, match="does not match module condition"):
        score_quantification(conflicting, _module_settings())


def test_scoring_requires_annotation_while_fasta_remains_independent() -> None:
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

    unannotated = _adata(annotated=False)
    with pytest.raises(ValueError, match="Run 'apb annotate' first"):
        score_quantification(unannotated, _module_settings())

    annotate_obs(unannotated, annotation)
    score_quantification(unannotated, _module_settings())
    validate_peptides_against_fasta(
        unannotated,
        fasta,
        FastaValidationConfig(backend="ahocorapy"),
    )

    assert unannotated.uns["anndata_proteomics"]["proteobench"]["scores"]["nr_feature"] == 3
    assert "fasta_validation" in unannotated.varm


def test_role_resolution_reports_missing_canonical_protein_column() -> None:
    adata = _adata()
    rule_document = _rule().model_dump(mode="json", by_alias=True)
    rule_document["column_roles"]["protein_assignment"] = "Missing_Proteins"
    rule_document["columns"]["var"]["select"]["Missing_Proteins"] = "Missing.Proteins"
    adata.uns["anndata_proteomics"]["rule_json"] = json.dumps(rule_document)
    with pytest.raises(ValueError, match="missing var column"):
        resolve_roles(adata)


def test_config_defaults_match_golden_json_projection(tmp_path: Path) -> None:
    module_path = tmp_path / "module.toml"
    module_path.write_text(_module_toml())

    module = load_module_settings(module_path)

    assert module.general.default_cutoff_min_feature == 1
    assert module.general.max_nr_observed == 6


def test_cli_scores_annotated_h5ad_and_describe_exposes_scores(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "converted.h5ad"
    _adata().write_h5ad(input_path)
    module_path = tmp_path / "module.toml"
    module_path.write_text(_module_toml())
    output_path = tmp_path / "scored.h5ad"

    result = proteobench_cmd(
        input_path,
        module_path,
        output=output_path,
    )

    assert result == 0
    restored = ad.read_h5ad(output_path)
    assert restored.uns["anndata_proteomics"]["proteobench"]["scores"]["nr_feature"] == 3
    assert describe(restored)["proteobench"]["scores"]["nr_feature"] == 3
    assert "proteobench" not in ad.read_h5ad(input_path).uns


def test_cli_rejects_unannotated_h5ad(tmp_path: Path) -> None:
    input_path = tmp_path / "converted.h5ad"
    unannotated = _adata(annotated=False)
    unannotated.write_h5ad(input_path)
    module_path = tmp_path / "module.toml"
    module_path.write_text(_module_toml())

    with pytest.raises(ValueError, match="Run 'apb annotate' first"):
        proteobench_cmd(
            input_path,
            module_path,
            output=tmp_path / "scored.h5ad",
        )


def test_cli_scores_every_modality_in_mudata(tmp_path: Path) -> None:
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
    output_path = tmp_path / "scored.h5mu"

    proteobench_cmd(
        input_path,
        module_path,
        output=output_path,
    )

    with mudata.set_options(pull_on_update=False):
        restored = mudata.read_h5mu(output_path)
    # Scores are per level: every modality keeps its own, and the container itself stays untouched.
    for name in ("ion", "protein"):
        assert (
            restored.mod[name].uns["anndata_proteomics"]["proteobench"]["scores"]["nr_feature"] == 3
        )
        assert "proteobench" in restored.mod[name].varm
    assert "proteobench" not in restored.uns


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


def test_configuration_models_reject_inconsistent_contracts() -> None:
    with pytest.raises(ValidationError, match="must not exceed"):
        ModuleGeneral(
            min_count_multispec=1,
            level="ion",
            default_cutoff_min_feature=4,
            max_nr_observed=3,
        )

    base = _module_settings().model_dump(by_alias=True)
    invalid_documents = [
        {**base, "species_mapper": {"a": "HUMAN", "b": "HUMAN"}},
        {**base, "species_mapper": {"a": "OTHER"}},
        {
            **base,
            "samples": [
                base["samples"][0],
                {**base["samples"][1], "raw_file": base["samples"][0]["raw_file"]},
            ],
        },
        {
            **base,
            "samples": [
                base["samples"][0],
                {**base["samples"][1], "sample_name": base["samples"][0]["sample_name"]},
            ],
        },
        {
            **base,
            "samples": [
                {**sample, "condition": "A"}
                for sample in cast(list[dict[str, Any]], base["samples"])
            ],
        },
    ]
    for document in invalid_documents:
        with pytest.raises(ValidationError):
            ModuleSettings.model_validate(document)


def test_role_target_and_optional_role_validation() -> None:
    roles = resolve.ResolvedRoles(proteins="Protein_Ids")
    stored = roles.as_dict()
    assert stored["Sample name"] == "obs:sample_name"
    assert stored["Condition"] == "obs:condition"
    # The feature axis is the var index at every level, not a level-specific var column.
    assert stored["feature"] == "var_names"

    # An AnnData is its own only target; a MuData yields every modality, in order.
    adata = _adata()
    assert resolve.resolve_targets(adata) == [adata]
    other = _adata()
    assert resolve.resolve_targets(SimpleNamespace(mod={"a": adata, "b": other})) == [adata, other]
    with pytest.raises(ValueError, match="no modality to score"):
        resolve.resolve_targets(SimpleNamespace(mod={}))

    missing_rule = _adata()
    missing_rule.uns["anndata_proteomics"]["rule_json"] = {}
    with pytest.raises(ValueError, match="no string"):
        resolve_roles(missing_rule)

    no_roles = _adata()
    rule_document = _rule().model_dump(mode="json", by_alias=True)
    rule_document["column_roles"] = {}
    no_roles.uns["anndata_proteomics"]["rule_json"] = json.dumps(rule_document)
    with pytest.raises(ValueError, match="protein_assignment"):
        resolve_roles(no_roles)

    no_assignment = _adata()
    rule_document = _rule().model_dump(mode="json", by_alias=True)
    rule_document["column_roles"] = {"fasta_accessions": "Protein_Ids"}
    no_assignment.uns["anndata_proteomics"]["rule_json"] = json.dumps(rule_document)
    with pytest.raises(ValueError, match="protein_assignment"):
        resolve_roles(no_assignment)


def test_alignment_and_intermediate_identity_guards() -> None:
    module = _module_settings()
    _, roles = resolve_roles(_adata())

    unknown = _adata()
    unknown.obs["sample_name"] = ["unknown", "A2", "B1", "B2"]
    with pytest.raises(ValueError, match="does not match"):
        align_runs(unknown, module)

    repeated = _adata()
    repeated.obs["sample_name"] = ["A1", "A1", "B1", "B2"]
    with pytest.raises(ValueError, match="one-to-one"):
        align_runs(repeated, module)

    incomplete = _adata()[[0, 1, 2], :].copy()
    with pytest.raises(ValueError, match="sample alignment is incomplete"):
        align_runs(incomplete, module)

    conflicting = _adata()
    conflicting.obs["condition"] = ["B", "A", "B", "B"]
    with pytest.raises(ValueError, match="does not match module condition"):
        align_runs(conflicting, module)

    missing_annotation = _adata(annotated=False)
    with pytest.raises(ValueError, match="Run 'apb annotate' first"):
        align_runs(missing_annotation, module)

    incomplete_annotation = _adata()
    cast(pd.DataFrame, incomplete_annotation.obs).loc["run_A1", "sample_name"] = None
    with pytest.raises(ValueError, match="complete sample annotation"):
        align_runs(incomplete_annotation, module)

    # var_names *is* the feature axis, so its uniqueness check is the feature-identity check.
    duplicate_names = _adata()
    duplicate_names.var_names = ["duplicate"] * duplicate_names.n_vars
    with pytest.raises(ValueError, match="unique var_names"):
        compute_intermediate(
            duplicate_names,
            module,
            roles,
            align_runs(_adata(), module),
            level="ion",
        )


def test_legacy_feature_column_follows_the_scored_level() -> None:
    """ProteoBench names its feature column per level; levels it has no module for use the level."""
    assert _intermediate(_adata()).legacy.columns[0] == "precursor ion"
    assert _intermediate(_adata(), level="peptidoform").legacy.columns[0] == "peptidoform"
    assert _intermediate(_adata(), level="fragment").legacy.columns[0] == "fragment"
    assert _intermediate(_adata(), level="protein").legacy.columns[0] == "protein"


def test_intermediate_low_level_edge_paths() -> None:
    collapsed = intermediate._collapse_positive_matrix(
        np.asarray([[0.0, np.nan]]),
        np.asarray([0, 1]),
        2,
        np.asarray([True, True], dtype=np.bool_),
        np.float64,
    )
    assert np.isnan(collapsed).all()
    centers = intermediate._empirical_centers(
        np.asarray([1.0]),
        np.asarray(["HUMAN"], dtype=object),
        np.asarray([False], dtype=np.bool_),
    )
    assert np.isnan(centers["median"]).all()
    assert intermediate._contaminants(pd.Series(["P1", "Cont_P2"])).tolist() == [False, True]
    assert intermediate._is_float32_backed(np.asarray([np.nan])) is False


def test_mapping_helpers_cover_empty_tokens() -> None:
    mapped_accession = mapping.map_reported_proteins(pd.Series(["Cont_P00722"]))
    assert mapped_accession.proteins.tolist() == ["sp|Cont_P00722|BGAL_ECOLI"]
    assert mapped_accession.matched_token_occurrences == 1

    mapped = mapping.map_reported_proteins(pd.Series([";P1_UNKNOWN,,"]))
    assert mapped.proteins.tolist() == ["P1_UNKNOWN"]
    assert mapped.unmatched_token_occurrences == 1


def test_metric_and_pipeline_edge_paths() -> None:
    result = _intermediate(_adata())
    with pytest.raises(ValueError, match="default cutoff"):
        build_scores(result.legacy, "hash", default_cutoff=0)
    empty = pd.DataFrame()
    assert np.isnan(compute_roc_auc(empty))
    missing_ratios = pd.DataFrame(columns=["species", "log2_A_vs_B", "log2_expectedRatio"])
    assert np.isnan(compute_roc_auc(missing_ratios))
    with pytest.raises(ValueError, match="unsupported aggregation"):
        metrics._absolute_aggregate("mode")
    assert metrics._json_compatible((np.int64(2), np.float64(np.nan), np.bool_(True))) == [
        2,
        None,
        True,
    ]

    collision = _adata()
    collision.uns["anndata_proteomics"]["proteobench"] = {"scores": {}}
    with pytest.raises(ValueError, match=r"\['scores'\]"):
        score_quantification(collision, _module_settings())
