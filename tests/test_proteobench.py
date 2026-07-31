"""ProteoBench matrix scoring and storage contracts."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, cast

import anndata as ad
import mudata
import numpy as np
import pandas as pd
import pytest
from description_support import describe_anndata as describe
from mudata import MuData
from numpy.typing import NDArray
from pydantic import ValidationError
from scipy import sparse

from anndata_proteomics.adapters.anndata import fasta as fasta_adapter
from anndata_proteomics.adapters.anndata import proteobench as proteobench_adapter
from anndata_proteomics.adapters.anndata.annotation import (
    read_observation_frames,
    write_sample_annotation,
)
from anndata_proteomics.adapters.anndata.proteobench import (
    resolve_roles,
)
from anndata_proteomics.annotation.loader import (
    AnnotationTable,
    InMemoryAnnotationOrigin,
)
from anndata_proteomics.proteobench import intermediate, mapping, metrics
from anndata_proteomics.proteobench.config import (
    ExpectedRatio,
    ModuleGeneral,
    ModuleSettings,
    SampleSettings,
    load_module_settings,
)
from anndata_proteomics.proteobench.contracts import QuantMatrix
from anndata_proteomics.proteobench.intermediate import align_runs, compute_intermediate
from anndata_proteomics.proteobench.mapping import (
    render_proteobench_features,
)
from anndata_proteomics.proteobench.metrics import ScoreConfig, build_scores, compute_roc_auc
from anndata_proteomics.proteobench.pipeline import ProteoBenchResult
from anndata_proteomics.rules.schema import ParseRule, QuantificationLevel
from anndata_proteomics.scripts.cli import proteobench as proteobench_cmd
from anndata_proteomics.workflows import fasta as fasta_workflow
from anndata_proteomics.workflows import proteobench as proteobench_workflow
from anndata_proteomics.workflows.proteobench import ProteoBenchLevelInput
from anndata_proteomics.workflows.sample_annotation import run_sample_annotation

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


def _matrix_values() -> NDArray[np.float64]:
    return np.asarray(
        [
            [10, 20, 5, 10, 10, 10],
            [10, 20, 5, 10, 10, 0],
            [10, 10, 20, 10, 10, np.nan],
            [10, 10, 20, 10, 10, -1],
        ],
        dtype=np.float64,
    )


def _observations(*, annotated: bool) -> pd.DataFrame:
    obs_data = {"Run": ["run_A1", "run_A2", "run_B1", "run_B2"]}
    if annotated:
        obs_data.update(
            {
                "sample_name": ["A1", "A2", "B1", "B2"],
                "condition": ["A", "A", "B", "B"],
            }
        )
    return pd.DataFrame(obs_data, index=["run_A1", "run_A2", "run_B1", "run_B2"])


def _features() -> pd.DataFrame:
    return pd.DataFrame(
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
        index=["H/2", "Y/2", "E/2", "C/2", "M/2", "N/2"],
    )


def _adata(*, sparse_x: bool = False, annotated: bool = True) -> ad.AnnData:
    values = _matrix_values()
    matrix = sparse.csr_matrix(np.nan_to_num(values, nan=0.0)) if sparse_x else values
    result = ad.AnnData(X=matrix, obs=_observations(annotated=annotated), var=_features())
    result.uns["anndata_proteomics"] = {
        "quantification_level": "ion",
        "software_name": "DIA-NN",
        "rule_json": json.dumps(_rule().model_dump(mode="json", by_alias=True)),
    }
    return result


def _score_quantification(
    obj: ad.AnnData | MuData,
    module_settings: ModuleSettings,
) -> ad.AnnData | MuData:
    """Compose the AnnData adapter with the backend-independent scoring workflow."""
    targets = tuple(proteobench_adapter.resolve_targets(obj))
    extracted = tuple(proteobench_adapter.read_level(target) for target in targets)
    results = proteobench_workflow.score_levels(
        tuple(level.calculation for level in extracted),
        module_settings,
    )
    for target, level, result in zip(targets, extracted, results, strict=True):
        proteobench_adapter.store_result(target, result, level.roles)
    return obj


def _plain_intermediate(
    matrix: QuantMatrix,
    level: QuantificationLevel,
) -> intermediate.IntermediateResult:
    module = _module_settings()
    features = _features()
    design = align_runs(_observations(annotated=True), module)
    return compute_intermediate(
        matrix,
        features.index,
        features["Protein_Ids"],
        module,
        design,
        level,
    )


def test_matrix_intermediate_matches_hand_computed_values() -> None:
    result = _plain_intermediate(_matrix_values(), "ion")
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
    result = _plain_intermediate(_matrix_values(), "ion")
    expected = pd.read_csv(GOLDEN_LEGACY_INTERMEDIATE, index_col=0)

    pd.testing.assert_frame_equal(result.legacy, expected, check_dtype=False)
    assert result.intermediate_hash == GOLDEN_LEGACY_INTERMEDIATE_HASH


def test_workflow_scores_plain_inputs_without_a_storage_container() -> None:
    features = _features()
    (result,) = proteobench_workflow.score_levels(
        (
            ProteoBenchLevelInput(
                observations=_observations(annotated=True),
                matrix=_matrix_values(),
                feature_ids=features.index,
                reported_proteins=features["Protein_Ids"],
                level="ion",
            ),
        ),
        _module_settings(),
    )

    assert result.intermediate.intermediate_hash == GOLDEN_LEGACY_INTERMEDIATE_HASH


def test_in_memory_boundaries_drive_proteobench_workflow() -> None:
    features = _features()
    source = {
        "ion": ProteoBenchLevelInput(
            observations=_observations(annotated=True),
            matrix=_matrix_values(),
            feature_ids=features.index,
            reported_proteins=features["Protein_Ids"],
            level="ion",
        )
    }
    persisted: dict[QuantificationLevel, intermediate.IntermediateResult] = {}

    def read_levels() -> tuple[ProteoBenchLevelInput, ...]:
        return tuple(source.values())

    def write_results(results: tuple[ProteoBenchResult, ...]) -> None:
        for inputs, result in zip(source.values(), results, strict=True):
            persisted[inputs.level] = result.intermediate

    calculated = proteobench_workflow.score_levels(read_levels(), _module_settings())
    write_results(calculated)

    assert persisted["ion"].intermediate_hash == GOLDEN_LEGACY_INTERMEDIATE_HASH


def test_proteobench_adapter_extracts_typed_level_separately() -> None:
    target = _adata()

    extracted = proteobench_adapter.read_level(target)

    assert extracted.calculation.observations is target.obs
    assert extracted.calculation.matrix is target.X
    assert extracted.calculation.feature_ids.equals(target.var_names)
    assert extracted.calculation.reported_proteins.equals(target.var["Protein_Ids"])
    assert extracted.calculation.level == "ion"
    assert extracted.roles.proteins == "Protein_Ids"


def test_proteobench_adapter_persists_typed_result_separately() -> None:
    features = _features()
    (calculated,) = proteobench_workflow.score_levels(
        (
            ProteoBenchLevelInput(
                observations=_observations(annotated=True),
                matrix=_matrix_values(),
                feature_ids=features.index,
                reported_proteins=features["Protein_Ids"],
                level="ion",
            ),
        ),
        _module_settings(),
    )
    target = _adata()

    proteobench_adapter.store_result(
        target,
        calculated,
        proteobench_adapter.ResolvedRoles(proteins="Protein_Ids"),
    )

    stored = target.varm["proteobench"]
    assert isinstance(stored, pd.DataFrame)
    assert stored.equals(calculated.intermediate.varm)
    assert target.uns["anndata_proteomics"]["proteobench"]["scores"]["nr_feature"] == 3


def test_dense_and_sparse_intermediates_are_equal() -> None:
    dense = _plain_intermediate(_matrix_values(), "ion")
    sparse_result = _plain_intermediate(
        sparse.csr_matrix(np.nan_to_num(_matrix_values(), nan=0.0)),
        "ion",
    )
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
    result = _plain_intermediate(_matrix_values(), "ion")
    scores = build_scores(result.legacy, result.intermediate_hash, ScoreConfig())

    assert list(scores.results) == ["1", "2", "3", "4", "5", "6"]
    assert scores.results["1"].root["nr_feature"] == 3
    assert scores.results["1"].root["roc_auc"] == 1.0
    assert scores.results["5"].root["nr_feature"] == 0
    assert np.isnan(scores.results["5"].root["CV_median"])
    assert scores.nr_feature == 3
    assert scores.intermediate_hash == result.intermediate_hash


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

    returned = _score_quantification(adata, _module_settings())
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
    _score_quantification(adata, _module_settings())
    with pytest.raises(ValueError, match="refusing to overwrite"):
        _score_quantification(adata, _module_settings())

    conflicting = _adata()
    conflicting.obs["condition"] = ["B", "A", "B", "B"]
    with pytest.raises(ValueError, match="does not match module condition"):
        _score_quantification(conflicting, _module_settings())


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
        _score_quantification(unannotated, _module_settings())

    annotation_result = run_sample_annotation(
        read_observation_frames(unannotated),
        annotation,
        InMemoryAnnotationOrigin(),
    )
    write_sample_annotation(unannotated, annotation_result)
    _score_quantification(unannotated, _module_settings())
    fasta_config = fasta_adapter.AnnDataPeptideFastaConfig()
    fasta_input = fasta_adapter.read_peptide_anndata_input(unannotated, fasta_config)
    fasta_result = fasta_workflow.validate_peptide_levels(
        (fasta_input,),
        fasta,
        fasta_config.matching,
    )
    fasta_adapter.write_peptide_validation(
        unannotated,
        fasta_result.levels[fasta_input.name],
        fasta_config.sequence_field,
    )
    fasta_adapter.write_anndata_fasta_config(unannotated, fasta_result.fasta_config)

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
    summary = describe(restored)
    summary_proteobench = summary["proteobench"]
    assert isinstance(summary_proteobench, dict)
    summary_scores = summary_proteobench["scores"]
    assert isinstance(summary_scores, dict)
    assert summary_scores["nr_feature"] == 3
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
    roles = proteobench_adapter.ResolvedRoles(proteins="Protein_Ids")
    stored = roles.as_dict()
    assert stored["Sample name"] == "obs:sample_name"
    assert stored["Condition"] == "obs:condition"
    # The feature axis is the var index at every level, not a level-specific var column.
    assert stored["feature"] == "var_names"

    # An AnnData is its own only target; a MuData yields every modality, in order.
    adata = _adata()
    assert proteobench_adapter.resolve_targets(adata) == [adata]
    other = _adata()
    # Distinct feature names keep the MuData global var axis unique.
    other.var_names = [f"other:{name}" for name in other.var_names]
    with mudata.set_options(pull_on_update=False):
        combined = MuData({"a": adata, "b": other})
        assert proteobench_adapter.resolve_targets(combined) == [adata, other]
        with pytest.raises(ValueError, match="no modality to score"):
            proteobench_adapter.resolve_targets(MuData({}))

    missing_rule = _adata()
    missing_rule.uns["anndata_proteomics"]["rule_json"] = {}
    # The namespace owner reports the offending key and the type it actually found.
    with pytest.raises(TypeError, match="rule_json'\\] must be a JSON string; got dict"):
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
        align_runs(proteobench_adapter.extract_observations(unknown), module)

    repeated = _adata()
    repeated.obs["sample_name"] = ["A1", "A1", "B1", "B2"]
    with pytest.raises(ValueError, match="one-to-one"):
        align_runs(proteobench_adapter.extract_observations(repeated), module)

    incomplete = _adata()[[0, 1, 2], :].copy()
    with pytest.raises(ValueError, match="sample alignment is incomplete"):
        align_runs(proteobench_adapter.extract_observations(incomplete), module)

    conflicting = _adata()
    conflicting.obs["condition"] = ["B", "A", "B", "B"]
    with pytest.raises(ValueError, match="does not match module condition"):
        align_runs(proteobench_adapter.extract_observations(conflicting), module)

    missing_annotation = _adata(annotated=False)
    with pytest.raises(ValueError, match="Run 'apb annotate' first"):
        align_runs(proteobench_adapter.extract_observations(missing_annotation), module)

    incomplete_annotation = _adata()
    cast(pd.DataFrame, incomplete_annotation.obs).loc["run_A1", "sample_name"] = None
    with pytest.raises(ValueError, match="complete sample annotation"):
        align_runs(proteobench_adapter.extract_observations(incomplete_annotation), module)

    # var_names *is* the feature axis, so its uniqueness check is the feature-identity check.
    duplicate_names = _adata()
    duplicate_names.var_names = ["duplicate"] * duplicate_names.n_vars
    with pytest.raises(ValueError, match="unique var_names"):
        compute_intermediate(
            proteobench_adapter.extract_quant_matrix(duplicate_names),
            duplicate_names.var_names.copy(),
            duplicate_names.var[roles.proteins].copy(),
            module,
            align_runs(proteobench_adapter.extract_observations(_adata()), module),
            "ion",
        )


def test_legacy_feature_column_follows_the_scored_level() -> None:
    """ProteoBench names its feature column per level; levels it has no module for use the level."""
    assert _plain_intermediate(_matrix_values(), "ion").legacy.columns[0] == "precursor ion"
    assert _plain_intermediate(_matrix_values(), "peptidoform").legacy.columns[0] == "peptidoform"
    assert _plain_intermediate(_matrix_values(), "fragment").legacy.columns[0] == "fragment"
    assert _plain_intermediate(_matrix_values(), "protein").legacy.columns[0] == "protein"


def test_intermediate_low_level_edge_paths() -> None:
    collapsed = intermediate._collapse_positive_matrix(
        np.asarray([[0.0, np.nan]]),
        np.asarray([0, 1], dtype=np.intp),
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


def test_mapping_preserves_missing_assignments_as_explicit_missing_values() -> None:
    mapped = mapping.map_reported_proteins(pd.Series([pd.NA], dtype="string"))

    assert pd.isna(mapped.proteins.iloc[0])
    assert mapped.matched_token_occurrences == 0
    assert mapped.unmatched_token_occurrences == 0


def test_metric_and_pipeline_edge_paths() -> None:
    with pytest.raises(ValueError, match="default_cutoff"):
        ScoreConfig(default_cutoff=7, max_nr_observed=6)
    empty = pd.DataFrame()
    assert np.isnan(compute_roc_auc(empty))
    missing_ratios = pd.DataFrame(columns=["species", "log2_A_vs_B", "log2_expectedRatio"])
    assert np.isnan(compute_roc_auc(missing_ratios))
    with pytest.raises(ValueError, match="unsupported aggregation"):
        metrics._absolute_aggregate("mode")

    collision = _adata()
    collision.uns["anndata_proteomics"]["proteobench"] = {"scores": {}}
    with pytest.raises(ValueError, match=r"\['scores'\]"):
        _score_quantification(collision, _module_settings())
