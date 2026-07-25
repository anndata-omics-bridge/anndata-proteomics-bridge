"""Golden ProteoBench compatibility test using the gitignored canonical cache."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from anndata_proteomics.converters.pipeline import convert_level, param_version
from anndata_proteomics.proteobench.config import load_module_settings
from anndata_proteomics.proteobench.intermediate import align_runs, compute_intermediate
from anndata_proteomics.proteobench.metrics import build_scores
from anndata_proteomics.proteobench.resolve import resolve_roles
from anndata_proteomics.readers.dispatch import read_table
from anndata_proteomics.readers.summary import describe
from anndata_proteomics.rules.schema import ParseRule

ROOT = Path(__file__).parents[1]
FIXTURE = (
    ROOT
    / "test_data_download/json_dir/Results_quant_ion_DIA_Astral"
    / "269bb8310dc3f501834aaab5ca6fe72791c426e3"
)
GOLDEN_JSON = (
    ROOT
    / "test_data_download/json_dir/Results_quant_ion_DIA_Astral"
    / "Results_quant_ion_DIA_Astral-main"
    / "269bb8310dc3f501834aaab5ca6fe72791c426e3.json"
)
MODULE_TOML = ROOT / "test_data_download/annotations/dia_astral.toml"
REQUIRED = (FIXTURE / "input_file.parquet", FIXTURE / "param_0..txt", GOLDEN_JSON, MODULE_TOML)
DDA_FIXTURE = (
    ROOT
    / "test_data_download/json_dir/Results_quant_ion_DDA_Astral"
    / "300beac4bd267751972cf484bb1cdee2fda0b3a4"
)
DDA_REQUIRED = (DDA_FIXTURE / "input_file.parquet", DDA_FIXTURE / "param_0..txt")


@dataclass(frozen=True)
class GoldenConversionCase:
    name: str
    collection: str
    fixture_hash: str
    software: str
    input_name: str
    params_name: str
    module_name: str
    protein_role: str
    wide_intensity_suffix: bool = False

    @property
    def fixture(self) -> Path:
        return ROOT / "test_data_download/json_dir" / self.collection / self.fixture_hash

    @property
    def golden_json(self) -> Path:
        return (
            ROOT
            / "test_data_download/json_dir"
            / self.collection
            / f"{self.collection}-main"
            / f"{self.fixture_hash}.json"
        )

    @property
    def required(self) -> tuple[Path, ...]:
        return (
            self.fixture / self.input_name,
            self.fixture / self.params_name,
            self.fixture / "result_performance.csv",
            self.golden_json,
            ROOT / "test_data_download/annotations" / self.module_name,
        )


GOLDEN_CONVERSION_CASES = (
    GoldenConversionCase(
        name="fragpipe-qexactive",
        collection="Results_quant_ion_DDA",
        fixture_hash="b9f217a22df2e15b1144830498fcb65f1dbcff57",
        software="fragpipe",
        input_name="input_file.tsv",
        params_name="param_0..workflow",
        module_name="dda_qexactive.toml",
        protein_role="Protein",
        wide_intensity_suffix=True,
    ),
    GoldenConversionCase(
        name="fragpipe-astral",
        collection="Results_quant_ion_DDA_Astral",
        fixture_hash="54db0c0be0a6f8ac40c38af944212c81466d748f",
        software="fragpipe",
        input_name="input_file.tsv",
        params_name="param_0..workflow",
        module_name="dda_astral.toml",
        protein_role="Protein",
        wide_intensity_suffix=True,
    ),
    GoldenConversionCase(
        name="maxquant-astral",
        collection="Results_quant_ion_DDA_Astral",
        fixture_hash="e5a709af0d04b71d0572029df66753a06195a0f6",
        software="maxquant",
        input_name="input_file.txt",
        params_name="param_0..xml",
        module_name="dda_astral.toml",
        protein_role="Proteins",
    ),
)


@pytest.mark.integration
@pytest.mark.skipif(
    not all(path.exists() for path in REQUIRED),
    reason="canonical DIA-NN ProteoBench cache is absent",
)
def test_diann_dia_astral_uses_precursor_normalised() -> None:
    target = _load_plain_converted_fixture()
    module = load_module_settings(MODULE_TOML)
    rule, roles = resolve_roles(target, module)
    design = align_runs(target, rule, roles, module)

    result = compute_intermediate(target, module, roles, design)
    expected = pd.read_csv(FIXTURE / "result_performance.csv")

    assert rule.axis.x_layer == "Precursor_Normalised"
    assert np.array_equal(
        np.asarray(cast(Any, target.X)),
        np.asarray(cast(Any, target.layers["Precursor_Normalised"])),
        equal_nan=True,
    )
    assert describe(target)["column_mapping"]["X"] == {
        "layer": "Precursor_Normalised",
        "source": "Precursor.Normalised",
        "source_kind": "column",
    }
    assert result.legacy.columns.tolist() == expected.columns.tolist()
    assert result.legacy["precursor ion"].is_unique

    golden = json.loads(GOLDEN_JSON.read_text(encoding="utf-8"))
    scores = build_scores(
        result.legacy,
        result.intermediate_hash,
        default_cutoff=module.general.default_cutoff_min_feature,
        max_nr_observed=module.general.max_nr_observed,
    )
    assert scores["nr_feature"] == 111_679
    assert scores["proteobench_version"] == golden["proteobench_version"]
    assert len(scores["intermediate_hash"]) == 40
    assert scores["results"].keys() == golden["results"].keys()
    assert {threshold: values["nr_feature"] for threshold, values in scores["results"].items()} == {
        "1": 111_679,
        "2": 100_607,
        "3": 94_507,
        "4": 87_871,
        "5": 80_570,
        "6": 70_101,
    }


@pytest.mark.integration
@pytest.mark.skipif(
    not all(path.exists() for path in DDA_REQUIRED),
    reason="canonical DIA-NN DDA ProteoBench cache is absent",
)
def test_diann_dda_astral_uses_ms1_normalised() -> None:
    data = read_table(DDA_FIXTURE / "input_file.parquet")
    params = DDA_FIXTURE / "param_0..txt"

    target = convert_level(
        data,
        "diann",
        "ion",
        param_version(params, "diann"),
        params_path=params,
    )

    metadata = target.uns["anndata_proteomics"]
    rule = ParseRule.model_validate_json(metadata["rule_json"])
    assert rule.axis.x_layer == "Ms1_Normalised"
    assert np.array_equal(
        np.asarray(cast(Any, target.X)),
        np.asarray(cast(Any, target.layers["Ms1_Normalised"])),
        equal_nan=True,
    )
    assert describe(target)["column_mapping"]["X"] == {
        "layer": "Ms1_Normalised",
        "source": "Ms1.Normalised",
        "source_kind": "column",
    }


@pytest.mark.integration
@pytest.mark.parametrize(
    "case",
    GOLDEN_CONVERSION_CASES,
    ids=lambda case: case.name,
)
def test_conversion_protein_completion_matches_proteobench_golden(
    case: GoldenConversionCase,
) -> None:
    if not all(path.exists() for path in case.required):
        pytest.skip(f"canonical {case.name} ProteoBench cache is absent")

    params = case.fixture / case.params_name
    data = read_table(case.fixture / case.input_name)
    target = convert_level(
        data,
        case.software,
        "ion",
        param_version(params, case.software),
        params_path=params,
    )
    _assert_representative_completed_protein(target, data, case)
    module = load_module_settings(ROOT / "test_data_download/annotations" / case.module_name)
    rule, roles = resolve_roles(target, module)
    assert roles.proteins == case.protein_role

    design = align_runs(target, rule, roles, module)
    result = compute_intermediate(target, module, roles, design)
    expected_intermediate = pd.read_csv(case.fixture / "result_performance.csv")
    if case.wide_intensity_suffix:
        expected_intermediate = expected_intermediate.rename(
            columns={
                column: column.removesuffix(" Intensity")
                for column in expected_intermediate.columns
            }
        )
    assert result.legacy.columns.tolist() == expected_intermediate.columns.tolist()
    assert result.legacy["precursor ion"].is_unique

    golden = json.loads(case.golden_json.read_text(encoding="utf-8"))
    scores = build_scores(
        result.legacy,
        result.intermediate_hash,
        default_cutoff=module.general.default_cutoff_min_feature,
        max_nr_observed=module.general.max_nr_observed,
    )
    if case.name == "fragpipe-qexactive":
        assert {
            threshold: values["nr_feature"] for threshold, values in scores["results"].items()
        } == {
            "1": 53_074,
            "2": 42_828,
            "3": 35_274,
            "4": 28_557,
            "5": 22_626,
            "6": 15_846,
        }
    else:
        _assert_scores_close(scores, golden)


def _load_plain_converted_fixture() -> ad.AnnData:
    data = read_table(FIXTURE / "input_file.parquet")
    params = FIXTURE / "param_0..txt"
    version = param_version(params, "diann")
    return convert_level(data, "diann", "ion", version, params_path=params)


def _assert_score_mapping_close(actual: dict[str, Any], expected: dict[str, Any]) -> None:
    assert actual.keys() == expected.keys()
    for key, expected_value in expected.items():
        actual_value = actual[key]
        if isinstance(expected_value, dict):
            _assert_score_mapping_close(actual_value, expected_value)
        elif isinstance(expected_value, int) and not isinstance(expected_value, bool):
            assert actual_value == expected_value
        elif isinstance(expected_value, float):
            assert actual_value == pytest.approx(expected_value, abs=2e-5)
        else:
            assert actual_value == expected_value


def _assert_scores_close(actual: dict[str, Any], expected: dict[str, Any]) -> None:
    assert actual["nr_feature"] == expected["nr_feature"]
    assert actual["proteobench_version"] == expected["proteobench_version"]
    assert len(actual["intermediate_hash"]) == 40
    assert actual["results"].keys() == expected["results"].keys()
    _assert_score_mapping_close(actual["results"], expected["results"])
    for key in (
        "median_abs_epsilon_global",
        "mean_abs_epsilon_global",
        "median_abs_epsilon_eq_species",
        "mean_abs_epsilon_eq_species",
        "median_abs_epsilon_precision_global",
        "mean_abs_epsilon_precision_global",
        "median_abs_epsilon_precision_eq_species",
        "mean_abs_epsilon_precision_eq_species",
    ):
        assert actual[key] == pytest.approx(expected[key], abs=2e-5)


def _assert_representative_completed_protein(
    target: ad.AnnData,
    raw: pd.DataFrame,
    case: GoldenConversionCase,
) -> None:
    if case.software == "fragpipe":
        mapped = raw["Mapped Proteins"]
        row = raw.loc[mapped.notna() & mapped.astype(str).ne("")].iloc[0]
        expected = f"{row['Protein']},{row['Mapped Proteins']}"
        modified_column = "Modified_Sequence"
        raw_modified_column = "Modified Sequence"
        charge_column = "Charge"
    else:
        row = raw.loc[raw["Proteins"].isna() & raw["Leading proteins"].notna()].iloc[0]
        expected = row["Leading proteins"]
        modified_column = "Modified_Sequence"
        raw_modified_column = "Modified sequence"
        charge_column = "Charge"

    matching = target.var[
        target.var[modified_column].astype(str).eq(str(row[raw_modified_column]))
        & target.var[charge_column].astype(str).eq(str(row[charge_column]))
    ]
    assert expected in matching[case.protein_role].tolist()
