"""Golden ProteoBench compatibility test using the gitignored canonical cache."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from anndata_proteomics.converters.pipeline import convert_level, param_version
from anndata_proteomics.proteobench.config import load_module_settings, load_tool_settings
from anndata_proteomics.proteobench.intermediate import align_runs, compute_intermediate
from anndata_proteomics.proteobench.metrics import build_scores
from anndata_proteomics.proteobench.resolve import resolve_roles
from anndata_proteomics.readers.dispatch import read_table

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
TOOL_TOML = ROOT / "tests/data/proteobench/parse_settings_diann.toml"
CONVERTED = (
    ROOT.parent / "apb_studio/apb_outputs/Results_quant_ion_DIA_Astral" / "diann-269bb831/ion.h5ad"
)
REQUIRED = (FIXTURE / "input_file.parquet", FIXTURE / "param_0..txt", GOLDEN_JSON, MODULE_TOML)


@pytest.mark.integration
@pytest.mark.skipif(
    not all(path.exists() for path in REQUIRED),
    reason="canonical DIA-NN ProteoBench cache is absent",
)
def test_diann_astral_intermediate_and_scores_match_golden() -> None:
    target = _load_plain_converted_fixture()
    module = load_module_settings(MODULE_TOML)
    tool = load_tool_settings(TOOL_TOML)
    rule, roles = resolve_roles(target, module, tool)
    design = align_runs(target, rule, roles, module, tool)

    result = compute_intermediate(target, module, tool, roles, design)
    expected = pd.read_csv(FIXTURE / "result_performance.csv")

    assert result.legacy.columns.tolist() == expected.columns.tolist()
    assert result.legacy.shape == expected.shape
    feature = "precursor ion"
    assert result.legacy[feature].tolist() == expected[feature].tolist()
    for column in expected.columns.drop(feature):
        if expected[column].dtype.kind in "biufc":
            np.testing.assert_allclose(
                result.legacy[column].to_numpy(dtype=float),
                expected[column].to_numpy(dtype=float),
                rtol=5e-5,
                atol=5e-4,
                equal_nan=True,
            )
        else:
            assert result.legacy[column].tolist() == expected[column].tolist()

    golden = json.loads(GOLDEN_JSON.read_text(encoding="utf-8"))
    scores = build_scores(
        result.legacy,
        result.intermediate_hash,
        default_cutoff=module.general.default_cutoff_min_feature,
        max_nr_observed=module.general.max_nr_observed,
    )
    assert scores["nr_feature"] == golden["nr_feature"]
    assert scores["proteobench_version"] == golden["proteobench_version"]
    assert len(scores["intermediate_hash"]) == 40
    assert scores["results"].keys() == golden["results"].keys()
    assert {threshold: values["nr_feature"] for threshold, values in scores["results"].items()} == {
        threshold: values["nr_feature"] for threshold, values in golden["results"].items()
    }
    _assert_score_mapping_close(scores["results"], golden["results"])
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
        assert scores[key] == pytest.approx(golden[key], abs=2e-5)


def _load_plain_converted_fixture() -> ad.AnnData:
    if CONVERTED.exists():
        return ad.read_h5ad(CONVERTED)
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
