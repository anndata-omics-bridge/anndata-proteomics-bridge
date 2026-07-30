"""Cached-corpus regressions for FragPipe and PEAKS version coverage."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import pytest

from anndata_proteomics._matrix_types import named_layers
from anndata_proteomics.converters.pipeline import (
    convert_level,
    resolve_parameters,
)
from anndata_proteomics.params.model import Parameters
from anndata_proteomics.readers.dispatch import read_table

ROOT = Path(__file__).resolve().parents[1] / "test_data_download" / "json_dir"


@dataclass(frozen=True)
class CorpusCase:
    """One cached vendor input whose conversion contract is pinned here."""

    slug: str
    relative_input: str
    expected_version_status: str

    @property
    def input_path(self) -> Path:
        return ROOT / self.relative_input


CASES = [
    CorpusCase(
        "fragpipe",
        "Results_quant_ion_DDA_Astral/54db0c0be0a6f8ac40c38af944212c81466d748f/input_file.tsv",
        "present",
    ),
    CorpusCase(
        "fragpipe",
        "Results_quant_ion_DDA_Astral/0c36dceb48ce85c9a56fb3d30c4f65ace4ed0aaf/input_file.tsv",
        "present",
    ),
    CorpusCase(
        "fragpipe",
        "Results_quant_ion_DDA/b9f217a22df2e15b1144830498fcb65f1dbcff57/input_file.tsv",
        "present",
    ),
    CorpusCase(
        "fragpipe",
        "Results_quant_ion_DDA/45486140efcbe205e2485f1ef4d668ec3d79fb99/input_file.txt",
        "present",
    ),
    CorpusCase(
        "fragpipe",
        "Results_quant_ion_DDA/28b0c3b9853a5b60c9e47428b8a51b4898083523/input_file.tsv",
        "present",
    ),
    CorpusCase(
        "peaks",
        "Results_quant_ion_DDA/9d1361331b165d6cc779ccf614419eb77057f573/input_file.csv",
        "present",
    ),
    CorpusCase(
        "peaks",
        "Results_quant_ion_DIA_diaPASEF/5691d485356c1abcd8efd2320a1666752a19f50b/input_file.csv",
        "missing",
    ),
    CorpusCase(
        "peaks",
        "Results_quant_ion_DIA_diaPASEF/806987cbdec8eb026bb422b100759e4191bed9d5/input_file.csv",
        "present",
    ),
    CorpusCase(
        "peaks",
        "Results_quant_ion_DIA_AIF/b5fddd9b5d27918e8d31ec07bcf599cbd214027a/input_file.txt",
        "missing",
    ),
    CorpusCase(
        "peaks",
        "Results_quant_ion_DIA_ZenoTOF/d08f2cb6e471eaa24bddaef80ee04cdbfb4f8d79/input_file.csv",
        "present",
    ),
    CorpusCase(
        "peaks",
        "Results_quant_ion_DIA_Astral/b342522520beee8fa15304c6683eed7a2f68e01f/input_file.csv",
        "present",
    ),
    CorpusCase(
        "peaks",
        "Results_quant_ion_DIA_Astral/aa1d53e859aae10d0f5b991862bc2dedbf2f79f6/input_file.csv",
        "present",
    ),
]


@pytest.mark.parametrize(
    "case",
    CASES,
    ids=lambda case: f"{case.slug}-{case.input_path.parent.name[:8]}",
)
def test_cached_fragpipe_and_peaks_inputs_convert_with_exact_run_axis(
    case: CorpusCase,
) -> None:
    if not case.input_path.exists():
        pytest.skip(f"cached ProteoBench input is absent: {case.input_path}")
    param_paths = sorted(case.input_path.parent.glob("param_0.*"))
    if not param_paths:
        pytest.skip(f"cached parameter file is absent beside {case.input_path}")
    param_path = param_paths[0]
    resolution = resolve_parameters(param_path, case.slug)
    frame = read_table(case.input_path)

    result = convert_level(
        frame,
        case.slug,
        "ion",
        resolution.version,
        params_path=param_path,
        parameter_resolution=resolution,
    )

    assert result.n_obs == 6
    assert all(name.startswith("LFQ_") for name in result.obs_names)
    assert not any(
        name.startswith(("Group ", "Condition ")) or name == "Best" for name in result.obs_names
    )
    assert not any(name.endswith((".raw", "_raw")) for name in result.obs_names)
    metadata = result.uns["anndata_proteomics"]
    assert metadata["search_parameters_version_status"] == case.expected_version_status
    expected_method = "columns" if case.expected_version_status == "missing" else "software_version"
    assert metadata["rule_selection_method"] == expected_method
    stored = json.loads(metadata["search_parameters"])
    assert set(stored) == set(Parameters.model_fields)
    if case.expected_version_status == "missing":
        assert stored["software_version"] is None
    if case.slug == "fragpipe":
        assert set(named_layers(result)) == {
            "Intensity",
            "Spectral_Count",
            "Apex_Retention_Time",
            "Match_Type",
        }
    else:
        assert {"Normalized_Area", "Sample_Mz", "Sample_RT_Mean"} <= set(named_layers(result))


def test_peaks_aif_txt_is_detected_as_comma_delimited() -> None:
    input_path = (
        ROOT
        / "Results_quant_ion_DIA_AIF"
        / "b5fddd9b5d27918e8d31ec07bcf599cbd214027a"
        / "input_file.txt"
    )
    if not input_path.exists():
        pytest.skip(f"cached ProteoBench input is absent: {input_path}")

    assert len(read_table(input_path).columns) == 38
