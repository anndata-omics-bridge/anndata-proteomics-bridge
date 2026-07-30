"""FragPipe workflow-file parser equivalence tests."""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
import pytest

from anndata_proteomics.params.model import Parameters
from anndata_proteomics.params.parsers.fragpipe import _read_workflow, extract_params

PROTEOBENCH_PARAMS = Path(__file__).resolve().parent / "params"

CASES = [
    "fragpipe.workflow",
    "fragpipe_older.workflow",
    "fragpipe_win_paths.workflow",
    "fragpipe_v22.workflow",
    "fragpipe_fdr_test.workflow",
    "fragpipe-version.workflow",
    "fragpipe_v23_noMBR.workflow",
]


def _expected(name: str) -> Parameters:
    csv = PROTEOBENCH_PARAMS / f"{Path(name).stem}_extracted_params.csv"
    df = pd.read_csv(csv, header=0, index_col=0)
    return Parameters.from_series(df.iloc[:, 0])


def _normalize(value: object) -> object:
    if value is None or value == "" or (isinstance(value, float) and math.isnan(value)):
        return None
    return value


@pytest.mark.parametrize("workflow_name", CASES)
def test_fragpipe_matches_proteobench(workflow_name: str):
    workflow = PROTEOBENCH_PARAMS / workflow_name
    expected_csv = PROTEOBENCH_PARAMS / f"{Path(workflow_name).stem}_extracted_params.csv"
    if not workflow.exists() or not expected_csv.exists():
        pytest.skip("ProteoBench fixture missing")
    params = extract_params(workflow).to_series()
    expected = _expected(workflow_name).to_series()

    fields = [
        "software_name",
        "software_version",
        "search_engine",
        "search_engine_version",
        "enzyme",
        "semi_enzymatic",
        "allowed_miscleavages",
        "fixed_mods",
        "variable_mods",
        "max_mods",
        "min_peptide_length",
        "max_peptide_length",
        "precursor_mass_tolerance",
        "fragment_mass_tolerance",
        "ident_fdr_psm",
        "ident_fdr_protein",
        "enable_match_between_runs",
        "min_precursor_charge",
        "protein_inference",
    ]
    mismatches = []
    for f in fields:
        a = _normalize(params.get(f))
        e = _normalize(expected.get(f))
        if str(a) != str(e):
            mismatches.append((f, a, e))
    assert not mismatches, f"Mismatched fields in {workflow_name}: {mismatches}"


def test_fragpipe_exposes_embedded_diann_quantification_version() -> None:
    params = extract_params(PROTEOBENCH_PARAMS / "fragpipe.workflow")

    assert params.software_name == "FragPipe"
    assert params.software_version == "23.0"
    assert params.quantification_software == "DIA-NN"
    assert params.quantification_software_version == "1.8.2 beta 8"


def test_fragpipe_without_diann_quantification_has_no_embedded_quantifier() -> None:
    params = extract_params(PROTEOBENCH_PARAMS / "fragpipe_v23_noMBR.workflow")

    assert params.quantification_software is None
    assert params.quantification_software_version is None


def test_fragpipe_reads_diann_version_from_legacy_executable_path() -> None:
    workflow = (
        "# FragPipe (22.0) runtime properties\n"
        "fragpipe-config.bin-diann=C\\:\\\\tools\\\\diann\\\\1.8.2_beta_8\\\\win\\\\DiaNN.exe\n"
    )

    assert _read_workflow(workflow)[3] == "1.8.2 beta 8"
