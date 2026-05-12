"""FragPipe workflow-file parser equivalence tests."""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
import pytest

from anndata_proteomics.params.fragpipe import extract_params
from anndata_proteomics.params.model import Parameters

PROTEOBENCH_PARAMS = Path("/Users/wolski/projects/anndata_bridge/ProteoBench/test/params")

CASES = [
    "fragpipe.workflow",
    "fragpipe_older.workflow",
    "fragpipe_win_paths.workflow",
    "fragpipe_v22.workflow",
    "fragpipe_fdr_test.workflow",
    "fragpipe-version.workflow",
]


def _expected(name: str) -> Parameters:
    csv = PROTEOBENCH_PARAMS / f"{Path(name).stem}_extracted_params.csv"
    df = pd.read_csv(csv, header=0, index_col=0)
    return Parameters.from_series(df.iloc[:, 0])


def _normalize(v):
    if v is None or v == "" or (isinstance(v, float) and math.isnan(v)):
        return None
    return v


@pytest.mark.parametrize("workflow_name", CASES)
def test_fragpipe_matches_proteobench(workflow_name):
    workflow = PROTEOBENCH_PARAMS / workflow_name
    expected_csv = PROTEOBENCH_PARAMS / f"{Path(workflow_name).stem}_extracted_params.csv"
    if not workflow.exists() or not expected_csv.exists():
        pytest.skip("ProteoBench fixture missing")
    params = extract_params(workflow).model_dump()
    expected = _expected(workflow_name).model_dump()

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
