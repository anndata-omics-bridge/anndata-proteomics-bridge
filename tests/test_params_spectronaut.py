"""Spectronaut parser equivalence tests."""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
import pytest

from anndata_proteomics.params.model import Parameters
from anndata_proteomics.params.spectronaut import extract_params

PROTEOBENCH_PARAMS = Path("/Users/wolski/projects/anndata_bridge/ProteoBench/test/params")

CASES = [
    "spectronaut_Experiment1_ExperimentSetupOverview_BGS_Factory_Settings.txt",
    "Spectronaut_dynamic.txt",
    "Spectronaut_static.txt",
    "Spectronaut_relative.txt",
]


def _normalize(v):
    if v is None or v == "" or (isinstance(v, float) and math.isnan(v)):
        return None
    return v


@pytest.mark.parametrize("txt_name", CASES)
def test_spectronaut_matches_proteobench(txt_name):
    txt = PROTEOBENCH_PARAMS / txt_name
    csv = txt.with_suffix(".csv")
    if not txt.exists() or not csv.exists():
        pytest.skip("ProteoBench fixture missing")

    params = extract_params(txt).model_dump()
    df = pd.read_csv(csv, header=0, index_col=0)
    expected = Parameters.from_series(df.iloc[:, 0]).model_dump()

    fields = [
        "software_name",
        "software_version",
        "search_engine",
        "ident_fdr_psm",
        "ident_fdr_protein",
        "enable_match_between_runs",
        "precursor_mass_tolerance",
        "fragment_mass_tolerance",
        "enzyme",
        "allowed_miscleavages",
        "min_peptide_length",
        "max_peptide_length",
        "fixed_mods",
        "variable_mods",
        "max_mods",
        "quantification_method",
        "protein_inference",
        "abundance_normalization_ions",
    ]
    mismatches = []
    for f in fields:
        a = _normalize(params.get(f))
        e = _normalize(expected.get(f))
        if str(a) != str(e):
            mismatches.append((f, a, e))
    assert not mismatches, f"Mismatched fields in {txt_name}: {mismatches}"
