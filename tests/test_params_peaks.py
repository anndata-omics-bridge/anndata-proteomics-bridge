"""PEAKS parser equivalence tests."""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
import pytest

from anndata_proteomics.params.model import Parameters
from anndata_proteomics.params.peaks import extract_params

PROTEOBENCH_PARAMS = Path("/Users/wolski/projects/anndata_bridge/ProteoBench/test/params")

CASES = [
    "PEAKS_parameters.txt",
    "PEAKS_parameters_DDA.txt",
    "PEAKS_parameters_DIA.txt",
    "PEAKS_parameters_DDA_new.txt",
    "PEAKS_diaPASEF.txt",
]


def _normalize(v):
    if v is None or v == "" or (isinstance(v, float) and math.isnan(v)):
        return None
    return v


@pytest.mark.parametrize("txt_name", CASES)
def test_peaks_matches_proteobench(txt_name):
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
        "search_engine_version",
        "ident_fdr_psm",
        "ident_fdr_peptide",
        "ident_fdr_protein",
        "enable_match_between_runs",
        "precursor_mass_tolerance",
        "fragment_mass_tolerance",
        "enzyme",
        "semi_enzymatic",
        "allowed_miscleavages",
        "min_peptide_length",
        "max_peptide_length",
        "fixed_mods",
        "variable_mods",
        "max_mods",
        "min_precursor_charge",
        "max_precursor_charge",
        "quantification_method",
        "abundance_normalization_ions",
    ]
    mismatches = []
    for f in fields:
        a = _normalize(params.get(f))
        e = _normalize(expected.get(f))
        if str(a) != str(e):
            mismatches.append((f, a, e))
    assert not mismatches, f"Mismatched fields in {txt_name}: {mismatches}"
