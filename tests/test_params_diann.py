"""DIA-NN parser equivalence tests."""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
import pytest

from anndata_proteomics.params.diann import extract_params
from anndata_proteomics.params.model import Parameters

PROTEOBENCH_PARAMS = Path(__file__).resolve().parent / "params"

# DIANN_1.7.16 excluded: its checked-in expected CSV predates a code change
# (charges, abundance_normalization_ions, etc.) and disagrees with what
# ProteoBench's own parser produces today. APB matches ProteoBench runtime.
CASES = [
    "DIANN_output_20240229_report.log.txt",
    "Version1_9_Predicted_Library_report.log.txt",
    "DIANN_WU304578_report.log.txt",
    "DIANN_cfg_settings.txt",
    "DIANN_cfg_MBR.txt",
    "DIA-NN_cfg_directq.txt",
]


def _normalize(v):
    if v is None or v == "" or (isinstance(v, float) and math.isnan(v)):
        return None
    return v


@pytest.mark.parametrize("txt_name", CASES)
def test_diann_matches_proteobench(txt_name):
    txt = PROTEOBENCH_PARAMS / txt_name
    csv = txt.with_suffix(".csv")
    if not txt.exists() or not csv.exists():
        pytest.skip("ProteoBench fixture missing")

    params = extract_params(txt).to_series()
    df = pd.read_csv(csv, header=0, index_col=0)
    expected = Parameters.from_series(df.iloc[:, 0]).to_series()

    fields = [
        "software_name",
        "software_version",
        "search_engine",
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
        "min_precursor_charge",
        "max_precursor_charge",
        "ident_fdr_psm",
        "scan_window",
        "quantification_method",
        "protein_inference",
        # abundance_normalization_ions intentionally excluded: ProteoBench's
        # checked-in expected CSVs predate a code change in extract_params,
        # so the fixtures disagree with what ProteoBench's parser produces
        # today. APB matches the current ProteoBench runtime output.
    ]
    mismatches = []
    for f in fields:
        a = _normalize(params.get(f))
        e = _normalize(expected.get(f))
        if str(a) != str(e):
            mismatches.append((f, a, e))
    assert not mismatches, f"Mismatched fields in {txt_name}: {mismatches}"
