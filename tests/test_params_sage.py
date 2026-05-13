"""Sage parser equivalence tests against ProteoBench fixtures."""

from __future__ import annotations

import io
from pathlib import Path

import pandas as pd
import pytest

from anndata_proteomics.params.model import Parameters
from anndata_proteomics.params.sage import extract_params

PROTEOBENCH_PARAMS = Path("/Users/wolski/projects/anndata_bridge/ProteoBench/test/params")
SAGE_PARAMETERFILE = PROTEOBENCH_PARAMS / "sage_parameterfile.json"
SAGE_PARAMETERFILE_CSV = PROTEOBENCH_PARAMS / "sage_parameterfile.csv"
SAGE_RESULTS = PROTEOBENCH_PARAMS / "sage_results.json"
SAGE_RESULTS_CSV = PROTEOBENCH_PARAMS / "sage_results.csv"


def _expected_series(csv: Path) -> pd.Series:
    df = pd.read_csv(csv, header=0, index_col=0)
    return df.iloc[:, 0]


def _assert_matches_expected(params: Parameters, expected_csv: Path) -> None:
    expected = Parameters.from_series(_expected_series(expected_csv))
    e = expected.to_series()
    a = params.to_series()
    fields_to_check = [
        "software_name",
        "software_version",
        "search_engine",
        "search_engine_version",
        "enzyme",
        "semi_enzymatic",
        "allowed_miscleavages",
        "min_peptide_length",
        "max_peptide_length",
        "max_mods",
        "min_precursor_charge",
        "max_precursor_charge",
        "precursor_mass_tolerance",
        "fragment_mass_tolerance",
        "enable_match_between_runs",
    ]
    mismatches = []
    for field in fields_to_check:
        if str(a.get(field)) != str(e.get(field)):
            mismatches.append((field, a.get(field), e.get(field)))
    assert not mismatches, f"Mismatched fields: {mismatches}"


@pytest.mark.skipif(not SAGE_PARAMETERFILE.exists(), reason="ProteoBench fixture missing")
def test_sage_parameterfile_matches_proteobench_csv():
    params = extract_params(SAGE_PARAMETERFILE)
    _assert_matches_expected(params, SAGE_PARAMETERFILE_CSV)


@pytest.mark.skipif(not SAGE_RESULTS.exists(), reason="ProteoBench fixture missing")
def test_sage_results_matches_proteobench_csv():
    params = extract_params(SAGE_RESULTS)
    _assert_matches_expected(params, SAGE_RESULTS_CSV)


def test_sage_accepts_filelike_object(tmp_path):
    payload = b"""{
        "version": "0.14.6",
        "database": {
            "enzyme": {
                "missed_cleavages": 1, "min_len": 7, "max_len": 50,
                "cleave_at": "KR", "restrict": null, "c_terminal": true,
                "semi_enzymatic": null
            },
            "static_mods": {"C": 57.02146},
            "variable_mods": {"M": [15.9949]},
            "max_variable_mods": 3
        },
        "precursor_tol": {"ppm": [-20.0, 20.0]},
        "fragment_tol": {"ppm": [-20.0, 20.0]},
        "precursor_charge": [1, 7]
    }"""
    buf = io.BytesIO(payload)
    params = extract_params(buf)
    assert params.software_name == "Sage"
    assert params.software_version == "0.14.6"
    assert params.precursor_mass_tolerance is not None
    assert params.precursor_mass_tolerance.to_legacy() == "20 ppm"
    # restrict is present (null) → enzyme stays raw (ProteoBench parity)
    assert params.enzyme == "KR"


def test_sage_trypsin_p_when_restrict_missing():
    payload = b"""{
        "version": "0.14.6",
        "database": {
            "enzyme": {
                "missed_cleavages": 1, "min_len": 7, "max_len": 50,
                "cleave_at": "KR", "semi_enzymatic": null
            },
            "static_mods": {}, "variable_mods": {}, "max_variable_mods": 3
        },
        "precursor_tol": {"ppm": [-20.0, 20.0]},
        "fragment_tol": {"ppm": [-20.0, 20.0]},
        "precursor_charge": [1, 7]
    }"""
    params = extract_params(io.BytesIO(payload))
    assert params.enzyme == "Trypsin/P"


def test_sage_trypsin_with_p_restrict():
    payload = b"""{
        "version": "0.14.6",
        "database": {
            "enzyme": {
                "missed_cleavages": 1, "min_len": 7, "max_len": 50,
                "cleave_at": "KR", "restrict": "P",
                "semi_enzymatic": null
            },
            "static_mods": {}, "variable_mods": {}, "max_variable_mods": 3
        },
        "precursor_tol": {"ppm": [-20.0, 20.0]},
        "fragment_tol": {"ppm": [-20.0, 20.0]},
        "precursor_charge": [1, 7]
    }"""
    params = extract_params(io.BytesIO(payload))
    assert params.enzyme == "Trypsin"


def test_sage_semi_enzymatic_true():
    payload = b"""{
        "version": "0.14.6",
        "database": {
            "enzyme": {
                "missed_cleavages": 1, "min_len": 7, "max_len": 50,
                "cleave_at": "KR", "semi_enzymatic": true
            },
            "static_mods": {}, "variable_mods": {}, "max_variable_mods": 3
        },
        "precursor_tol": {"ppm": [-20.0, 20.0]},
        "fragment_tol": {"ppm": [-20.0, 20.0]},
        "precursor_charge": [1, 7]
    }"""
    params = extract_params(io.BytesIO(payload))
    assert params.semi_enzymatic is True
