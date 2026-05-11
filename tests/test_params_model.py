"""Tests for ``anndata_proteomics.params.model.Parameters``."""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

from anndata_proteomics.params.model import Parameters


def test_construct_empty_has_all_none():
    p = Parameters()
    dumped = p.model_dump()
    assert dumped["software_name"] is None
    assert dumped["enzyme"] is None
    assert dumped["scan_window"] is None


def test_extra_fields_are_preserved():
    p = Parameters(software_name="Sage", vendor_specific_thing=42)
    assert p.software_name == "Sage"
    assert p.model_dump()["vendor_specific_thing"] == 42


def test_to_series_roundtrip():
    p = Parameters(
        software_name="Sage",
        software_version="0.14.6",
        enzyme="KR",
        min_peptide_length=7,
        max_peptide_length=50,
    )
    s = p.to_series()
    assert s["software_name"] == "Sage"
    assert s["min_peptide_length"] == 7
    assert s["enzyme"] == "KR"
    assert s["abundance_normalization_ions"] is None


def test_from_series_treats_empty_and_nan_as_none():
    s = pd.Series(
        {
            "software_name": "Sage",
            "enzyme": "",
            "max_peptide_length": math.nan,
            "fixed_mods": "{'C': 57.02146}",
        }
    )
    p = Parameters.from_series(s)
    assert p.software_name == "Sage"
    assert p.enzyme is None
    assert p.max_peptide_length is None
    assert p.fixed_mods == "{'C': 57.02146}"


def test_from_series_preserves_literal_none_string():
    # ProteoBench writes literal "None" strings for unset fields; preserve them.
    s = pd.Series({"software_name": "Sage", "ident_fdr_psm": "None"})
    p = Parameters.from_series(s)
    assert p.ident_fdr_psm == "None"


PROTEOBENCH_PARAMS = Path("/Users/wolski/projects/anndata_bridge/ProteoBench/test/params")


def test_from_series_reads_sage_expected_csv():
    csv = PROTEOBENCH_PARAMS / "sage_parameterfile.csv"
    if not csv.exists():
        return  # ProteoBench fixture missing; skip without xfail noise.
    df = pd.read_csv(csv, header=0, index_col=0)
    series = df.iloc[:, 0]
    p = Parameters.from_series(series)
    assert p.software_name == "Sage"
    assert p.enzyme == "KR"
    assert str(p.min_peptide_length) == "7"
