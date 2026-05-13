"""Tests for ``anndata_proteomics.params.model.Parameters``."""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
import pytest
from pydantic import ValidationError

from anndata_proteomics.params.model import MassTolerance, Parameters, Probability


def test_construct_empty_has_all_none():
    p = Parameters()
    dumped = p.model_dump()
    assert dumped["software_name"] is None
    assert dumped["enzyme"] is None
    assert dumped["scan_window"] is None


def test_extra_fields_are_rejected():
    with pytest.raises(ValidationError):
        Parameters(software_name="Sage", vendor_specific_thing=42)


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
    assert p.fixed_mods[0].source == "{'C': 57.02146}"
    assert p.to_series()["fixed_mods"] == "{'C': 57.02146}"


def test_from_series_treats_literal_none_string_as_none():
    s = pd.Series({"software_name": "Sage", "ident_fdr_psm": "None"})
    p = Parameters.from_series(s)
    assert p.ident_fdr_psm is None


def test_probability_rejects_invalid_values():
    with pytest.raises(ValidationError):
        Probability(value=1.2)


def test_parameters_reject_negative_mz():
    with pytest.raises(ValidationError):
        Parameters(min_precursor_mz=-1)


def test_parameters_reject_invalid_charge():
    with pytest.raises(ValidationError):
        Parameters(min_precursor_charge=0)


def test_parameters_reject_invalid_range_ordering():
    with pytest.raises(ValidationError):
        Parameters(min_precursor_mz=900, max_precursor_mz=300)


def test_mass_tolerance_normalises_signed_range_to_half_width():
    parsed = MassTolerance.parse("[-20 ppm, 20 ppm]")
    assert parsed.mode == "absolute"
    assert parsed.value == 20.0
    assert parsed.unit == "ppm"


def test_mass_tolerance_rejects_asymmetric_range():
    with pytest.raises(ValueError, match="asymmetric"):
        MassTolerance.parse("[-10 ppm, 30 ppm]")


def test_mass_tolerance_accepts_only_ppm_or_da_units():
    assert MassTolerance.parse("20 ppm").unit == "ppm"
    assert MassTolerance.parse("0.5 Da").unit == "Da"
    assert MassTolerance.parse("20 Th").unit == "Da"
    with pytest.raises(ValidationError):
        MassTolerance(mode="absolute", value=20)
    with pytest.raises(ValueError):
        MassTolerance.parse("20 kg")


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
    assert p.min_peptide_length == 7
