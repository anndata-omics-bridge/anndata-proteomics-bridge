"""MetaMorpheus parser equivalence tests."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from anndata_proteomics.params.metamorpheus import extract_params
from anndata_proteomics.params.model import Parameters

PROTEOBENCH_PARAMS = Path("/Users/wolski/projects/anndata_bridge/ProteoBench/test/params")
TOML_FILE = PROTEOBENCH_PARAMS / "metamorpheus_search_task_config.toml"
VERSION_FILE = PROTEOBENCH_PARAMS / "metamorpheus_version_result.txt"
EXPECTED_CSV = PROTEOBENCH_PARAMS / "metamorpheus_parameters.csv"


def _expected() -> Parameters:
    df = pd.read_csv(EXPECTED_CSV, header=0, index_col=0)
    return Parameters.from_series(df.iloc[:, 0])


@pytest.mark.skipif(not TOML_FILE.exists(), reason="ProteoBench fixture missing")
def test_metamorpheus_matches_proteobench():
    params = extract_params(TOML_FILE, VERSION_FILE).to_series()
    expected = _expected().to_series()
    fields = [
        "software_name",
        "software_version",
        "search_engine",
        "enzyme",
        "allowed_miscleavages",
        "fixed_mods",
        "variable_mods",
        "precursor_mass_tolerance",
        "fragment_mass_tolerance",
        "min_peptide_length",
        "max_peptide_length",
        "max_mods",
        "min_precursor_charge",
        "max_precursor_charge",
        "enable_match_between_runs",
        "quantification_method",
        "ident_fdr_psm",
    ]
    mismatches = []
    for f in fields:
        if str(params.get(f)) != str(expected.get(f)):
            mismatches.append((f, params.get(f), expected.get(f)))
    assert not mismatches, f"Mismatched fields: {mismatches}"


@pytest.mark.skipif(not TOML_FILE.exists(), reason="ProteoBench fixture missing")
def test_metamorpheus_input_order_insensitive():
    direct = extract_params(TOML_FILE, VERSION_FILE).to_series()
    reversed_ = extract_params(VERSION_FILE, TOML_FILE).to_series()
    assert direct.equals(reversed_)
