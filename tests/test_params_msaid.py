"""MSAID parameter parser equivalence tests."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from anndata_proteomics.params.model import Parameters
from anndata_proteomics.params.msaid import extract_params

PROTEOBENCH_PARAMS = Path("/Users/wolski/projects/anndata_bridge/ProteoBench/test/params")


def _expected_tsv(path: Path) -> Parameters:
    df = pd.read_csv(path, sep="\t", header=0, index_col=0)
    return Parameters.from_series(df.iloc[:, 0])


def test_msaid_matches_proteobench():
    csv = PROTEOBENCH_PARAMS / "MSAID_default_params.csv"
    expected_tsv = PROTEOBENCH_PARAMS / "MSAID_default_params.tsv"
    if not csv.exists() or not expected_tsv.exists():
        pytest.skip("ProteoBench fixture missing")

    params = extract_params(csv).model_dump()
    expected = _expected_tsv(expected_tsv).model_dump()
    # semi_enzymatic intentionally excluded: ProteoBench's dynamic dataclass only
    # populates fields declared in the JSON template, and DIA_ion.json omits
    # semi_enzymatic, so the expected TSV has it blank. APB always emits it.
    fields = [
        "software_name",
        "search_engine",
        "search_engine_version",
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
        "quantification_method",
        "enable_match_between_runs",
    ]
    mismatches = []
    for f in fields:
        if str(params.get(f)) != str(expected.get(f)):
            mismatches.append((f, params.get(f), expected.get(f)))
    assert not mismatches, f"Mismatched fields: {mismatches}"
