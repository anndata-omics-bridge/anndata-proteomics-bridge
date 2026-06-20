"""MaxQuant XML parser equivalence tests."""

from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd
import pytest

from anndata_proteomics.params.maxquant import extract_params
from anndata_proteomics.params.model import Parameters

PROTEOBENCH_PARAMS = Path(__file__).resolve().parent / "params"

CASES = [
    ("mqpar1.5.3.30_MBR.xml", "mqpar1.5.3.30_MBR_sel.json"),
    ("mqpar_MQ1.6.3.3_MBR.xml", "mqpar_MQ1.6.3.3_MBR_sel.json"),
    ("mqpar_MQ2.1.3.0_noMBR.xml", "mqpar_MQ2.1.3.0_noMBR_sel.json"),
    ("mqpar_mq2.6.2.0_1mc_MBR.xml", "mqpar_mq2.6.2.0_1mc_MBR_sel.json"),
]


def _normalize(v):
    """Treat NaN, None, and 'NaN' as equivalent for comparison."""
    if v is None:
        return None
    if isinstance(v, float) and math.isnan(v):
        return None
    return v


@pytest.mark.parametrize(("xml_name", "expected_name"), CASES)
def test_maxquant_matches_proteobench(xml_name, expected_name):
    xml_path = PROTEOBENCH_PARAMS / xml_name
    expected_path = PROTEOBENCH_PARAMS / expected_name
    if not xml_path.exists() or not expected_path.exists():
        pytest.skip("ProteoBench fixture missing")
    # ProteoBench wrote literal "NaN" tokens to JSON; allow them.
    expected_raw = json.loads(expected_path.read_text().replace("NaN", "null"))
    # Round-trip the expected values through the shared model so tolerance and
    # modification fields are serialized identically on both sides (the model
    # canonicalizes mass-tolerance dicts and the modification separator). This
    # mirrors the DIA-NN / Spectronaut equivalence tests; it does not weaken the
    # comparison, only removes formatting differences the model owns.
    expected = Parameters.from_series(pd.Series(expected_raw)).to_series()
    params = extract_params(xml_path).to_series()

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
        "semi_enzymatic",
        "allowed_miscleavages",
        "min_peptide_length",
        "fixed_mods",
        "variable_mods",
        "max_mods",
        "max_precursor_charge",
    ]
    mismatches = []
    for field in fields:
        actual = _normalize(params.get(field))
        expected_value = _normalize(expected.get(field))
        if actual != expected_value:
            mismatches.append((field, actual, expected_value))
    assert not mismatches, f"Mismatched fields: {mismatches}"
