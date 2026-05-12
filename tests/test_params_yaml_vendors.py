"""Equivalence tests for YAML-based parameter parsers (AlphaPept, WOMBAT)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from anndata_proteomics.params.alphapept import extract_params as alphapept_extract
from anndata_proteomics.params.model import Parameters
from anndata_proteomics.params.wombat import extract_params as wombat_extract

PROTEOBENCH_PARAMS = Path("/Users/wolski/projects/anndata_bridge/ProteoBench/test/params")

COMMON_FIELDS = [
    "software_name",
    "software_version",
    "search_engine",
    "enzyme",
    "allowed_miscleavages",
    "max_mods",
    "min_peptide_length",
    "max_peptide_length",
    "min_precursor_charge",
    "max_precursor_charge",
    "precursor_mass_tolerance",
    "fragment_mass_tolerance",
    "enable_match_between_runs",
]


def _expected(csv: Path) -> Parameters:
    df = pd.read_csv(csv, header=0, index_col=0)
    return Parameters.from_series(df.iloc[:, 0])


def _compare(params: Parameters, csv: Path, extra: list[str] = ()) -> None:
    expected = _expected(csv)
    e = expected.to_series()
    a = params.to_series()
    mismatches = []
    for field in list(COMMON_FIELDS) + list(extra):
        if str(a.get(field)) != str(e.get(field)):
            mismatches.append((field, a.get(field), e.get(field)))
    assert not mismatches, f"Mismatched fields: {mismatches}"


@pytest.mark.parametrize(
    "yaml_name",
    ["alphapept_0.4.9.yaml", "alphapept_0.4.9_unnormalized.yaml"],
)
def test_alphapept_matches_proteobench(yaml_name):
    yaml_path = PROTEOBENCH_PARAMS / yaml_name
    csv_path = yaml_path.with_suffix(".csv")
    if not yaml_path.exists():
        pytest.skip("ProteoBench fixture missing")
    params = alphapept_extract(yaml_path)
    _compare(params, csv_path, extra=["ident_fdr_psm", "ident_fdr_protein", "fixed_mods", "variable_mods"])


def test_wombat_matches_proteobench():
    yaml_path = PROTEOBENCH_PARAMS / "wombat_params.yaml"
    csv_path = PROTEOBENCH_PARAMS / "wombat_params.csv"
    if not yaml_path.exists():
        pytest.skip("ProteoBench fixture missing")
    params = wombat_extract(yaml_path)
    _compare(
        params,
        csv_path,
        extra=[
            "ident_fdr_psm",
            "ident_fdr_peptide",
            "ident_fdr_protein",
            "fixed_mods",
            "variable_mods",
            "abundance_normalization_ions",
        ],
    )
