"""End-to-end: every packaged rule with a [modifications] block parses tokens correctly."""

from __future__ import annotations

import tomllib
from pathlib import Path

import pandas as pd
import pytest

from anndata_proteomics.modifications.pipeline import apply_modifications
from anndata_proteomics.rules.schema import ParseRule

PARSING_RULES = (
    Path(__file__).parent.parent / "src" / "anndata_proteomics" / "parsing_rules"
)


def _load(rel: str) -> ParseRule:
    return ParseRule(**tomllib.loads((PARSING_RULES / rel).read_text()))


@pytest.mark.parametrize(
    ("rule_path", "modified_sequence", "expected_proforma"),
    [
        (
            "diann/parse_diann_ion.toml",
            "(UniMod:1)AAPEPTIDE",
            "[UNIMOD:1]-AAPEPTIDE",
        ),
        (
            "diann/parse_diann_ion.toml",
            "PEPM(UniMod:35)TIDE",
            "PEPM[UNIMOD:35]TIDE",
        ),
        (
            "fragpipe/parse_fragpipe_ion_1.toml",
            "PEPM[15.9949]TIDE",
            "PEPM[UNIMOD:35]TIDE",
        ),
        (
            "fragpipe/parse_fragpipe_ion_1.toml",
            "PEPC[57.0215]TIDE",
            "PEPC[UNIMOD:4]TIDE",
        ),
        (
            "maxquant/parse_maxquant_ion_1.toml",
            "_(ac)PEPTM(ox)IDE_",
            "[UNIMOD:1]-PEPTM[UNIMOD:35]IDE",
        ),
        (
            "maxquant/parse_maxquant_ion_1.toml",
            "_(Acetyl (Protein N-term))PEPTM(Oxidation (M))IDE_",
            "[UNIMOD:1]-PEPTM[UNIMOD:35]IDE",
        ),
        (
            "peaks/parse_peaks_ion_1.toml",
            "PEPM(+15.99)TIDE",
            "PEPM[UNIMOD:35]TIDE",
        ),
        (
            "spectronaut/parse_spectronaut_ion_1.toml",
            "_[Acetyl (Protein N-term)]PEPTM[Oxidation (M)]IDE_",
            "[UNIMOD:1]-PEPTM[UNIMOD:35]IDE",
        ),
        (
            "wombat/parse_wombat_peptidoform_1.toml",
            "[Acetyl]-PEPTM[Oxidation]IDE",
            "[UNIMOD:1]-PEPTM[UNIMOD:35]IDE",
        ),
    ],
)
def test_packaged_rule_token_mapping(rule_path, modified_sequence, expected_proforma):
    rule = _load(rule_path)
    assert rule.modifications is not None, f"{rule_path} lacks a [modifications] block"
    df = pd.DataFrame({rule.modifications.source_column: [modified_sequence]})
    result = apply_modifications(df, rule.modifications)
    assert result.loc[0, "proforma_sequence"] == expected_proforma
