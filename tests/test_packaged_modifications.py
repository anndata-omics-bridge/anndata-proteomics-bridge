"""End-to-end tests for packaged document-level modification mappings."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from anndata_proteomics.modifications.pipeline import apply_modifications
from anndata_proteomics.rules.loader import load_rule_from_path
from anndata_proteomics.rules.schema import ParseRule, QuantificationLevel

PARSING_RULES = Path(__file__).parent.parent / "src" / "anndata_proteomics" / "parsing_rules"


def _load(relative: str, level: QuantificationLevel) -> ParseRule:
    return load_rule_from_path(PARSING_RULES / relative, level)


@pytest.mark.parametrize(
    ("rule_path", "level", "modified_sequence", "expected_proforma"),
    [
        ("diann/v1/rules.json", "ion", "(UniMod:1)AAPEPTIDE", "[UNIMOD:1]-AAPEPTIDE"),
        ("diann/v1/rules.json", "ion", "PEPM(UniMod:35)TIDE", "PEPM[UNIMOD:35]TIDE"),
        ("fragpipe/rules.json", "ion", "PEPM[15.9949]TIDE", "PEPM[UNIMOD:35]TIDE"),
        ("fragpipe/rules.json", "ion", "PEPC[57.0215]TIDE", "PEPC[UNIMOD:4]TIDE"),
        (
            "maxquant/rules.json",
            "ion",
            "_(ac)PEPTM(ox)IDE_",
            "[UNIMOD:1]-PEPTM[UNIMOD:35]IDE",
        ),
        (
            "maxquant/rules.json",
            "ion",
            "_(Acetyl (Protein N-term))PEPTM(Oxidation (M))IDE_",
            "[UNIMOD:1]-PEPTM[UNIMOD:35]IDE",
        ),
        ("peaks/rules.json", "ion", "PEPM(+15.99)TIDE", "PEPM[UNIMOD:35]TIDE"),
        (
            "spectronaut/rules.json",
            "ion",
            "_[Acetyl (Protein N-term)]PEPTM[Oxidation (M)]IDE_",
            "[UNIMOD:1]-PEPTM[UNIMOD:35]IDE",
        ),
        (
            "diann/v1/rules.json",
            "fragment",
            "PEPM(UniMod:35)TIDE",
            "PEPM[UNIMOD:35]TIDE",
        ),
        (
            "spectronaut/rules.json",
            "fragment",
            "PEPM[Oxidation (M)]TIDE",
            "PEPM[UNIMOD:35]TIDE",
        ),
        (
            "wombat/rules.json",
            "ion",
            "[Acetyl]-PEPTM[Oxidation]IDE",
            "[UNIMOD:1]-PEPTM[UNIMOD:35]IDE",
        ),
    ],
)
def test_packaged_rule_token_mapping(
    rule_path: str,
    level: QuantificationLevel,
    modified_sequence: str,
    expected_proforma: str,
) -> None:
    rule = _load(rule_path, level)
    assert rule.modifications is not None
    # These cases all feed an inline modified sequence; site_list rules take parallel
    # name/site columns instead and are covered by tests/test_modifications_site_list.py.
    assert rule.modifications.parser == "token_regex"
    frame = pd.DataFrame({rule.modifications.source_column: [modified_sequence]})
    result = apply_modifications(frame, rule.modifications)
    assert result.loc[0, "proforma_sequence"] == expected_proforma
