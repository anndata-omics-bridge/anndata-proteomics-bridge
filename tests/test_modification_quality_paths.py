"""Edge-path coverage for modification normalization and SDRF rendering."""

from __future__ import annotations

import pandas as pd
import pytest
from pydantic import ValidationError

from anndata_proteomics.modifications import apply_rules
from anndata_proteomics.modifications.apply_rules import (
    MapEntry,
    ModificationRule,
    apply_rule,
)
from anndata_proteomics.modifications.model import (
    ModificationOccurrence,
    SearchedModification,
)
from anndata_proteomics.modifications.pipeline import apply_modifications
from anndata_proteomics.modifications.proforma import render_proforma
from anndata_proteomics.modifications.sdrf import from_sdrf_value, to_sdrf_value
from anndata_proteomics.rules.schema import (
    ModificationMapEntry,
    TokenRegexModifications,
)


def test_target_matching_and_numeric_entry_without_mass() -> None:
    assert apply_rules._target_matches(None, None, "Anywhere")
    assert not apply_rules._target_matches(["C-term"], None, "N-term")
    assert (
        apply_rules._match_entry(
            [MapEntry(token="name", name="Named")],
            "1.25",
            "M",
            "Anywhere",
            False,
        )
        is None
    )
    assert (
        apply_rules._match_entry(
            [
                MapEntry(
                    token="1.25",
                    name="Mass",
                    mass_delta=1.25,
                    position="N-term",
                )
            ],
            "1.25",
            "M",
            "Anywhere",
            False,
        )
        is None
    )
    assert (
        apply_rules._match_entry(
            [MapEntry(token="named", name="Named", position="N-term")],
            "named",
            "M",
            "Anywhere",
            False,
        )
        is None
    )


def test_before_residue_without_a_following_residue_and_unknown_positions() -> None:
    terminal = apply_rule(
        "PEP(x)",
        ModificationRule(
            source_column="Modified",
            token_pattern=r"\(([^)]+)\)",
            token_position="before_residue",
            unknown_policy="preserve",
        ),
    )
    assert terminal.stripped_sequence == "PEP"
    assert terminal.unknown_tokens == ["x"]

    nterm = apply_rule(
        "(n)PEP",
        ModificationRule(
            source_column="Modified",
            token_pattern=r"\(([^)]+)\)",
            token_position="after_residue",
            unknown_policy="preserve",
        ),
    )
    assert nterm.proforma_sequence == "[n]-PEP"

    unlocalized = apply_rule(
        "1(x)2",
        ModificationRule(
            source_column="Modified",
            token_pattern=r"\(([^)]+)\)",
            token_position="before_residue",
            unknown_policy="preserve",
        ),
    )
    assert unlocalized.stripped_sequence == ""
    assert unlocalized.unknown_tokens == ["x"]

    cterm = apply_rule(
        "PEP(x)",
        ModificationRule(
            source_column="Modified",
            token_pattern=r"\(([^)]+)\)",
            token_position="after_residue",
            entries=(MapEntry(token="x", name="C terminal", position="C-term"),),
        ),
    )
    assert cterm.proforma_sequence == "PEP-[C terminal]"

    unknown_cterm = apply_rule(
        "PEP(y)",
        ModificationRule(
            source_column="Modified",
            token_pattern=r"\(([^)]+)\)",
            token_position="after_residue",
            unknown_policy="preserve",
        ),
    )
    assert unknown_cterm.proforma_sequence == "PEP-[y]"


def test_modification_models_reject_invalid_identity_and_position() -> None:
    with pytest.raises(ValidationError, match="accession"):
        SearchedModification(name="bad", accession="invalid")
    with pytest.raises(ValidationError, match="accession"):
        ModificationOccurrence(name="bad", accession="invalid")
    with pytest.raises(ValidationError, match="non-negative"):
        ModificationOccurrence(name="bad", sequence_index=-1)


def test_apply_modifications_requires_source_column() -> None:
    settings = TokenRegexModifications(
        parser="token_regex",
        source_column="Modified",
        token_pattern=r"\(([^)]+)\)",
        map=[ModificationMapEntry(token="ox", accession="UNIMOD:35")],
    )
    with pytest.raises(KeyError, match="source_column"):
        apply_modifications(pd.DataFrame({"Other": ["PEP"]}), settings)


def test_proforma_renders_terminal_and_unlocalized_occurrences() -> None:
    rendered = render_proforma(
        "PEP",
        [
            ModificationOccurrence(name="C terminal", position="C-term"),
            ModificationOccurrence(name="ignored", position="unknown"),
        ],
        unknown_tokens={-1: "unknown-n", 3: "unknown-c"},
    )
    assert rendered == "[unknown-n]-PEP-[C terminal][unknown-c]"


def test_sdrf_skips_bare_tokens_and_rejects_unknown_modification_type() -> None:
    parsed = from_sdrf_value("ignored;NT=Oxidation")
    assert parsed.name == "Oxidation"
    assert (
        to_sdrf_value(SearchedModification(name="No position", position=None)) == "NT=No position"
    )
    with pytest.raises(ValueError, match="unknown MT"):
        from_sdrf_value("NT=Oxidation;MT=impossible")
