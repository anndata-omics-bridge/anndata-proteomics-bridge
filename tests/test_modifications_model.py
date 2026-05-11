"""Tests for modification models, SDRF, ProForma, and token application."""

from __future__ import annotations

import pytest

from anndata_proteomics.modifications.apply_rules import (
    MapEntry,
    ModificationRule,
    apply_rule,
)
from anndata_proteomics.modifications.model import (
    ModificationOccurrence,
    ModType,
    SearchedModification,
)
from anndata_proteomics.modifications.proforma import render_proforma
from anndata_proteomics.modifications.sdrf import from_sdrf_value, to_sdrf_value


# --- SDRF -------------------------------------------------------------------


def test_sdrf_render_canonical_order():
    mod = SearchedModification(
        name="Oxidation",
        accession="UNIMOD:35",
        mod_type=ModType.variable,
        target="M",
        position="Anywhere",
    )
    assert to_sdrf_value(mod) == "NT=Oxidation;AC=UNIMOD:35;MT=variable;TA=M;PP=Anywhere"


def test_sdrf_render_omits_empty_fields():
    mod = SearchedModification(name="Unknown")
    assert to_sdrf_value(mod) == "NT=Unknown;PP=Anywhere"


def test_sdrf_parse_order_insensitive():
    parsed = from_sdrf_value("PP=Anywhere;MT=variable;NT=Oxidation;AC=UNIMOD:35;TA=M")
    assert parsed.name == "Oxidation"
    assert parsed.accession == "UNIMOD:35"
    assert parsed.mod_type is ModType.variable
    assert parsed.target == "M"
    assert parsed.position == "Anywhere"


def test_sdrf_roundtrip():
    canonical = "NT=Carbamidomethyl;AC=UNIMOD:4;MT=fixed;TA=C;PP=Anywhere"
    assert to_sdrf_value(from_sdrf_value(canonical)) == canonical


def test_sdrf_missing_nt_raises():
    with pytest.raises(ValueError):
        from_sdrf_value("AC=UNIMOD:35;TA=M")


# --- ProForma rendering ----------------------------------------------------


def test_proforma_residue_with_accession():
    occ = [ModificationOccurrence(name="Oxidation", accession="UNIMOD:35", sequence_index=0, position="Anywhere")]
    assert render_proforma("MPEPTIDE", occ) == "M[UNIMOD:35]PEPTIDE"


def test_proforma_nterm_and_internal():
    occ = [
        ModificationOccurrence(name="Acetyl", accession="UNIMOD:1", position="N-term"),
        ModificationOccurrence(name="Oxidation", accession="UNIMOD:35", sequence_index=3, position="Anywhere"),
    ]
    assert render_proforma("PEPMIDE", occ) == "[UNIMOD:1]-PEPM[UNIMOD:35]IDE"


def test_proforma_falls_back_to_name():
    occ = [ModificationOccurrence(name="Oxidation", sequence_index=0)]
    assert render_proforma("M", occ) == "M[Oxidation]"


def test_proforma_unknown_token_preserved():
    s = render_proforma("MPEPTIDE", [], unknown_tokens={0: "weirdmass"})
    assert s == "M[weirdmass]PEPTIDE"


# --- apply_rule -------------------------------------------------------------


_OX_M = MapEntry(token="15.9949", name="Oxidation", accession="UNIMOD:35", target="M", position="Anywhere", mass_delta=15.9949)
_CAM_C = MapEntry(token="57.0215", name="Carbamidomethyl", accession="UNIMOD:4", target="C", position="Anywhere", mass_delta=57.02146)
_AC_NT = MapEntry(token="42.0106", name="Acetyl", accession="UNIMOD:1", target="N-term", position="N-term", mass_delta=42.0106)


def _rule(**overrides) -> ModificationRule:
    base = dict(
        source_column="Modified Sequence",
        token_pattern=r"\[([^\]]+)\]",
        token_position="after_residue",
        entries=(_OX_M, _CAM_C, _AC_NT),
    )
    base.update(overrides)
    return ModificationRule(**base)


def test_apply_rule_fragpipe_style_numeric_token():
    result = apply_rule("PEPM[15.9949]TIDE", _rule())
    assert result.stripped_sequence == "PEPMTIDE"
    assert result.proforma_sequence == "PEPM[UNIMOD:35]TIDE"
    assert len(result.occurrences) == 1
    assert result.occurrences[0].sequence_index == 3
    assert result.occurrences[0].target_residue == "M"


def test_apply_rule_nterm_acetyl():
    result = apply_rule("[42.0106]PEPTIDE", _rule())
    assert result.stripped_sequence == "PEPTIDE"
    assert result.proforma_sequence == "[UNIMOD:1]-PEPTIDE"
    assert result.occurrences[0].position == "N-term"


def test_apply_rule_maxquant_style_parens():
    rule = _rule(
        token_pattern=r"\(([^)]+)\)",
        entries=(
            MapEntry(token="ox", name="Oxidation", accession="UNIMOD:35", target="M", position="Anywhere"),
            MapEntry(token="ac", name="Acetyl", accession="UNIMOD:1", target="N-term", position="N-term"),
        ),
    )
    result = apply_rule("_(ac)PEPTM(ox)IDE_", rule)
    assert result.stripped_sequence == "PEPTMIDE"
    assert result.proforma_sequence == "[UNIMOD:1]-PEPTM[UNIMOD:35]IDE"
    assert {occ.position for occ in result.occurrences} == {"N-term", "Anywhere"}


def test_apply_rule_unknown_preserve():
    result = apply_rule("PEPM[99.9999]TIDE", _rule())
    assert "99.9999" in result.unknown_tokens
    assert result.proforma_sequence == "PEPM[99.9999]TIDE"


def test_apply_rule_unknown_error():
    rule = _rule(unknown_policy="error")
    with pytest.raises(ValueError):
        apply_rule("PEPM[99.9999]TIDE", rule)


def test_apply_rule_unknown_drop():
    rule = _rule(unknown_policy="drop")
    result = apply_rule("PEPM[99.9999]TIDE", rule)
    assert result.unknown_tokens == []
    assert result.proforma_sequence == "PEPMTIDE"


def test_apply_rule_mass_disambiguated_by_target():
    # Two entries share mass; target residue picks the right one.
    rule = _rule(
        entries=(
            MapEntry(token="79.9663", name="Phospho-S", accession="UNIMOD:21", target="S", position="Anywhere", mass_delta=79.9663),
            MapEntry(token="79.9663", name="Sulfo-Y", accession="UNIMOD:40", target="Y", position="Anywhere", mass_delta=79.9663),
        ),
    )
    on_s = apply_rule("PEPS[79.9663]TIDE", rule)
    on_y = apply_rule("PEPY[79.9663]TIDE", rule)
    assert on_s.occurrences[0].accession == "UNIMOD:21"
    assert on_y.occurrences[0].accession == "UNIMOD:40"


def test_apply_rule_alphapept_before_residue_lowercase():
    rule = _rule(
        token_pattern=r"([a-z]+)",
        token_position="before_residue",
        entries=(
            MapEntry(token="ox", name="Oxidation", accession="UNIMOD:35", target="M", position="Anywhere"),
        ),
    )
    result = apply_rule("PEPToxMIDE", rule)
    assert result.stripped_sequence == "PEPTMIDE"
    assert result.proforma_sequence == "PEPTM[UNIMOD:35]IDE"
    assert result.occurrences[0].sequence_index == 4
