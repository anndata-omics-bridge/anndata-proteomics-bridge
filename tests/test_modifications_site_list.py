"""The ``site_list`` modification parser (alphabase-style parallel columns)."""

from __future__ import annotations

import pandas as pd
import pytest
from pydantic import ValidationError

from anndata_proteomics.modifications.pipeline import apply_modifications
from anndata_proteomics.rules.schema import SiteListModifications

_MAP = [
    {"token": "Carbamidomethyl@C", "accession": "UNIMOD:4"},
    {"token": "Oxidation@M", "accession": "UNIMOD:35"},
    {"token": "Acetyl@Protein_N-term", "accession": "UNIMOD:1"},
]


def _rule(**overrides: object) -> SiteListModifications:
    return SiteListModifications.model_validate(
        {
            "parser": "site_list",
            "sequence_column": "sequence",
            "modification_column": "mods",
            "site_column": "mod_sites",
            "map": _MAP,
        }
        | overrides
    )


def _apply(sequence: list[str], mods: list[str | None], sites: list[str | None]) -> pd.DataFrame:
    frame = pd.DataFrame({"sequence": sequence, "mods": mods, "mod_sites": sites})
    return apply_modifications(frame, _rule())


def test_sites_pair_index_wise_not_in_sorted_order() -> None:
    """AlphaDIA writes ``Oxidation@M;Carbamidomethyl@C`` with sites ``9;2``.

    Pairing by sorted site instead of by index would swap the two accessions.
    """
    out = _apply(["TCSSFIAAMER"], ["Oxidation@M;Carbamidomethyl@C"], ["9;2"])

    assert out.loc[0, "proforma_sequence"] == "TC[UNIMOD:4]SSFIAAM[UNIMOD:35]ER"
    assert out.loc[0, "stripped_sequence"] == "TCSSFIAAMER"


def test_site_zero_is_the_n_terminus() -> None:
    out = _apply(["MDESSGGK"], ["Acetyl@Protein_N-term"], ["0"])

    assert out.loc[0, "proforma_sequence"] == "[UNIMOD:1]-MDESSGGK"


def test_unmodified_rows_render_the_bare_sequence() -> None:
    out = _apply(["PEPTIDE", "PEPTIDE"], [None, ""], [None, ""])

    assert out["proforma_sequence"].tolist() == ["PEPTIDE", "PEPTIDE"]
    assert out["unknown_mod_tokens"].tolist() == [[], []]


def test_same_sequence_and_charge_keep_distinct_peptidoforms() -> None:
    """The reason this parser is required rather than optional.

    Without it both rows render ``TCSSFIAAMER`` and collapse into one feature,
    silently summing an oxidised and a non-oxidised precursor.
    """
    out = _apply(
        ["TCSSFIAAMER", "TCSSFIAAMER"],
        ["Oxidation@M;Carbamidomethyl@C", "Carbamidomethyl@C"],
        ["9;2", "2"],
    )

    assert out["proforma_sequence"].nunique() == 2


def test_repeated_token_localizes_to_each_site() -> None:
    out = _apply(["CAC"], ["Carbamidomethyl@C;Carbamidomethyl@C"], ["1;3"])

    assert out.loc[0, "proforma_sequence"] == "C[UNIMOD:4]AC[UNIMOD:4]"


def test_unknown_token_is_preserved_verbatim_by_default() -> None:
    out = _apply(["PEPTIDE"], ["Phospho@S"], ["3"])

    assert out.loc[0, "proforma_sequence"] == "PEP[Phospho@S]TIDE"
    assert out.loc[0, "unknown_mod_tokens"] == ["Phospho@S"]


def test_unknown_policy_drop_omits_the_token() -> None:
    frame = pd.DataFrame({"sequence": ["PEPTIDE"], "mods": ["Phospho@S"], "mod_sites": ["3"]})
    out = apply_modifications(frame, _rule(unknown_policy="drop"))

    assert out.loc[0, "proforma_sequence"] == "PEPTIDE"


def test_mismatched_token_and_site_counts_raise() -> None:
    with pytest.raises(ValueError, match="length mismatch"):
        _apply(["PEPTIDE"], ["Oxidation@M;Carbamidomethyl@C"], ["3"])


def test_non_integer_site_raises() -> None:
    with pytest.raises(ValueError, match="non-integer modification site"):
        _apply(["PEPTIDE"], ["Oxidation@M"], ["middle"])


def test_missing_source_column_names_the_parser() -> None:
    frame = pd.DataFrame({"sequence": ["PEPTIDE"], "mods": ["Oxidation@M"]})
    with pytest.raises(KeyError, match="needs column"):
        apply_modifications(frame, _rule())


def test_site_base_zero_shifts_residue_indexing() -> None:
    frame = pd.DataFrame({"sequence": ["CAC"], "mods": ["Carbamidomethyl@C"], "mod_sites": ["2"]})
    out = apply_modifications(frame, _rule(site_base=0))

    assert out.loc[0, "proforma_sequence"] == "CAC[UNIMOD:4]"


def test_parser_discriminator_rejects_token_regex_fields() -> None:
    with pytest.raises(ValidationError):
        SiteListModifications.model_validate(
            {
                "parser": "site_list",
                "sequence_column": "sequence",
                "modification_column": "mods",
                "site_column": "mod_sites",
                "token_pattern": r"\((.*?)\)",
                "map": _MAP,
            }
        )
