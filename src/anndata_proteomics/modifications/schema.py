"""Parsing-rule models for vendor modification encodings."""

from __future__ import annotations

from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field

TokenPosition = Literal[
    "before_residue", "after_residue", "n_term", "c_term", "embedded", "unknown"
]
UnknownPolicy = Literal["preserve", "drop", "error"]


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid")


class ModificationMapEntry(_Strict):
    """User-facing JSON entry: a vendor token plus the Unimod accession.

    ``name``, ``target``, ``position`` and ``mass_delta`` are NOT carried
    on the entry itself — they are filled at rule-load time from
    ``modifications/unimod_registry.toml``. This keeps the per-tool rule files
    free of duplicated canonical data and guarantees that all tools agree
    on what e.g. ``UNIMOD:35`` means.
    """

    token: str
    accession: str


class TokenRegexModifications(_Strict):
    """Extract vendor modification tokens with a regex and map them to Unimod.

    For vendors that write modifications inline in the sequence itself, e.g.
    ``"PEPM[15.9949]TIDE"`` or ``"_(ac)PEPTIDEM(ox)_"``.
    """

    parser: Literal["token_regex"]
    source_column: str
    token_pattern: str
    token_position: TokenPosition = "after_residue"
    case_sensitive: bool = False
    unknown_policy: UnknownPolicy = "preserve"
    output_column: str = "proforma_sequence"
    map: list[ModificationMapEntry] = Field(min_length=1)

    @property
    def source_columns(self) -> tuple[str, ...]:
        """Input columns this parser reads."""
        return (self.source_column,)


class SiteListModifications(_Strict):
    """Read modifications from parallel name and site columns beside a bare sequence.

    The alphabase layout, used by AlphaDIA::

        sequence     mods                            mod_sites
        TCSSFIAAMER  Oxidation@M;Carbamidomethyl@C   9;2

    ``mods`` and ``mod_sites`` are ``delimiter``-joined and paired **index-wise**,
    not sorted — above, ``Oxidation@M`` sits at 9 and ``Carbamidomethyl@C`` at 2.
    Sites are ``site_base``-indexed into the sequence; site ``0`` denotes the
    N-terminus regardless of ``site_base`` (AlphaDIA writes ``Acetyl@Protein_N-term``
    that way). Tokens are matched exactly against ``map``, since a name like
    ``Oxidation@M`` already carries its own target.
    """

    parser: Literal["site_list"]
    sequence_column: str
    modification_column: str
    site_column: str
    delimiter: str = ";"
    site_base: int = Field(default=1, ge=0, le=1)
    case_sensitive: bool = False
    unknown_policy: UnknownPolicy = "preserve"
    output_column: str = "proforma_sequence"
    map: list[ModificationMapEntry] = Field(min_length=1)

    @property
    def source_columns(self) -> tuple[str, ...]:
        """Input columns this parser reads."""
        return (self.sequence_column, self.modification_column, self.site_column)


Modifications = Annotated[
    TokenRegexModifications | SiteListModifications,
    Field(discriminator="parser"),
]
