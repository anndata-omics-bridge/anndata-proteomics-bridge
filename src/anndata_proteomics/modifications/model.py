"""Modification identity models."""

from __future__ import annotations

from enum import Enum
from typing import Any

from pydantic import BaseModel, ConfigDict, Field


class ModType(str, Enum):
    """Whether a modification was searched as fixed, variable, or unknown."""

    fixed = "fixed"
    variable = "variable"
    unknown = "unknown"


class SearchedModification(BaseModel):
    """A modification declared in a search-engine parameter file.

    Used for SDRF metadata export (``comment[modification parameters]``).
    Carries no sequence localization — that lives on
    :class:`ModificationOccurrence`.
    """

    model_config = ConfigDict(extra="forbid")

    name: str
    accession: str | None = None
    mod_type: ModType = ModType.unknown
    target: str | None = None
    position: str | None = "Anywhere"
    mass_delta: float | None = None
    source: str | None = None


class ModificationOccurrence(BaseModel):
    """A localized modification on a peptide.

    Used to build ProForma sequence strings from vendor modified-sequence
    columns.
    """

    model_config = ConfigDict(extra="forbid")

    name: str
    accession: str | None = None
    target_residue: str | None = None
    sequence_index: int | None = None
    position: str | None = None
    mass_delta: float | None = None
    source_token: str | None = None


class ModifiedSequence(BaseModel):
    """A modified peptide as observed in a quantification result row."""

    model_config = ConfigDict(extra="forbid")

    stripped_sequence: str
    proforma_sequence: str
    occurrences: list[ModificationOccurrence] = Field(default_factory=list)
    source_sequence: str | None = None
    unknown_tokens: list[str] = Field(default_factory=list)
