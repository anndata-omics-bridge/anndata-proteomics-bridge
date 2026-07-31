"""Built-in canonical-modification registry.

The TOML data file (``unimod_registry.toml``, sibling of this module) is
the single source of truth for ``name``, ``target``, ``position`` and
``mass_delta`` of each supported modification. Per-tool parsing-rule
Parsing rules reference modifications by accession only; the runtime resolves
them via this registry, raising an error if the accession is unknown.
"""

from __future__ import annotations

import math
import tomllib
from functools import lru_cache
from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field


class UnimodEntry(BaseModel):
    """One canonical modification record."""

    model_config = ConfigDict(extra="forbid")

    accession: str
    name: str
    aliases: list[str] = Field(default_factory=list)
    target: list[str]  # allowed residues/termini, e.g. ["S", "T", "Y"] for Phospho
    position: str
    mass_delta: float


class UnimodRegistry(BaseModel):
    """Top-level shape of the registry TOML file."""

    model_config = ConfigDict(extra="forbid")

    entries: list[UnimodEntry] = Field(min_length=1)


_REGISTRY_TOML = Path(__file__).with_name("unimod_registry.toml")


@lru_cache(maxsize=1)
def load_registry() -> dict[str, UnimodEntry]:
    """Load the bundled registry as ``{accession: UnimodEntry}``.

    Cached after the first call so re-loads in tests are free.
    """
    data = tomllib.loads(_REGISTRY_TOML.read_text(encoding="utf-8"))
    parsed = UnimodRegistry(**data)
    by_accession: dict[str, UnimodEntry] = {}
    for entry in parsed.entries:
        if entry.accession in by_accession:
            raise ValueError(f"duplicate accession in unimod_registry.toml: {entry.accession!r}")
        by_accession[entry.accession] = entry
    return by_accession


def resolve(accession: str) -> UnimodEntry:
    """Return the canonical record for ``accession`` or raise ``KeyError``."""
    registry = load_registry()
    try:
        return registry[accession]
    except KeyError:
        raise KeyError(
            f"accession {accession!r} not found in unimod_registry.toml; "
            f"add it there before referencing it from a parsing rule"
        ) from None


def find_by_name(name: str) -> UnimodEntry | None:
    """Find a canonical modification by accession, name, or shared synonym.

    Unknown names return ``None`` so parameter parsers can preserve
    vendor-specific vocabulary instead of forcing it into the small bundled
    registry.
    """
    normalized = name.strip().casefold()
    registry = load_registry()

    for accession, entry in registry.items():
        known_names = {
            accession.casefold(),
            entry.name.casefold(),
            *(alias.casefold() for alias in entry.aliases),
        }
        if normalized in known_names:
            return entry

    return None


def find_by_mass(
    mass_delta: float,
    *,
    tolerance: float = 0.001,
) -> UnimodEntry | None:
    """Find one canonical modification by monoisotopic mass within tolerance.

    Unknown masses return ``None``. An ambiguous match raises because silently
    choosing one identity would corrupt the search-parameter record.
    """
    if tolerance < 0:
        raise ValueError("mass tolerance must be non-negative")

    matches = [
        entry
        for entry in load_registry().values()
        if math.isclose(
            entry.mass_delta,
            mass_delta,
            rel_tol=0,
            abs_tol=tolerance,
        )
    ]
    if len(matches) > 1:
        accessions = ", ".join(entry.accession for entry in matches)
        raise ValueError(
            f"mass delta {mass_delta} is ambiguous within {tolerance} Da: {accessions}"
        )
    return matches[0] if matches else None
