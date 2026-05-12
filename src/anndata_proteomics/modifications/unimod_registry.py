"""Built-in canonical-modification registry.

The TOML data file (``unimod_registry.toml``, sibling of this module) is
the single source of truth for ``name``, ``target``, ``position`` and
``mass_delta`` of each supported modification. Per-tool parsing-rule
TOMLs reference modifications by accession only; the runtime resolves
them via this registry, raising an error if the accession is unknown.
"""

from __future__ import annotations

import tomllib
from functools import lru_cache
from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field


class UnimodEntry(BaseModel):
    """One canonical modification record."""

    model_config = ConfigDict(extra="forbid")

    accession: str
    name: str
    target: str
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
            raise ValueError(
                f"duplicate accession in unimod_registry.toml: {entry.accession!r}"
            )
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
            f"add it there before referencing it from a parsing-rule TOML"
        ) from None
