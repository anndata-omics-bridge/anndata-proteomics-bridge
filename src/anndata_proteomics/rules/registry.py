"""Locate packaged parsing-rule TOMLs by (software, quantification_level, file_version)."""

from __future__ import annotations

from collections.abc import Iterator
from importlib import resources
from pathlib import Path


class RuleNotFound(LookupError):
    """Raised when (software, quantification_level, file_version) has no packaged TOML."""


def packaged_rules_root() -> Path:
    """Filesystem path to the parsing_rules/ directory inside the package."""
    traversable = resources.files("anndata_proteomics") / "parsing_rules"
    return Path(str(traversable))


def iter_packaged_rules() -> Iterator[Path]:
    """Yield every packaged parse_*.toml under parsing_rules/<vendor>/, sorted."""
    yield from sorted(packaged_rules_root().glob("*/parse_*.toml"))


def find_rule(
    software: str, quantification_level: str, file_version: str = "1"
) -> Path:
    """Resolve (software, level, version) to a packaged TOML path.

    The vendor directory and the filename both use the lowercase software token
    (e.g. "diann", "fragpipe") — not the human-readable software_name from the TOML.
    """
    root = packaged_rules_root()
    vendor_dir = root / software
    candidate = (
        vendor_dir / f"parse_{software}_{quantification_level}_{file_version}.toml"
    )
    if candidate.exists():
        return candidate
    available = (
        sorted(p.name for p in vendor_dir.glob("parse_*.toml"))
        if vendor_dir.exists()
        else None
    )
    if available is None:
        raise RuleNotFound(
            f"No packaged rules for software={software!r}; "
            f"vendor folder {vendor_dir} does not exist."
        )
    raise RuleNotFound(
        f"No packaged rule at {candidate.name} in {software}/; "
        f"available in {software}/: {available}"
    )
