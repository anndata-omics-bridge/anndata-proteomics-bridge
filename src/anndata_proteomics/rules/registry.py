"""Discover packaged software-version rule documents and their level locators."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from importlib import resources
from pathlib import Path

from anndata_proteomics.rules.schema import QuantificationLevel


class RuleNotFound(LookupError):
    """Raised when no packaged document level covers a requested software version."""


@dataclass(frozen=True, order=True, slots=True)
class RuleLocator:
    """Operational address of one level inside a software-version document."""

    path: Path
    level: QuantificationLevel


def packaged_rules_root() -> Path:
    """Filesystem path to the package's parsing-rule documents."""
    traversable = resources.files("anndata_proteomics") / "parsing_rules"
    return Path(str(traversable))


def iter_packaged_documents() -> Iterator[Path]:
    """Yield one self-contained JSON document per software-version grouping."""
    root = packaged_rules_root()
    paths = set(root.glob("*/rules.json")) | set(root.glob("*/v*/rules.json"))
    yield from sorted(paths)


def document_vendor(path: Path | str) -> str:
    """Return the vendor folder slug for a packaged document path."""
    resolved = Path(path).resolve()
    root = packaged_rules_root().resolve()
    try:
        return resolved.relative_to(root).parts[0]
    except (ValueError, IndexError) as exc:
        raise ValueError(f"not a packaged rule document: {resolved}") from exc


def document_paths_for_software(software: str) -> tuple[Path, ...]:
    """Return every packaged version document for a vendor slug."""
    vendor = packaged_rules_root() / software
    if not vendor.is_dir():
        return ()
    paths = set(vendor.glob("rules.json")) | set(vendor.glob("v*/rules.json"))
    return tuple(sorted(paths))


def iter_packaged_rules() -> Iterator[RuleLocator]:
    """Yield every packaged document-level pair in stable order."""
    from anndata_proteomics.rules.loader import load_rule_document

    for path in iter_packaged_documents():
        document = load_rule_document(path)
        for level in document.levels:
            yield RuleLocator(path=path, level=level)


def find_rule(
    software: str,
    level: QuantificationLevel,
    version: str | None = None,
) -> RuleLocator:
    """Resolve a packaged level locator or raise :class:`RuleNotFound`."""
    from anndata_proteomics.rules.loader import resolve_rule_locator

    locator = resolve_rule_locator(software, level, version)
    if locator is None:
        raise RuleNotFound(
            f"no packaged rule for software={software!r} level={level!r} version={version!r}"
        )
    return locator
