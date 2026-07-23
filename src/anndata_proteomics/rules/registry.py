"""Discover packaged software-version rule documents and their level locators."""

from __future__ import annotations

from collections.abc import Iterator

from anndata_proteomics.rules._discovery import (
    RuleLocator as RuleLocator,
    document_paths_for_software as document_paths_for_software,
    document_vendor as document_vendor,
    iter_packaged_documents,
    packaged_rules_root as packaged_rules_root,
)
from anndata_proteomics.rules.schema import QuantificationLevel


class RuleNotFound(LookupError):
    """Raised when no packaged document level covers a requested software version."""


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
