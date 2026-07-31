"""Discover packaged software-version rule documents and their level locators."""

from __future__ import annotations

from collections.abc import Iterator

from anndata_proteomics.rules._discovery import (
    RuleLocator as RuleLocator,
)
from anndata_proteomics.rules._discovery import (
    document_paths_for_software as document_paths_for_software,
)
from anndata_proteomics.rules._discovery import (
    document_vendor as document_vendor,
)
from anndata_proteomics.rules._discovery import (
    iter_packaged_documents as iter_packaged_documents,
)
from anndata_proteomics.rules._discovery import (
    packaged_rules_root as packaged_rules_root,
)
from anndata_proteomics.rules.loader import (
    AmbiguousRuleLocators,
    RuleLocatorResolution,
    RuleLocatorUnavailable,
    load_rule_document,
    resolve_rule_locator_for_version,
    resolve_rule_locator_without_version,
)
from anndata_proteomics.rules.schema import QuantificationLevel


class RuleNotFound(LookupError):
    """Raised when no packaged document level covers a requested software version."""


class AmbiguousRuleError(LookupError):
    """Raised when several packaged locators satisfy lookup evidence."""


def iter_packaged_rules() -> Iterator[RuleLocator]:
    """Yield every packaged document-level pair in stable order."""
    for path in iter_packaged_documents():
        document = load_rule_document(path)
        for level in document.levels:
            yield RuleLocator(path=path, level=level)


def find_rule(
    software: str,
    level: QuantificationLevel,
) -> RuleLocator:
    """Require one unambiguous packaged locator without version evidence."""
    return _require_locator(resolve_rule_locator_without_version(software, level))


def find_rule_for_version(
    software: str,
    level: QuantificationLevel,
    version: str,
) -> RuleLocator:
    """Require one packaged locator for a concrete software version."""
    return _require_locator(resolve_rule_locator_for_version(software, level, version))


def _require_locator(resolution: RuleLocatorResolution) -> RuleLocator:
    if isinstance(resolution, RuleLocatorUnavailable):
        raise RuleNotFound(resolution.reason)
    if isinstance(resolution, AmbiguousRuleLocators):
        raise AmbiguousRuleError(resolution.reason)
    return resolution
