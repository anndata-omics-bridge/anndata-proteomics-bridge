"""Load and validate self-contained software-version JSON rule documents."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from pydantic import ValidationError

from anndata_proteomics.rules._discovery import (
    RuleLocator,
    document_paths_for_software,
)
from anndata_proteomics.rules.schema import (
    ParseRule,
    ParseRuleDocument,
    QuantificationLevel,
)

JsonObject = dict[str, Any]


class RuleDocumentError(ValueError):
    """Raised when a JSON rule document or requested level is invalid."""


def _reject_json_constant(value: str) -> None:
    """Reject JavaScript constants that are not part of standard JSON."""
    raise RuleDocumentError(f"non-standard JSON value {value!r} is not allowed")


def _object_without_duplicates(pairs: list[tuple[str, Any]]) -> JsonObject:
    """Build one JSON object while rejecting ambiguous duplicate keys."""
    result: JsonObject = {}
    for key, value in pairs:
        if key in result:
            raise RuleDocumentError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def parse_rule_source(text: str, *, path: Path | str = "<memory>") -> JsonObject:
    """Parse one strict JSON object without yet applying the Pydantic schema."""
    source = Path(path)
    try:
        data = json.loads(
            text,
            object_pairs_hook=_object_without_duplicates,
            parse_constant=_reject_json_constant,
        )
    except (json.JSONDecodeError, RuleDocumentError) as exc:
        exc.add_note(f"in {source}")
        raise
    if not isinstance(data, dict):
        raise RuleDocumentError(f"rule document root must be a JSON object; in {source}")
    return data


def validate_rule_source(
    data: JsonObject,
    *,
    path: Path | str = "<memory>",
) -> ParseRuleDocument:
    """Validate a source object and every effective level through Pydantic."""
    source = Path(path)
    try:
        document = ParseRuleDocument.model_validate(data)
    except ValidationError as exc:
        exc.add_note(f"in {source}")
        raise
    for level in document.levels:
        try:
            document.effective_rule(level)
        except ValidationError as exc:
            exc.add_note(f"in {source}; level: {level}")
            raise
    return document


def parse_rule_document(
    text: str,
    *,
    path: Path | str = "<memory>",
) -> ParseRuleDocument:
    """Parse and fully validate one software-version JSON rule document."""
    return validate_rule_source(parse_rule_source(text, path=path), path=path)


def read_rule_document(path: Path | str) -> JsonObject:
    """Read one strict JSON source object without Pydantic normalization."""
    source = Path(path)
    if not source.exists():
        raise FileNotFoundError(source)
    return parse_rule_source(source.read_text(encoding="utf-8"), path=source)


def load_rule_document(path: Path | str) -> ParseRuleDocument:
    """Load one document and validate every declared level."""
    source = Path(path)
    if not source.exists():
        raise FileNotFoundError(source)
    return parse_rule_document(source.read_text(encoding="utf-8"), path=source)


def load_rule(
    source: RuleLocator | Path | str,
    level: QuantificationLevel | None = None,
) -> ParseRule:
    """Load one effective level from a document or locator.

    A path may omit ``level`` only when its document declares exactly one level.
    """
    if isinstance(source, RuleLocator):
        path = source.path
        selected_level = source.level
    else:
        path = Path(source)
        selected_level = level
    document = load_rule_document(path)
    if selected_level is None:
        if len(document.levels) != 1:
            raise RuleDocumentError(
                f"{path} contains levels {list(document.levels)}; select one explicitly"
            )
        selected_level = next(iter(document.levels))
    try:
        return document.effective_rule(selected_level)
    except KeyError as exc:
        raise RuleDocumentError(
            f"{path} has no level {selected_level!r}; available: {list(document.levels)}"
        ) from exc


def load_rules(path: Path | str) -> dict[QuantificationLevel, ParseRule]:
    """Load every effective rule from one document."""
    return load_rule_document(path).effective_rules()


def _software_version_matches(pattern: str, version: str) -> bool:
    """Return whether a document's version regex matches a parsed version."""
    try:
        return re.search(pattern, version) is not None
    except re.error as exc:
        raise ValueError(f"invalid software_version regex {pattern!r}") from exc


def _rules_equivalent_without_version(rules: list[ParseRule]) -> bool:
    """Return whether version documents expose the same conversion rule body."""
    normalized = []
    for rule in rules:
        data = rule.model_dump(by_alias=True, mode="json")
        data.pop("software_version", None)
        data.pop("file_version", None)
        normalized.append(data)
    return all(item == normalized[0] for item in normalized[1:])


def resolve_rule_locator(
    software: str,
    level: QuantificationLevel,
    version: str | None,
) -> RuleLocator | None:
    """Resolve a document-level pair using the existing software-version regexes."""
    candidates: list[tuple[Path, ParseRuleDocument]] = []
    for path in document_paths_for_software(software):
        document = load_rule_document(path)
        if level not in document.levels:
            continue
        if version is not None and not _software_version_matches(
            document.software_version, version
        ):
            continue
        candidates.append((path, document))
    if len(candidates) == 1:
        return RuleLocator(path=candidates[0][0], level=level)
    if version is None and candidates:
        rules = [document.effective_rule(level) for _, document in candidates]
        if _rules_equivalent_without_version(rules):
            return RuleLocator(path=candidates[0][0], level=level)
    return None


def load_packaged_rule(
    software: str,
    quantification_level: QuantificationLevel,
    version: str | None = None,
) -> ParseRule:
    """Load a packaged effective rule for software, level, and optional version."""
    locator = resolve_rule_locator(software, quantification_level, version)
    if locator is None:
        raise ValueError(
            f"no packaged rule for software={software!r} "
            f"level={quantification_level!r} version={version!r}"
        )
    return load_rule(locator)


def resolve_rule_for_version(
    software: str,
    quantification_level: QuantificationLevel,
    version: str | None,
) -> ParseRule | None:
    """Return the effective level covered by ``version``, or ``None``."""
    locator = resolve_rule_locator(software, quantification_level, version)
    return None if locator is None else load_rule(locator)
