"""Load and validate self-contained software-version JSON rule documents."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Never, cast

from pydantic import ValidationError

from anndata_proteomics.params.model import Parameters
from anndata_proteomics.rules._discovery import (
    RuleLocator,
    document_paths_for_software,
)
from anndata_proteomics.rules.schema import (
    ParseRule,
    ParseRuleDocument,
    QuantificationLevel,
    RuleCompositionError,
)

type JsonValue = str | int | float | bool | None | list[JsonValue] | dict[str, JsonValue]
type JsonObject = dict[str, JsonValue]


class RuleDocumentError(ValueError):
    """Raised when a JSON rule document or requested level is invalid."""


@dataclass(frozen=True, slots=True)
class RuleLocatorUnavailable:
    """No unique packaged document locator satisfies the requested evidence."""

    reason: str


@dataclass(frozen=True, slots=True)
class AmbiguousRuleLocators:
    """Several packaged document locators satisfy the requested evidence."""

    reason: str
    locators: tuple[RuleLocator, ...]


type RuleLocatorResolution = RuleLocator | RuleLocatorUnavailable | AmbiguousRuleLocators


@dataclass(frozen=True, slots=True)
class PackagedRuleUnavailable:
    """No unique packaged effective rule satisfies the requested evidence."""

    reason: str


@dataclass(frozen=True, slots=True)
class AmbiguousPackagedRules:
    """Several packaged effective rules satisfy the requested evidence."""

    reason: str


type PackagedRuleResolution = ParseRule | PackagedRuleUnavailable | AmbiguousPackagedRules


def _reject_json_constant(value: str) -> Never:
    """Reject JavaScript constants that are not part of standard JSON."""
    raise RuleDocumentError(f"non-standard JSON value {value!r} is not allowed")


def _object_without_duplicates(pairs: list[tuple[str, JsonValue]]) -> JsonObject:
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
        data = cast(
            JsonValue,
            json.loads(
                text,
                object_pairs_hook=_object_without_duplicates,
                parse_constant=_reject_json_constant,
            ),
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
            document.validate_effective_rule_variants(level)
        except (ValidationError, RuleCompositionError) as exc:
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


def load_rule(locator: RuleLocator) -> ParseRule:
    """Load one unparameterized effective rule from an explicit locator."""
    document = load_rule_document(locator.path)
    _require_document_level(document, locator.path, locator.level)
    return document.effective_rule(locator.level)


def load_parameterized_rule(locator: RuleLocator, parameters: Parameters) -> ParseRule:
    """Load one parameterized effective rule from an explicit locator."""
    document = load_rule_document(locator.path)
    _require_document_level(document, locator.path, locator.level)
    return document.parameterized_effective_rule(locator.level, parameters)


def load_rule_from_path(
    path: Path | str,
    level: QuantificationLevel,
) -> ParseRule:
    """Load one unparameterized level from an explicit document path."""
    source = Path(path)
    document = load_rule_document(source)
    _require_document_level(document, source, level)
    return document.effective_rule(level)


def load_parameterized_rule_from_path(
    path: Path | str,
    level: QuantificationLevel,
    parameters: Parameters,
) -> ParseRule:
    """Load one parameterized level from an explicit document path."""
    source = Path(path)
    document = load_rule_document(source)
    _require_document_level(document, source, level)
    return document.parameterized_effective_rule(level, parameters)


def load_single_level_rule(path: Path | str) -> ParseRule:
    """Load the sole unparameterized level from a single-level document."""
    source = Path(path)
    document = load_rule_document(source)
    level = _require_single_level(document, source)
    return document.effective_rule(level)


def load_single_parameterized_rule(
    path: Path | str,
    parameters: Parameters,
) -> ParseRule:
    """Load the sole parameterized level from a single-level document."""
    source = Path(path)
    document = load_rule_document(source)
    level = _require_single_level(document, source)
    return document.parameterized_effective_rule(level, parameters)


def load_rules(path: Path | str) -> dict[QuantificationLevel, ParseRule]:
    """Load every unparameterized effective rule from one document."""
    return load_rule_document(path).effective_rules()


def load_parameterized_rules(
    path: Path | str,
    parameters: Parameters,
) -> dict[QuantificationLevel, ParseRule]:
    """Load every parameterized effective rule from one document."""
    return load_rule_document(path).parameterized_effective_rules(parameters)


def _require_document_level(
    document: ParseRuleDocument,
    path: Path,
    level: QuantificationLevel,
) -> None:
    if level not in document.levels:
        raise RuleDocumentError(
            f"{path} has no level {level!r}; available: {list(document.levels)}"
        )


def _require_single_level(
    document: ParseRuleDocument,
    path: Path,
) -> QuantificationLevel:
    if len(document.levels) != 1:
        raise RuleDocumentError(
            f"{path} contains levels {list(document.levels)}; select one explicitly"
        )
    return next(iter(document.levels))


def software_version_matches(pattern: str, version: str) -> bool:
    """Return whether a document's version regex matches a parsed version."""
    try:
        return re.search(pattern, version) is not None
    except re.error as exc:
        raise ValueError(f"invalid software_version regex {pattern!r}") from exc


def _rules_equivalent_without_version(rules: list[ParseRule]) -> bool:
    """Return whether version documents expose the same conversion rule body."""
    normalized = [
        rule.model_dump_json(
            by_alias=True,
            exclude={"software_version", "file_version"},
        )
        for rule in rules
    ]
    return all(item == normalized[0] for item in normalized[1:])


def resolve_rule_locator_for_version(
    software: str,
    level: QuantificationLevel,
    version: str,
) -> RuleLocatorResolution:
    """Resolve one ungated packaged locator for a concrete software version."""
    candidates = _unparameterized_candidates(software, level)
    matching = _matching_version_candidates(candidates, version)
    return _single_locator_resolution(software, level, f"version={version!r}", matching)


def resolve_parameterized_rule_locator_for_version(
    software: str,
    level: QuantificationLevel,
    version: str,
    parameters: Parameters,
) -> RuleLocatorResolution:
    """Resolve one parameter-compatible packaged locator for a concrete version."""
    candidates = _parameterized_candidates(software, level, parameters)
    matching = _matching_version_candidates(candidates, version)
    return _single_locator_resolution(software, level, f"version={version!r}", matching)


def resolve_rule_locator_without_version(
    software: str,
    level: QuantificationLevel,
) -> RuleLocatorResolution:
    """Resolve an unambiguous ungated locator without software-version evidence."""
    candidates = _unparameterized_candidates(software, level)
    return _versionless_unparameterized_resolution(software, level, candidates)


def resolve_parameterized_rule_locator_without_version(
    software: str,
    level: QuantificationLevel,
    parameters: Parameters,
) -> RuleLocatorResolution:
    """Resolve an unambiguous parameter-compatible locator without version evidence."""
    candidates = _parameterized_candidates(software, level, parameters)
    return _versionless_parameterized_resolution(software, level, parameters, candidates)


type _RuleDocumentCandidate = tuple[Path, ParseRuleDocument]


def _unparameterized_candidates(
    software: str,
    level: QuantificationLevel,
) -> list[_RuleDocumentCandidate]:
    candidates: list[_RuleDocumentCandidate] = []
    for path in document_paths_for_software(software):
        document = load_rule_document(path)
        if document.level_is_unconditionally_available(level):
            candidates.append((path, document))
    return candidates


def _parameterized_candidates(
    software: str,
    level: QuantificationLevel,
    parameters: Parameters,
) -> list[_RuleDocumentCandidate]:
    candidates: list[_RuleDocumentCandidate] = []
    for path in document_paths_for_software(software):
        document = load_rule_document(path)
        if document.level_is_available_for(level, parameters):
            candidates.append((path, document))
    return candidates


def _matching_version_candidates(
    candidates: list[_RuleDocumentCandidate],
    version: str,
) -> list[_RuleDocumentCandidate]:
    return [
        candidate
        for candidate in candidates
        if software_version_matches(candidate[1].software_version, version)
    ]


def _single_locator_resolution(
    software: str,
    level: QuantificationLevel,
    evidence: str,
    candidates: list[_RuleDocumentCandidate],
) -> RuleLocatorResolution:
    if len(candidates) == 1:
        return RuleLocator(path=candidates[0][0], level=level)
    if not candidates:
        return RuleLocatorUnavailable(
            f"no packaged rule for software={software!r} level={level!r} {evidence}"
        )
    locators = tuple(RuleLocator(path=path, level=level) for path, _document in candidates)
    return AmbiguousRuleLocators(
        reason=(
            f"{len(candidates)} packaged rules match software={software!r} "
            f"level={level!r} {evidence}"
        ),
        locators=locators,
    )


def _versionless_unparameterized_resolution(
    software: str,
    level: QuantificationLevel,
    candidates: list[_RuleDocumentCandidate],
) -> RuleLocatorResolution:
    exact = _single_locator_resolution(software, level, "without version evidence", candidates)
    if isinstance(exact, RuleLocator | RuleLocatorUnavailable):
        return exact
    rules = [document.effective_rule(level) for _, document in candidates]
    if _rules_equivalent_without_version(rules):
        return RuleLocator(path=candidates[0][0], level=level)
    return exact


def _versionless_parameterized_resolution(
    software: str,
    level: QuantificationLevel,
    parameters: Parameters,
    candidates: list[_RuleDocumentCandidate],
) -> RuleLocatorResolution:
    exact = _single_locator_resolution(software, level, "without version evidence", candidates)
    if isinstance(exact, RuleLocator | RuleLocatorUnavailable):
        return exact
    rules = [document.parameterized_effective_rule(level, parameters) for _, document in candidates]
    if _rules_equivalent_without_version(rules):
        return RuleLocator(path=candidates[0][0], level=level)
    return exact


def load_packaged_rule(
    software: str,
    quantification_level: QuantificationLevel,
) -> ParseRule:
    """Require one unparameterized packaged rule without version evidence."""
    return _require_packaged_rule(resolve_rule_without_version(software, quantification_level))


def load_packaged_rule_for_version(
    software: str,
    quantification_level: QuantificationLevel,
    version: str,
) -> ParseRule:
    """Require one unparameterized packaged rule for a concrete version."""
    return _require_packaged_rule(resolve_rule_for_version(software, quantification_level, version))


def load_parameterized_packaged_rule(
    software: str,
    quantification_level: QuantificationLevel,
    parameters: Parameters,
) -> ParseRule:
    """Require one parameterized packaged rule without version evidence."""
    return _require_packaged_rule(
        resolve_parameterized_rule_without_version(software, quantification_level, parameters)
    )


def load_parameterized_packaged_rule_for_version(
    software: str,
    quantification_level: QuantificationLevel,
    version: str,
    parameters: Parameters,
) -> ParseRule:
    """Require one parameterized packaged rule for a concrete version."""
    return _require_packaged_rule(
        resolve_parameterized_rule_for_version(
            software,
            quantification_level,
            version,
            parameters,
        )
    )


def resolve_rule_for_version(
    software: str,
    quantification_level: QuantificationLevel,
    version: str,
) -> PackagedRuleResolution:
    """Resolve an unparameterized rule for one concrete version."""
    resolution = resolve_rule_locator_for_version(software, quantification_level, version)
    if isinstance(resolution, RuleLocatorUnavailable):
        return PackagedRuleUnavailable(resolution.reason)
    if isinstance(resolution, AmbiguousRuleLocators):
        return AmbiguousPackagedRules(resolution.reason)
    return load_rule(resolution)


def resolve_parameterized_rule_for_version(
    software: str,
    quantification_level: QuantificationLevel,
    version: str,
    parameters: Parameters,
) -> PackagedRuleResolution:
    """Resolve a parameterized rule for one concrete version."""
    resolution = resolve_parameterized_rule_locator_for_version(
        software,
        quantification_level,
        version,
        parameters,
    )
    if isinstance(resolution, RuleLocatorUnavailable):
        return PackagedRuleUnavailable(resolution.reason)
    if isinstance(resolution, AmbiguousRuleLocators):
        return AmbiguousPackagedRules(resolution.reason)
    return load_parameterized_rule(resolution, parameters)


def resolve_rule_without_version(
    software: str,
    quantification_level: QuantificationLevel,
) -> PackagedRuleResolution:
    """Resolve an unparameterized rule without version evidence."""
    resolution = resolve_rule_locator_without_version(software, quantification_level)
    if isinstance(resolution, RuleLocatorUnavailable):
        return PackagedRuleUnavailable(resolution.reason)
    if isinstance(resolution, AmbiguousRuleLocators):
        return AmbiguousPackagedRules(resolution.reason)
    return load_rule(resolution)


def resolve_parameterized_rule_without_version(
    software: str,
    quantification_level: QuantificationLevel,
    parameters: Parameters,
) -> PackagedRuleResolution:
    """Resolve a parameterized rule without version evidence."""
    resolution = resolve_parameterized_rule_locator_without_version(
        software,
        quantification_level,
        parameters,
    )
    if isinstance(resolution, RuleLocatorUnavailable):
        return PackagedRuleUnavailable(resolution.reason)
    if isinstance(resolution, AmbiguousRuleLocators):
        return AmbiguousPackagedRules(resolution.reason)
    return load_parameterized_rule(resolution, parameters)


def _require_packaged_rule(resolution: PackagedRuleResolution) -> ParseRule:
    if isinstance(resolution, PackagedRuleUnavailable):
        raise ValueError(resolution.reason)
    if isinstance(resolution, AmbiguousPackagedRules):
        raise ValueError(resolution.reason)
    return resolution
