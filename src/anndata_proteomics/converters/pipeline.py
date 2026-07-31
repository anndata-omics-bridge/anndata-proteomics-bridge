"""Typed parsing-rule selection for vendor quantification tables.

This module knows about vendor headers, rule documents, and parsed search
parameters.  It does not know about AnnData, MuData, or persistence.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from anndata_proteomics.converters.recognize import matches
from anndata_proteomics.params.model import Parameters
from anndata_proteomics.params.registry import parse_params
from anndata_proteomics.rules.loader import (
    AmbiguousPackagedRules,
    PackagedRuleUnavailable,
    load_parameterized_rule,
    load_rule,
    load_rule_document,
    resolve_parameterized_rule_for_version,
    resolve_rule_for_version,
)
from anndata_proteomics.rules.registry import iter_packaged_rules
from anndata_proteomics.rules.schema import ParseRule, QuantificationLevel

LEVELS: tuple[QuantificationLevel, ...] = (
    "ion",
    "peptidoform",
    "peptide",
    "protein",
    "fragment",
)
MUDATA = "mudata"

ParameterVersionStatus = Literal["present", "missing"]
RuleSelectionMethod = Literal[
    "software_version",
    "columns",
    "rule_config",
]


class RuleUnavailableError(ValueError):
    """A requested quantification level is not available for this input."""


class AmbiguousRuleError(ValueError):
    """More than one parsing rule is valid for the supplied evidence."""


class NoConvertibleLevelsError(ValueError):
    """No quantification level is available for a compound conversion."""


@dataclass(frozen=True, slots=True)
class PresentRuleVersion:
    """A software version was present in parsed parameter metadata."""

    value: str


@dataclass(frozen=True, slots=True)
class MissingRuleVersion:
    """Parsed parameter metadata genuinely contained no applicable version."""


type RuleVersion = PresentRuleVersion | MissingRuleVersion


@dataclass(frozen=True, slots=True)
class ParameterResolution:
    """One parsed parameter file and its typed primary-version result."""

    source_path: Path
    parameters: Parameters
    version: RuleVersion


@dataclass(frozen=True, slots=True)
class RuleSelection:
    """One effective parsing rule and the evidence that selected it."""

    rule: ParseRule
    method: RuleSelectionMethod


@dataclass(frozen=True, slots=True)
class RuleUnavailable:
    """No parsing rule fits one requested quantification level."""

    reason: str


@dataclass(frozen=True, slots=True)
class AmbiguousRules:
    """Several parsing rules fit evidence that should select exactly one."""

    reason: str


type RuleLookup = RuleSelection | RuleUnavailable | AmbiguousRules


@dataclass(frozen=True, slots=True)
class RecognizedSoftware:
    """All matching packaged rules identify one vendor slug."""

    slug: str


@dataclass(frozen=True, slots=True)
class UnrecognizedSoftware:
    """No unique vendor slug can be inferred from the headers."""


type SoftwareRecognition = RecognizedSoftware | UnrecognizedSoftware


def software_slug(software_name: str) -> str:
    """Map a catalog software name such as ``DIA-NN`` to its rule slug."""
    return re.sub(r"[^a-z0-9]", "", software_name.lower())


def version_status(version: RuleVersion) -> ParameterVersionStatus:
    """Return the serialized status for a typed rule-version result."""
    if isinstance(version, PresentRuleVersion):
        return "present"
    return "missing"


def resolve_parameters(param_path: Path | str, slug: str) -> ParameterResolution:
    """Parse one parameter file once and retain its typed version result."""
    source_path = Path(param_path)
    parameters = parse_params(source_path, software=slug)
    parsed_version = parameters.software_version
    version: RuleVersion = (
        MissingRuleVersion() if parsed_version is None else PresentRuleVersion(parsed_version)
    )
    return ParameterResolution(
        source_path=source_path,
        parameters=parameters,
        version=version,
    )


def resolve_rule_version(
    parameter_resolution: ParameterResolution,
    rule_slug: str,
) -> RuleVersion:
    """Resolve the version belonging to a primary or quantification software."""
    parameters = parameter_resolution.parameters
    candidates = (
        (parameters.software_name, parameters.software_version),
        (
            parameters.quantification_software,
            parameters.quantification_software_version,
        ),
    )
    for software_name, version in candidates:
        if software_name is not None and software_slug(software_name) == rule_slug:
            if version is None:
                return MissingRuleVersion()
            return PresentRuleVersion(version)
    if not any(software_name for software_name, _version in candidates):
        return parameter_resolution.version
    return MissingRuleVersion()


def recognize_software(headers: Iterable[str]) -> SoftwareRecognition:
    """Infer a vendor slug when every matching packaged rule agrees on it."""
    header_set = set(headers)
    slugs = {
        software_slug(rule.software_name)
        for path in iter_packaged_rules()
        if matches(header_set, rule := load_rule(path))
    }
    if len(slugs) == 1:
        return RecognizedSoftware(next(iter(slugs)))
    return UnrecognizedSoftware()


def _rule_matches_headers(headers: Iterable[str], rule: ParseRule) -> bool:
    """Return whether all input columns required by ``rule`` are available."""
    header_set = set(headers)
    label_column = (
        rule.fragments.label_column
        if rule.fragments is not None and rule.fragments.label_strategy == "column"
        else None
    )
    return matches(header_set, rule) and (label_column is None or label_column in header_set)


def _column_matching_rules(
    slug: str,
    level: QuantificationLevel,
    headers: Iterable[str],
) -> list[ParseRule]:
    """Find unparameterized packaged variants compatible with the headers."""
    header_set = set(headers)
    candidates: list[ParseRule] = []
    for locator in iter_packaged_rules():
        if locator.level != level:
            continue
        document = load_rule_document(locator.path)
        if not document.level_is_unconditionally_available(level):
            continue
        rule = load_rule(locator)
        if software_slug(rule.software_name) == slug and _rule_matches_headers(header_set, rule):
            candidates.append(rule)
    return candidates


def _parameterized_column_matching_rules(
    slug: str,
    level: QuantificationLevel,
    headers: Iterable[str],
    parameters: Parameters,
) -> list[ParseRule]:
    """Find parameter-gated packaged variants compatible with the headers."""
    header_set = set(headers)
    candidates: list[ParseRule] = []
    for locator in iter_packaged_rules():
        if locator.level != level:
            continue
        document = load_rule_document(locator.path)
        if not document.level_is_available_for(level, parameters):
            continue
        rule = load_parameterized_rule(locator, parameters)
        if software_slug(rule.software_name) == slug and _rule_matches_headers(header_set, rule):
            candidates.append(rule)
    return candidates


def _lookup_by_columns(
    slug: str,
    level: QuantificationLevel,
    candidates: list[ParseRule],
) -> RuleLookup:
    """Classify a column-driven rule lookup without exception-based branching."""
    if len(candidates) == 1:
        return RuleSelection(candidates[0], "columns")
    if not candidates:
        return RuleUnavailable(
            f"{slug} {level}: no rule matches the file columns while software version is missing"
        )
    return AmbiguousRules(
        f"{slug} {level}: {len(candidates)} rules match the file columns while software "
        "version is missing; provide a version or explicit rule config"
    )


def find_rule_by_columns(
    slug: str,
    level: QuantificationLevel,
    headers: Iterable[str],
) -> RuleLookup:
    """Look up one rule from headers without parsed search parameters."""
    return _lookup_by_columns(slug, level, _column_matching_rules(slug, level, headers))


def find_parameterized_rule_by_columns(
    slug: str,
    level: QuantificationLevel,
    headers: Iterable[str],
    parameters: Parameters,
) -> RuleLookup:
    """Look up one parameter-gated rule from headers."""
    candidates = _parameterized_column_matching_rules(slug, level, headers, parameters)
    return _lookup_by_columns(slug, level, candidates)


def _lookup_resolved_rule(
    slug: str,
    level: QuantificationLevel,
    version: str,
    headers: Iterable[str],
    rule: ParseRule,
) -> RuleLookup:
    """Classify one version-resolved rule against vendor headers."""
    if not _rule_matches_headers(headers, rule):
        return RuleUnavailable(
            f"{slug} {level}: file columns don't match the rule for software version "
            f"{version!r} — verify the version / provide the right param file"
        )
    return RuleSelection(rule, "software_version")


def find_rule_for_version(
    slug: str,
    level: QuantificationLevel,
    version: str,
    headers: Iterable[str],
) -> RuleLookup:
    """Look up an unparameterized rule for one concrete software version."""
    rule = resolve_rule_for_version(slug, level, version)
    if isinstance(rule, PackagedRuleUnavailable):
        return RuleUnavailable(f"{slug} {level}: no rule covers software version {version!r}")
    if isinstance(rule, AmbiguousPackagedRules):
        return AmbiguousRules(rule.reason)
    return _lookup_resolved_rule(slug, level, version, headers, rule)


def find_parameterized_rule_for_version(
    slug: str,
    level: QuantificationLevel,
    version: str,
    headers: Iterable[str],
    parameters: Parameters,
) -> RuleLookup:
    """Look up a parameter-gated rule for one concrete software version."""
    rule = resolve_parameterized_rule_for_version(
        slug,
        level,
        version,
        parameters,
    )
    if isinstance(rule, PackagedRuleUnavailable):
        return RuleUnavailable(f"{slug} {level}: no rule covers software version {version!r}")
    if isinstance(rule, AmbiguousPackagedRules):
        return AmbiguousRules(rule.reason)
    return _lookup_resolved_rule(slug, level, version, headers, rule)


def require_rule_selection(lookup: RuleLookup) -> RuleSelection:
    """Return a selected rule or raise the lookup's precise hard failure."""
    if isinstance(lookup, RuleSelection):
        return lookup
    if isinstance(lookup, AmbiguousRules):
        raise AmbiguousRuleError(lookup.reason)
    raise RuleUnavailableError(lookup.reason)


def select_rule_by_columns(
    slug: str,
    level: QuantificationLevel,
    headers: Iterable[str],
) -> RuleSelection:
    """Require the unique unparameterized rule selected by headers."""
    return require_rule_selection(find_rule_by_columns(slug, level, headers))


def select_parameterized_rule_by_columns(
    slug: str,
    level: QuantificationLevel,
    headers: Iterable[str],
    parameters: Parameters,
) -> RuleSelection:
    """Require the unique parameter-gated rule selected by headers."""
    lookup = find_parameterized_rule_by_columns(slug, level, headers, parameters)
    return require_rule_selection(lookup)


def select_rule_for_version(
    slug: str,
    level: QuantificationLevel,
    version: str,
    headers: Iterable[str],
) -> RuleSelection:
    """Require an unparameterized rule selected by a concrete version."""
    return require_rule_selection(find_rule_for_version(slug, level, version, headers))


def select_parameterized_rule_for_version(
    slug: str,
    level: QuantificationLevel,
    version: str,
    headers: Iterable[str],
    parameters: Parameters,
) -> RuleSelection:
    """Require a parameter-gated rule selected by a concrete version."""
    lookup = find_parameterized_rule_for_version(slug, level, version, headers, parameters)
    return require_rule_selection(lookup)


def _available_rule_selections(
    lookup: Callable[[QuantificationLevel], RuleLookup],
) -> dict[QuantificationLevel, RuleSelection]:
    """Collect available selections while preserving ambiguous hard failures."""
    selections: dict[QuantificationLevel, RuleSelection] = {}
    for level in LEVELS:
        result = lookup(level)
        if isinstance(result, RuleSelection):
            selections[level] = result
        elif isinstance(result, AmbiguousRules):
            raise AmbiguousRuleError(result.reason)
    return selections


def available_rules_by_columns(
    slug: str,
    headers: Iterable[str],
) -> dict[QuantificationLevel, RuleSelection]:
    """Return every unparameterized level uniquely selected by headers."""
    header_set = set(headers)
    return _available_rule_selections(lambda level: find_rule_by_columns(slug, level, header_set))


def available_parameterized_rules_by_columns(
    slug: str,
    headers: Iterable[str],
    parameters: Parameters,
) -> dict[QuantificationLevel, RuleSelection]:
    """Return every parameter-gated level uniquely selected by headers."""
    header_set = set(headers)
    return _available_rule_selections(
        lambda level: find_parameterized_rule_by_columns(
            slug,
            level,
            header_set,
            parameters,
        )
    )


def available_rules_for_version(
    slug: str,
    version: str,
    headers: Iterable[str],
) -> dict[QuantificationLevel, RuleSelection]:
    """Return every unparameterized level available for one version."""
    header_set = set(headers)
    return _available_rule_selections(
        lambda level: find_rule_for_version(slug, level, version, header_set)
    )


def available_parameterized_rules_for_version(
    slug: str,
    version: str,
    headers: Iterable[str],
    parameters: Parameters,
) -> dict[QuantificationLevel, RuleSelection]:
    """Return every parameter-gated level available for one version."""
    header_set = set(headers)
    return _available_rule_selections(
        lambda level: find_parameterized_rule_for_version(
            slug,
            level,
            version,
            header_set,
            parameters,
        )
    )


def matching_rules(
    rules: Mapping[QuantificationLevel, ParseRule],
    headers: Iterable[str],
) -> dict[QuantificationLevel, ParseRule]:
    """Return explicit-document levels whose required columns match a table."""
    header_set = set(headers)
    return {level: rule for level, rule in rules.items() if _rule_matches_headers(header_set, rule)}
