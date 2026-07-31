"""Param-driven conversion orchestration: vendor file + params → AnnData / MuData.

The non-UI core that backs the ``apb convert`` CLI. Given a vendor table plus the software
version (parsed from the parameter file), it selects the matching parsing-rule variant, converts
each quantification level, and (for the default target) wraps them into a MuData on a shared run
axis. No GUI / marimo / test-data-cache dependency — this is plain library code.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import mudata
import pandas as pd
from anndata import AnnData
from loguru import logger
from mudata import MuData

from anndata_proteomics.converters import assemble
from anndata_proteomics.converters.recognize import matches
from anndata_proteomics.params.anndata_io import write_search_parameters
from anndata_proteomics.params.model import Parameters
from anndata_proteomics.params.registry import parse_params
from anndata_proteomics.rules.loader import (
    load_rule,
    load_rule_document,
    resolve_rule_for_version,
)
from anndata_proteomics.rules.registry import iter_packaged_rules
from anndata_proteomics.rules.schema import ParseRule, QuantificationLevel

# Quantification levels, coarse to fine. Not every vendor exposes every level.
LEVELS: tuple[QuantificationLevel, ...] = (
    "ion",
    "peptidoform",
    "peptide",
    "protein",
    "fragment",
)
MUDATA = "mudata"
# Per-level var_names prefix so modalities don't collide on the global axis.
_PREFIX = {
    "fragment": "frg:",
    "ion": "ion:",
    "peptidoform": "pfm:",
    "peptide": "pep:",
    "protein": "prt:",
}

ParameterVersionStatus = Literal["present", "missing"]
RuleSelectionMethod = Literal[
    "software_version",
    "columns",
    "rule_config",
]


class RuleUnavailableError(ValueError):
    """A requested quantification level is not available for this input."""


class NoConvertibleLevelsError(ValueError):
    """No quantification level is available for a compound conversion."""


@dataclass(frozen=True)
class ParameterResolution:
    """One successful parameter parse with explicit version availability."""

    source_path: Path
    parameters: Parameters
    version: str | None
    version_status: ParameterVersionStatus


def software_slug(software_name: str) -> str:
    """Map a catalog ``software_name`` (e.g. "DIA-NN") to a parsing-rule vendor slug ("diann")."""
    return re.sub(r"[^a-z0-9]", "", software_name.lower())


def resolve_rule_version(
    parameter_resolution: ParameterResolution,
    rule_slug: str,
) -> tuple[str | None, ParameterVersionStatus]:
    """Resolve a parsing-rule version from primary or quantification software metadata."""
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
            return version, "present" if version is not None else "missing"
    if not any(software_name for software_name, _version in candidates):
        return parameter_resolution.version, parameter_resolution.version_status
    return None, "missing"


def _effective_rule_version(
    rule_slug: str,
    version: str | None,
    parameter_resolution: ParameterResolution | None,
) -> tuple[str | None, ParameterVersionStatus | None]:
    if parameter_resolution is None:
        return version, None
    resolved_version, version_status = resolve_rule_version(
        parameter_resolution,
        rule_slug,
    )
    if resolved_version is not None:
        return resolved_version, version_status
    if version is not None and version_status == "missing":
        return version, "present"
    return version, version_status


def recognize_software(headers: Iterable[str]) -> str | None:
    """Infer the vendor slug from a file's column headers.

    Unlike ``recognize`` (which needs a *unique* rule match), a single vendor file can match
    several of its own level rules (e.g. a DIA-NN report matches ion/fragment/protein). This
    returns the vendor slug when every matching packaged rule belongs to the same vendor, else
    ``None`` (zero matches, or an ambiguous match spanning multiple vendors).
    """
    header_set = set(headers)
    slugs = {
        software_slug(rule.software_name)
        for path in iter_packaged_rules()
        if matches(header_set, rule := load_rule(path))
    }
    return next(iter(slugs)) if len(slugs) == 1 else None


def param_version(param_path: Path | None, slug: str) -> str | None:
    """Software version parsed from the parameter file, or ``None`` when absent."""
    if param_path is None:
        return None
    return resolve_parameters(param_path, slug).version


def resolve_parameters(param_path: Path | str, slug: str) -> ParameterResolution:
    """Parse parameters once and report whether the parsed version is present."""
    source_path = Path(param_path)
    parameters = parse_params(source_path, software=slug)
    version = parameters.software_version
    return ParameterResolution(
        source_path=source_path,
        parameters=parameters,
        version=version,
        version_status="present" if version is not None else "missing",
    )


def _rule_matches_headers(headers: Iterable[str], rule: ParseRule) -> bool:
    """Return whether all input columns required by ``rule`` are available."""
    header_set = set(headers)
    label_column = (
        rule.fragments.label_column
        if rule.fragments is not None and rule.fragments.label_strategy == "column"
        else None
    )
    return matches(header_set, rule) and (label_column is None or label_column in header_set)


def _column_matching_rule_variants(
    slug: str,
    level: QuantificationLevel,
    headers: Iterable[str],
    search_parameters: Parameters | None = None,
) -> list[ParseRule]:
    """Return packaged variants for ``slug``/``level`` that fit the input headers."""
    header_set = set(headers)
    candidates = []
    for locator in iter_packaged_rules():
        if locator.level != level:
            continue
        # Same availability gate the version resolver applies: two levels of one document can
        # share both the version regex and the header schema, so headers cannot separate them.
        if not load_rule_document(locator.path).level_is_available(level, search_parameters):
            continue
        rule = load_rule(locator, search_parameters=search_parameters)
        if software_slug(rule.software_name) != slug:
            continue
        if _rule_matches_headers(header_set, rule):
            candidates.append(rule)
    return candidates


def _select_rule(
    slug: str,
    level: QuantificationLevel,
    version: str | None,
    headers: Iterable[str],
    *,
    version_status: ParameterVersionStatus | None = None,
    search_parameters: Parameters | None = None,
) -> tuple[ParseRule, RuleSelectionMethod]:
    """Resolve one rule together with the method used to select it."""
    status = version_status or ("present" if version is not None else "missing")
    header_set = set(headers)

    if status == "missing":
        candidates = _column_matching_rule_variants(
            slug,
            level,
            header_set,
            search_parameters,
        )
        if len(candidates) == 1:
            return candidates[0], "columns"
        if not candidates:
            raise RuleUnavailableError(
                f"{slug} {level}: no rule matches the file columns while software "
                "version is missing"
            )
        raise ValueError(
            f"{slug} {level}: {len(candidates)} rules match the file columns while "
            "software version is missing; provide a version or explicit rule config"
        )

    rule = resolve_rule_for_version(
        slug,
        level,
        version,
        search_parameters=search_parameters,
    )
    if rule is None:
        raise RuleUnavailableError(f"{slug} {level}: no rule covers software version {version!r}")
    if not _rule_matches_headers(header_set, rule):
        raise RuleUnavailableError(
            f"{slug} {level}: file columns don't match the rule for software version "
            f"{version!r} — verify the version / provide the right param file"
        )
    return rule, "software_version"


def select_rule(
    slug: str,
    level: QuantificationLevel,
    version: str | None,
    headers: Iterable[str],
    *,
    version_status: ParameterVersionStatus | None = None,
    search_parameters: Parameters | None = None,
) -> ParseRule:
    """Resolve the rule for (slug, level) at ``version`` and validate it against ``headers``.

    A genuinely missing version selects the unique column-compatible variant. A present,
    uncovered version remains a hard failure.
    """
    return _select_rule(
        slug,
        level,
        version,
        headers,
        version_status=version_status,
        search_parameters=search_parameters,
    )[0]


def _available_rule_selections(
    slug: str,
    version: str | None,
    headers: Iterable[str],
    *,
    version_status: ParameterVersionStatus | None = None,
    search_parameters: Parameters | None = None,
) -> dict[QuantificationLevel, tuple[ParseRule, RuleSelectionMethod]]:
    """Resolve all available levels while propagating unexpected rule failures."""
    header_set = set(headers)
    selections: dict[
        QuantificationLevel,
        tuple[ParseRule, RuleSelectionMethod],
    ] = {}
    for level in LEVELS:
        try:
            selections[level] = _select_rule(
                slug,
                level,
                version,
                header_set,
                version_status=version_status,
                search_parameters=search_parameters,
            )
        except RuleUnavailableError:
            continue
    return selections


def convertible_levels(
    slug: str,
    version: str | None,
    headers: Iterable[str],
    *,
    version_status: ParameterVersionStatus | None = None,
    search_parameters: Parameters | None = None,
) -> list[QuantificationLevel]:
    """Levels whose version-selected rule both exists and matches this file's columns."""
    return list(
        _available_rule_selections(
            slug,
            version,
            headers,
            version_status=version_status,
            search_parameters=search_parameters,
        )
    )


def available_targets(
    slug: str,
    version: str | None,
    headers: Iterable[str],
    *,
    version_status: ParameterVersionStatus | None = None,
    search_parameters: Parameters | None = None,
) -> list[str]:
    """Convertible levels plus the all-level MuData target when any level resolves."""
    levels = convertible_levels(
        slug,
        version,
        headers,
        version_status=version_status,
        search_parameters=search_parameters,
    )
    targets = [str(level) for level in levels]
    if levels:
        targets.append(MUDATA)
    return targets


def matching_rules(
    rules: Mapping[QuantificationLevel, ParseRule],
    headers: Iterable[str],
) -> dict[QuantificationLevel, ParseRule]:
    """Return document levels whose required source columns match a vendor table."""
    header_set = set(headers)
    matched = {}
    for level, rule in rules.items():
        if _rule_matches_headers(header_set, rule):
            matched[level] = rule
    return matched


def convert_level(
    df: pd.DataFrame,
    slug: str,
    level: QuantificationLevel,
    version: str | None,
    *,
    params_path: Path | str | None = None,
    parameter_resolution: ParameterResolution | None = None,
    strict: bool = False,
) -> AnnData:
    if parameter_resolution is None and params_path is not None:
        parameter_resolution = resolve_parameters(params_path, slug)
    version, version_status = _effective_rule_version(
        slug,
        version,
        parameter_resolution,
    )
    search_parameters = (
        parameter_resolution.parameters if parameter_resolution is not None else None
    )
    rule, selection_method = _select_rule(
        slug,
        level,
        version,
        df.columns,
        version_status=version_status,
        search_parameters=search_parameters,
    )
    adata = assemble.convert(
        df,
        rule,
        params_path=None if parameter_resolution is not None else params_path,
        strict=strict,
    )
    if parameter_resolution is not None:
        attach_parameter_resolution(
            adata,
            parameter_resolution,
            selection_method=selection_method,
            warn_missing=True,
        )
    else:
        set_rule_selection_method(adata, selection_method)
    logger.info(f"  {level}: {adata.shape[0]} obs × {adata.shape[1]} var")
    return adata


def build_mudata(
    df: pd.DataFrame,
    slug: str,
    version: str | None,
    *,
    params_path: Path | str | None = None,
    parameter_resolution: ParameterResolution | None = None,
    strict: bool = False,
) -> MuData:
    """Build a MuData over the levels whose version-selected rule fits this file (shared run axis).

    Levels the version doesn't provide (e.g. fragment on DIA-NN 2.x) are skipped, not failed.
    """
    if parameter_resolution is None and params_path is not None:
        parameter_resolution = resolve_parameters(params_path, slug)
    version, version_status = _effective_rule_version(
        slug,
        version,
        parameter_resolution,
    )
    search_parameters = (
        parameter_resolution.parameters if parameter_resolution is not None else None
    )
    selections = _available_rule_selections(
        slug,
        version,
        df.columns,
        version_status=version_status,
        search_parameters=search_parameters,
    )
    resolvable = set(selections)
    if not resolvable:
        raise NoConvertibleLevelsError(
            f"{slug}: no level resolves for software version {version!r}"
        )
    skipped = [level for level in LEVELS if level not in resolvable]
    if skipped:
        logger.info(f"skipping levels not provided by software version {version!r}: {skipped}")
    rules: dict[QuantificationLevel, ParseRule] = {
        level: selection[0] for level, selection in selections.items()
    }
    selection_method = next(iter(selections.values()))[1]
    return build_mudata_from_rules(
        df,
        rules,
        params_path=params_path,
        parameter_resolution=parameter_resolution,
        rule_selection_method=selection_method,
        software=slug,
        strict=strict,
    )


def build_mudata_from_rules(
    df: pd.DataFrame,
    rules: Mapping[QuantificationLevel, ParseRule],
    *,
    params_path: Path | str | None = None,
    parameter_resolution: ParameterResolution | None = None,
    rule_selection_method: RuleSelectionMethod | None = None,
    software: str | None = None,
    strict: bool = False,
) -> MuData:
    """Build MuData from already selected effective rules."""
    if not rules:
        raise ValueError("no levels supplied")
    mods: dict[str, AnnData] = {}
    for level in LEVELS:
        if level not in rules:
            continue
        logger.info(f"converting level: {level}")
        adata = assemble.convert(
            df.copy(),
            rules[level],
            params_path=None if parameter_resolution is not None else params_path,
            strict=strict,
        )
        if parameter_resolution is not None:
            attach_parameter_resolution(
                adata,
                parameter_resolution,
                selection_method=rule_selection_method,
                warn_missing=False,
            )
        elif rule_selection_method is not None:
            set_rule_selection_method(adata, rule_selection_method)
        logger.info(f"  {level}: {adata.shape[0]} obs × {adata.shape[1]} var")
        adata.var_names = [_PREFIX[level] + str(v) for v in adata.var_names]
        mods[level] = adata
    # Adopt the mudata 0.4 default now (no auto-pull of per-modality obs/var into the global
    # frames); each modality keeps its own obs/var. Silences the 0.3 FutureWarning.
    with mudata.set_options(pull_on_update=False):
        md = MuData(mods, axis=0)
    if parameter_resolution is not None:
        attach_parameter_resolution(
            md,
            parameter_resolution,
            selection_method=rule_selection_method,
            warn_missing=True,
        )
    elif params_path is not None and software is not None:
        params = parse_params(params_path, software=software)
        write_search_parameters(md, params, source_path=str(params_path))
    if rule_selection_method is not None:
        set_rule_selection_method(md, rule_selection_method)
    logger.info(
        f"  MuData: {md.n_obs} obs × {sum(a.n_vars for a in mods.values())} var, {len(mods)} mods"
    )
    return md


def set_rule_selection_method(
    target: AnnData | MuData,
    selection_method: RuleSelectionMethod,
) -> None:
    """Store rule-selection provenance on one AnnData or MuData container."""
    target.uns.setdefault("anndata_proteomics", {})
    target.uns["anndata_proteomics"]["rule_selection_method"] = selection_method


def attach_parameter_resolution(
    target: AnnData | MuData,
    resolution: ParameterResolution,
    *,
    selection_method: RuleSelectionMethod | None,
    warn_missing: bool,
) -> None:
    """Attach one already-parsed parameter result and its selection provenance."""
    target.uns.setdefault("anndata_proteomics", {})
    metadata = target.uns["anndata_proteomics"]
    metadata["search_parameters_version_status"] = resolution.version_status
    metadata["search_parameters_path"] = str(resolution.source_path)
    if selection_method is not None:
        metadata["rule_selection_method"] = selection_method
    write_search_parameters(
        target,
        resolution.parameters,
        source_path=str(resolution.source_path),
    )
    if resolution.version_status == "missing" and warn_missing:
        logger.warning(
            "no software version in search parameters {}; selected rule by columns",
            resolution.source_path,
        )
