"""Backend-independent orchestration for quantification conversion."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass

import pandas as pd

from anndata_proteomics.converters._pieces import ConversionPieces
from anndata_proteomics.converters.assemble import convert_table
from anndata_proteomics.converters.pipeline import (
    NoConvertibleLevelsError,
    ParameterResolution,
    PresentRuleVersion,
    RuleSelection,
    available_parameterized_rules_by_columns,
    available_parameterized_rules_for_version,
    resolve_rule_version,
    select_parameterized_rule_by_columns,
    select_parameterized_rule_for_version,
)
from anndata_proteomics.rules.schema import QuantificationLevel


@dataclass(frozen=True, slots=True)
class LevelConversion:
    """The complete backend-neutral result for one selected level."""

    level: QuantificationLevel
    selection: RuleSelection
    pieces: ConversionPieces


def convert_selected_level(
    data: pd.DataFrame,
    level: QuantificationLevel,
    selection: RuleSelection,
    *,
    strict: bool = False,
) -> LevelConversion:
    """Convert one table using one already-selected effective rule."""
    pieces = convert_table(data, selection.rule, strict=strict)
    return LevelConversion(level=level, selection=selection, pieces=pieces)


def convert_selected_levels(
    data: pd.DataFrame,
    selections: Mapping[QuantificationLevel, RuleSelection],
    *,
    strict: bool = False,
) -> dict[QuantificationLevel, LevelConversion]:
    """Convert every explicitly selected level from one vendor table."""
    if not selections:
        raise NoConvertibleLevelsError("no levels supplied")
    return {
        level: convert_selected_level(data, level, selection, strict=strict)
        for level, selection in selections.items()
    }


def select_rule_from_parameters(
    headers: Iterable[str],
    slug: str,
    level: QuantificationLevel,
    resolution: ParameterResolution,
) -> RuleSelection:
    """Select one rule using one completely parsed parameter result."""
    version = resolve_rule_version(resolution, slug)
    if isinstance(version, PresentRuleVersion):
        return select_parameterized_rule_for_version(
            slug,
            level,
            version.value,
            headers,
            resolution.parameters,
        )
    return select_parameterized_rule_by_columns(
        slug,
        level,
        headers,
        resolution.parameters,
    )


def select_rules_from_parameters(
    headers: Iterable[str],
    slug: str,
    resolution: ParameterResolution,
) -> dict[QuantificationLevel, RuleSelection]:
    """Select every available level using one parsed parameter result."""
    version = resolve_rule_version(resolution, slug)
    if isinstance(version, PresentRuleVersion):
        return available_parameterized_rules_for_version(
            slug,
            version.value,
            headers,
            resolution.parameters,
        )
    return available_parameterized_rules_by_columns(
        slug,
        headers,
        resolution.parameters,
    )


def convert_level_from_parameters(
    data: pd.DataFrame,
    slug: str,
    level: QuantificationLevel,
    resolution: ParameterResolution,
    *,
    strict: bool = False,
) -> LevelConversion:
    """Select and convert one level from a parsed parameter result."""
    selection = select_rule_from_parameters(data.columns, slug, level, resolution)
    return convert_selected_level(data, level, selection, strict=strict)
