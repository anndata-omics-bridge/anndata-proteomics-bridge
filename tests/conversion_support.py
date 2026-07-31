"""Test composition for the pure conversion and AnnData adapter boundaries."""

from __future__ import annotations

import pandas as pd
from anndata import AnnData

from anndata_proteomics.adapters.anndata.conversion import to_anndata, write_parameter_resolution
from anndata_proteomics.converters.pipeline import (
    ParameterResolution,
    RuleSelection,
    select_rule_for_version,
)
from anndata_proteomics.rules.schema import ParseRule, QuantificationLevel
from anndata_proteomics.workflows.conversion import (
    convert_level_from_parameters,
    convert_selected_level,
)


def convert_to_anndata(
    data: pd.DataFrame,
    rule: ParseRule,
    *,
    strict: bool = False,
) -> AnnData:
    """Compose the conversion computation and AnnData persistence in tests."""
    conversion = convert_selected_level(
        data,
        rule.quantification_level,
        RuleSelection(rule, "rule_config"),
        strict=strict,
    )
    return to_anndata(conversion)


def convert_versioned_level_to_anndata(
    data: pd.DataFrame,
    slug: str,
    level: QuantificationLevel,
    version: str,
    *,
    strict: bool = False,
) -> AnnData:
    """Compose version-based selection, conversion, and AnnData persistence."""
    conversion = convert_selected_level(
        data,
        level,
        select_rule_for_version(slug, level, version, data.columns),
        strict=strict,
    )
    return to_anndata(conversion)


def convert_parameterized_level_to_anndata(
    data: pd.DataFrame,
    slug: str,
    level: QuantificationLevel,
    resolution: ParameterResolution,
    *,
    strict: bool = False,
) -> AnnData:
    """Compose parameter-based selection, conversion, and AnnData persistence."""
    conversion = convert_level_from_parameters(
        data,
        slug,
        level,
        resolution,
        strict=strict,
    )
    result = to_anndata(conversion)
    write_parameter_resolution(result, resolution)
    return result
