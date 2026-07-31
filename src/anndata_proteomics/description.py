"""Typed, backend-neutral construction of APB container descriptions."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from anndata_proteomics.params.model import Parameters
from anndata_proteomics.rules.schema import ColumnGroup, ParseRule, QuantificationLevel
from anndata_proteomics.serialization import JsonObject, JsonValue, to_json_compatible

_VIEW_SCHEMA_VERSION = "1"


@dataclass(frozen=True, slots=True)
class MissingSearchParameters:
    """No search parameters were stored in the APB namespace."""


@dataclass(frozen=True, slots=True)
class MissingRuleMetadata:
    """No conversion rule was stored in the APB namespace."""


@dataclass(frozen=True, slots=True)
class MissingQcMetadata:
    """No QC result was stored in the APB namespace."""


@dataclass(frozen=True, slots=True)
class MissingProteoBenchMetadata:
    """No ProteoBench score result was stored in the APB namespace."""


@dataclass(frozen=True, slots=True)
class MissingQuantificationLevel:
    """No quantification level was stored at this description scope."""


@dataclass(frozen=True, slots=True)
class MissingSoftwareName:
    """No quantification software name was stored at this description scope."""


type DescriptionQuantificationLevel = QuantificationLevel | MissingQuantificationLevel
type DescriptionSoftwareName = str | MissingSoftwareName


class DescriptionConversionMetadata(BaseModel):
    """Validated conversion fields extracted from the APB namespace."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    schema_version: str | None = Field(default=None, min_length=1)
    quantification_level: QuantificationLevel | None = None
    software_name: str | None = Field(default=None, min_length=1)
    rule_selection_method: Literal["software_version", "columns", "rule_config"] | None = None
    search_parameters_path: str | None = Field(default=None, min_length=1)
    search_parameters_version_status: Literal["present", "missing"] | None = None


@dataclass(frozen=True, slots=True)
class DescriptionMetadata:
    """Validated APB metadata required by description calculations."""

    quantification_level: DescriptionQuantificationLevel
    software_name: DescriptionSoftwareName
    conversion: DescriptionConversionMetadata
    search_parameters: Parameters | MissingSearchParameters
    rule: ParseRule | MissingRuleMetadata
    annotations: JsonObject
    qc: JsonValue | MissingQcMetadata
    proteobench: JsonObject | MissingProteoBenchMetadata


@dataclass(frozen=True, slots=True)
class AnnDataDescriptionSource:
    """The exact non-quantitative fields needed to describe one data matrix."""

    n_runs: int
    n_features: int
    layers: tuple[str, ...]
    metadata: DescriptionMetadata


@dataclass(frozen=True, slots=True)
class MuDataDescriptionSource:
    """The exact non-quantitative fields needed to describe a modality collection."""

    n_runs: int
    n_features: int
    modalities: Mapping[str, AnnDataDescriptionSource]
    metadata: DescriptionMetadata


def calculate_anndata_description(source: AnnDataDescriptionSource) -> JsonObject:
    """Calculate one AnnData description from extracted metadata."""
    metadata = source.metadata
    parameters = metadata.search_parameters
    software_version: JsonValue = None
    if isinstance(parameters, Parameters):
        software_version = to_json_compatible(parameters.software_version)

    layers: list[JsonValue] = []
    layers.extend(source.layers)
    quantification: JsonObject = {
        "n_runs": source.n_runs,
        "n_features": source.n_features,
        "level": _quantification_level_value(metadata.quantification_level),
        "software_name": _software_name_value(metadata.software_name),
        "software_version": software_version,
        "layers": layers,
    }
    result: JsonObject = {
        "schema_version": _VIEW_SCHEMA_VERSION,
        "container_type": "anndata",
        "quantification": quantification,
    }

    column_mapping = _column_mapping(source)
    if column_mapping:
        result["column_mapping"] = column_mapping
    _add_metadata_views(result, metadata)
    return result


def calculate_mudata_description(source: MuDataDescriptionSource) -> JsonObject:
    """Calculate one MuData description from extracted root and modality metadata."""
    modalities: list[JsonValue] = []
    modalities.extend(source.modalities)
    descriptions: JsonObject = {
        name: calculate_anndata_description(modality)
        for name, modality in source.modalities.items()
    }
    quantification: JsonObject = {
        "n_runs": source.n_runs,
        "n_features": source.n_features,
        "modalities": modalities,
    }
    result: JsonObject = {
        "schema_version": _VIEW_SCHEMA_VERSION,
        "container_type": "mudata",
        "quantification": quantification,
        "modalities": descriptions,
    }
    _add_metadata_views(result, source.metadata)
    return result


def _add_metadata_views(
    result: JsonObject,
    metadata: DescriptionMetadata,
) -> None:
    conversion = metadata.conversion.model_dump(mode="json", exclude_none=True)
    if conversion:
        result["conversion"] = to_json_compatible(conversion)

    parameters = metadata.search_parameters
    if isinstance(parameters, Parameters):
        result["search_parameters"] = to_json_compatible(parameters.model_dump(mode="json"))

    if metadata.annotations:
        result["annotations"] = metadata.annotations

    if not isinstance(metadata.qc, MissingQcMetadata):
        result["qc"] = metadata.qc

    if not isinstance(metadata.proteobench, MissingProteoBenchMetadata):
        result["proteobench"] = metadata.proteobench


def _column_mapping(source: AnnDataDescriptionSource) -> JsonObject:
    """Describe where the effective rule placed vendor and computed columns."""
    rule = source.metadata.rule
    if isinstance(rule, MissingRuleMetadata):
        return {}

    layers_by_name = {layer.name: layer for layer in rule.layers}
    materialized_layers = set(source.layers) | {rule.axis.x_layer}
    source_kind = "column" if rule.input_shape == "long" else "pattern"
    layers: JsonObject = {
        name: {
            "source": layer.source,
            "source_kind": source_kind,
        }
        for name, layer in layers_by_name.items()
        if name in materialized_layers
    }
    x_layer = layers_by_name[rule.axis.x_layer]
    x_mapping: JsonObject = {
        "layer": x_layer.name,
        "source": x_layer.source,
        "source_kind": source_kind,
    }
    return {
        "X": x_mapping,
        "layers": layers,
        "obs": _column_group_mapping(rule.columns.obs),
        "var": _column_group_mapping(rule.columns.var),
    }


def _column_group_mapping(group: ColumnGroup) -> JsonObject:
    """Map output column names to their vendor source or compute operation."""
    mapping: JsonObject = {}
    mapping.update(group.select)
    mapping.update({column.name: f"computed:{column.how}" for column in group.compute})
    return mapping


def _quantification_level_value(level: DescriptionQuantificationLevel) -> JsonValue:
    if isinstance(level, MissingQuantificationLevel):
        return None
    return level


def _software_name_value(software_name: DescriptionSoftwareName) -> JsonValue:
    if isinstance(software_name, MissingSoftwareName):
        return None
    return software_name
