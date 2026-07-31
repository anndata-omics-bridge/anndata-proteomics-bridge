"""Parse APB namespace values into exact description metadata."""

from __future__ import annotations

import json

from anndata_proteomics.description import (
    DescriptionConversionMetadata,
    DescriptionMetadata,
    MissingProteoBenchMetadata,
    MissingQcMetadata,
    MissingQuantificationLevel,
    MissingRuleMetadata,
    MissingSearchParameters,
    MissingSoftwareName,
)
from anndata_proteomics.params.model import Parameters
from anndata_proteomics.rules.schema import ParseRule
from anndata_proteomics.serialization import JsonObject, JsonValue, to_json_compatible

_CONVERSION_FIELDS = (
    "schema_version",
    "quantification_level",
    "software_name",
    "rule_selection_method",
    "search_parameters_path",
    "search_parameters_version_status",
)


def parse_description_metadata(namespace: JsonObject) -> DescriptionMetadata:
    """Validate and narrow one extracted APB namespace for description calculation."""
    conversion = DescriptionConversionMetadata.model_validate(
        {field: namespace[field] for field in _CONVERSION_FIELDS if field in namespace}
    )
    qc = _read_qc(namespace["qc"]) if "qc" in namespace else MissingQcMetadata()
    return DescriptionMetadata(
        quantification_level=(
            conversion.quantification_level
            if conversion.quantification_level is not None
            else MissingQuantificationLevel()
        ),
        software_name=(
            conversion.software_name
            if conversion.software_name is not None
            else MissingSoftwareName()
        ),
        conversion=conversion,
        search_parameters=_read_search_parameters(namespace.get("search_parameters")),
        rule=_read_rule(namespace.get("rule_json")),
        annotations=_annotation_provenance(
            namespace.get("obs_annotations_json"),
            namespace.get("var_annotations_json"),
            namespace.get("fasta_config"),
        ),
        qc=qc,
        proteobench=_read_proteobench(namespace.get("proteobench")),
    )


def _read_search_parameters(
    raw: JsonValue,
) -> Parameters | MissingSearchParameters:
    if raw is None:
        return MissingSearchParameters()
    if not isinstance(raw, str):
        raise TypeError("stored search parameters must be JSON text")
    return Parameters.model_validate(json.loads(raw))


def _read_rule(raw: JsonValue) -> ParseRule | MissingRuleMetadata:
    if raw is None:
        return MissingRuleMetadata()
    if not isinstance(raw, str):
        raise ValueError(
            "converted object has no string "
            "uns['anndata_proteomics']['rule_json']; "
            f"got {type(raw).__name__}"
        )
    return ParseRule.model_validate_json(raw)


def _annotation_provenance(
    obs_payload: JsonValue,
    var_payload: JsonValue,
    fasta_config_payload: JsonValue,
) -> JsonObject:
    result: JsonObject = {}
    obs = _decode_json_value(obs_payload)
    if obs:
        result["obs"] = obs
    var = _decode_json_value(var_payload)
    if var:
        result["var"] = var
    fasta_config = _decode_json_value(fasta_config_payload)
    if fasta_config is not None:
        result["fasta_config"] = fasta_config
    return result


def _read_qc(raw: JsonValue) -> JsonValue:
    return _decode_json_value(raw)


def _read_proteobench(
    raw: JsonValue,
) -> JsonObject | MissingProteoBenchMetadata:
    if not isinstance(raw, dict) or "scores" not in raw:
        return MissingProteoBenchMetadata()
    return raw


def _decode_json_value(value: JsonValue) -> JsonValue:
    if value is None:
        return None
    if isinstance(value, str):
        return to_json_compatible(json.loads(value))
    return to_json_compatible(value)
