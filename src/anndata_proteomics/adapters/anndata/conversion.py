"""AnnData and MuData persistence for typed conversion results."""

from __future__ import annotations

import json
from collections.abc import Mapping

import mudata
from anndata import AnnData
from mudata import MuData

from anndata_proteomics.adapters.anndata.namespace import update_namespace, write_namespace
from anndata_proteomics.adapters.anndata.params import write_search_parameters
from anndata_proteomics.converters.pipeline import (
    LEVELS,
    ParameterResolution,
    RuleSelectionMethod,
    version_status,
)
from anndata_proteomics.rules.schema import QuantificationLevel
from anndata_proteomics.workflows.conversion import LevelConversion

_PREFIX: Mapping[QuantificationLevel, str] = {
    "fragment": "frg:",
    "ion": "ion:",
    "peptidoform": "pfm:",
    "peptide": "pep:",
    "protein": "prt:",
}


def to_anndata(conversion: LevelConversion) -> AnnData:
    """Persist one backend-neutral level conversion as AnnData."""
    pieces = conversion.pieces
    rule = conversion.selection.rule
    adata = AnnData(
        X=pieces.X,
        obs=pieces.obs,
        var=pieces.var,
        layers=pieces.layers,
    )
    write_namespace(
        adata,
        {
            "rule_json": json.dumps(rule.model_dump(mode="json", by_alias=True)),
            "schema_version": rule.schema_version,
            "software_name": rule.software_name,
            "input_shape": rule.input_shape,
            "quantification_level": rule.quantification_level,
            "rule_selection_method": conversion.selection.method,
        },
    )
    return adata


def to_mudata(
    conversions: Mapping[QuantificationLevel, LevelConversion],
) -> MuData:
    """Persist selected level conversions as one shared-observation MuData."""
    if not conversions:
        raise ValueError("no level conversions supplied")
    modalities: dict[str, AnnData] = {}
    methods: set[RuleSelectionMethod] = set()
    for level in LEVELS:
        conversion = conversions.get(level)
        if conversion is None:
            continue
        if conversion.level != level:
            raise ValueError(
                f"conversion key {level!r} does not match result level {conversion.level!r}"
            )
        adata = to_anndata(conversion)
        adata.var_names = [_PREFIX[level] + str(name) for name in adata.var_names]
        modalities[level] = adata
        methods.add(conversion.selection.method)
    if len(methods) != 1:
        raise ValueError(f"MuData conversion requires one rule-selection method, got {methods}")
    with mudata.set_options(pull_on_update=False):
        result = MuData(modalities, axis=0)
    update_namespace(result, {"rule_selection_method": methods.pop()})
    return result


def write_parameter_resolution(
    target: AnnData | MuData,
    resolution: ParameterResolution,
) -> None:
    """Persist one already-parsed parameter result on one storage target."""
    update_namespace(
        target,
        {
            "search_parameters_version_status": version_status(resolution.version),
            "search_parameters_path": str(resolution.source_path),
        },
    )
    write_search_parameters(target, resolution.parameters)
