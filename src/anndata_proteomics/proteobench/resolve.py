"""Resolve ProteoBench roles against a converted APB object."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from anndata_proteomics.proteobench.config import ModuleSettings, ToolSettings
from anndata_proteomics.rules.schema import ParseRule


@dataclass(frozen=True)
class ResolvedRoles:
    """Concrete APB locations used by ProteoBench scoring."""

    proteins: str
    feature: str
    raw_file: str | None
    intensity: str = "X"
    reverse: str | None = None
    modification_source: str | None = None

    def as_dict(self) -> dict[str, str]:
        """Return a storage-safe ProteoBench-role mapping."""
        roles = {
            "Proteins": f"var:{self.proteins}",
            "feature": f"var:{self.feature}",
            "Intensity": self.intensity,
            "Raw file": f"obs:{self.raw_file}" if self.raw_file else "obs_names",
        }
        if self.reverse is not None:
            roles["Reverse"] = f"var:{self.reverse}"
        if self.modification_source is not None:
            roles["Modification source"] = f"var:{self.modification_source}"
        return roles


def resolve_target(obj: Any, level: str) -> Any:
    """Resolve the AnnData modality selected by a module level."""
    if hasattr(obj, "mod"):
        if level not in obj.mod:
            raise ValueError(
                f"MuData has no {level!r} modality; available modalities: {list(obj.mod)}"
            )
        return obj.mod[level]

    metadata = obj.uns.get("anndata_proteomics") or {}
    actual_level = metadata.get("quantification_level")
    if actual_level != level:
        raise ValueError(
            f"module level {level!r} does not match AnnData quantification level {actual_level!r}"
        )
    return obj


def resolve_roles(
    target: Any,
    module_settings: ModuleSettings,
    tool_settings: ToolSettings,
) -> tuple[ParseRule, ResolvedRoles]:
    """Resolve per-tool raw sources through the effective APB conversion rule."""
    metadata = target.uns.get("anndata_proteomics") or {}
    rule_json = metadata.get("rule_json")
    if not isinstance(rule_json, str):
        raise ValueError(
            "converted object has no string "
            "uns['anndata_proteomics']['rule_json']; rerun apb convert"
        )
    rule = ParseRule.model_validate_json(rule_json)
    if rule.quantification_level != module_settings.general.level:
        raise ValueError(
            f"module level {module_settings.general.level!r} does not match stored rule level "
            f"{rule.quantification_level!r}"
        )

    var_by_source = _invert_unique(rule.columns.var.select, axis="var")
    obs_by_source = _invert_unique(rule.columns.obs.select, axis="obs")

    protein_source = tool_settings.source_for("Proteins")
    assert protein_source is not None
    proteins = var_by_source.get(protein_source)
    if proteins is None or proteins not in target.var.columns:
        raise ValueError(
            f"ProteoBench maps raw column {protein_source!r} to 'Proteins', but the "
            "stored APB rule does not retain it in columns.var.select"
        )

    feature = {
        "ion": "ProForma_ion",
        "peptidoform": "ProForma_peptidoform",
    }[module_settings.general.level]
    if feature not in target.var.columns:
        raise ValueError(f"converted {module_settings.general.level} object lacks var[{feature!r}]")

    raw_file_source = tool_settings.source_for("Raw file")
    raw_file = obs_by_source.get(raw_file_source) if raw_file_source is not None else None
    if raw_file is not None and raw_file not in target.obs.columns:
        raise ValueError(f"stored rule resolves Raw file to missing obs column {raw_file!r}")

    reverse = None
    reverse_source = tool_settings.source_for("Reverse")
    if reverse_source is not None:
        reverse = var_by_source.get(reverse_source)
        if reverse is None or reverse not in target.var.columns:
            raise ValueError(
                f"per-tool settings require raw {reverse_source!r} for Reverse, but the "
                "stored APB rule does not retain it in columns.var.select"
            )

    modification_source = None
    modification_parser = tool_settings.modifications_parser
    if modification_parser is not None:
        raw_modification_source = (
            tool_settings.source_for(modification_parser.parse_column)
            or modification_parser.parse_column
        )
        modification_source = var_by_source.get(raw_modification_source)
        if modification_source is None or modification_source not in target.var.columns:
            raise ValueError(
                "per-tool settings parse modifications from raw column "
                f"{raw_modification_source!r}, but the stored APB rule does not retain it "
                "in columns.var.select"
            )

    _validate_intensity(rule, tool_settings)
    return rule, ResolvedRoles(
        proteins=proteins,
        feature=feature,
        raw_file=raw_file,
        reverse=reverse,
        modification_source=modification_source,
    )


def _invert_unique(mapping: dict[str, str], *, axis: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for output, source in mapping.items():
        if source in result:
            raise ValueError(
                f"stored APB rule maps raw {axis} source {source!r} to more than one output"
            )
        result[source] = output
    return result


def _validate_intensity(rule: ParseRule, tool_settings: ToolSettings) -> None:
    source = tool_settings.source_for("Intensity")
    if source is None:
        return
    x_layer = next(layer for layer in rule.layers if layer.name == rule.axis.x_layer)
    if rule.input_shape == "long" and x_layer.source != source:
        raise ValueError(
            f"ProteoBench maps {source!r} to Intensity, but APB X comes from {x_layer.source!r}"
        )
