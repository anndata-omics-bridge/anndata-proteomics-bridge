"""Resolve ProteoBench roles against a converted APB object."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from anndata_proteomics.proteobench.config import ModuleSettings
from anndata_proteomics.rules.schema import ParseRule


@dataclass(frozen=True)
class ResolvedRoles:
    """Concrete APB locations used by ProteoBench scoring."""

    proteins: str
    feature: str
    raw_file: str | None
    intensity: str = "X"

    def as_dict(self) -> dict[str, str]:
        """Return a storage-safe ProteoBench-role mapping."""
        roles = {
            "Proteins": f"var:{self.proteins}",
            "feature": f"var:{self.feature}",
            "Intensity": self.intensity,
            "Raw file": f"obs:{self.raw_file}" if self.raw_file else "obs_names",
        }
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
) -> tuple[ParseRule, ResolvedRoles]:
    """Resolve canonical scoring locations from the stored APB conversion rule."""
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

    if rule.column_roles is None:
        raise ValueError(
            "stored APB rule has no column_roles; rerun conversion with a rule "
            "that declares protein_accessions"
        )
    proteins = rule.column_roles.protein_accessions
    if proteins not in target.var.columns:
        raise ValueError(
            f"stored APB rule maps protein accessions to missing var column {proteins!r}"
        )

    feature = {
        "ion": "ProForma_ion",
        "peptidoform": "ProForma_peptidoform",
    }[module_settings.general.level]
    if feature not in target.var.columns:
        raise ValueError(f"converted {module_settings.general.level} object lacks var[{feature!r}]")

    raw_file = next(
        (column for column in rule.axis.obs_keys if column in target.obs.columns),
        None,
    )
    return rule, ResolvedRoles(
        proteins=proteins,
        feature=feature,
        raw_file=raw_file,
    )
