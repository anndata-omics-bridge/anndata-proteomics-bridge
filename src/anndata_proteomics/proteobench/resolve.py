"""Resolve ProteoBench roles against a converted APB object."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from anndata_proteomics.rules.schema import ParseRule

# The feature axis is the level's identity: `var_names` is built from the rule's `axis.var_keys`
# (one key as-is, several joined), so it is the ProteoBench feature at every quantification level.
FEATURE_AXIS = "var_names"


@dataclass(frozen=True)
class ResolvedRoles:
    """Concrete APB locations used by ProteoBench scoring."""

    proteins: str
    feature: str = FEATURE_AXIS
    intensity: str = "X"

    def as_dict(self) -> dict[str, str]:
        """Return a storage-safe ProteoBench-role mapping."""
        roles = {
            "Proteins": f"var:{self.proteins}",
            "feature": self.feature,
            "Intensity": self.intensity,
            "Sample name": "obs:sample_name",
            "Condition": "obs:condition",
        }
        return roles


def resolve_targets(obj: Any) -> list[Any]:
    """Return every quantitative AnnData a container holds.

    Scores are per level and belong in that level's own AnnData, so a MuData is scored modality by
    modality rather than through one modality selected by the module settings.
    """
    if hasattr(obj, "mod"):
        if not obj.mod:
            raise ValueError("MuData has no modality to score")
        return list(obj.mod.values())
    return [obj]


def resolve_roles(target: Any) -> tuple[ParseRule, ResolvedRoles]:
    """Resolve canonical scoring locations from the stored APB conversion rule."""
    metadata = target.uns.get("anndata_proteomics") or {}
    rule_json = metadata.get("rule_json")
    if not isinstance(rule_json, str):
        raise ValueError(
            "converted object has no string "
            "uns['anndata_proteomics']['rule_json']; rerun apb convert"
        )
    rule = ParseRule.model_validate_json(rule_json)
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

    return rule, ResolvedRoles(proteins=proteins)
