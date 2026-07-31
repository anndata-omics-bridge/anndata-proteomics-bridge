"""Resolve ProteoBench roles against a converted APB object."""

from __future__ import annotations

from dataclasses import dataclass

from anndata import AnnData
from mudata import MuData

from anndata_proteomics.rules.anndata_io import read_stored_rule
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


def resolve_targets(obj: AnnData | MuData) -> list[AnnData]:
    """Return every quantitative AnnData a container holds.

    Scores are per level and belong in that level's own AnnData, so a MuData is scored modality by
    modality rather than through one modality selected by the module settings.
    """
    if isinstance(obj, MuData):
        if not obj.mod:
            raise ValueError("MuData has no modality to score")
        return [_require_anndata(name, target) for name, target in obj.mod.items()]
    return [obj]


def _require_anndata(name: str, target: object) -> AnnData:
    """Reject a modality that is not a quantitative AnnData."""
    if not isinstance(target, AnnData):
        raise TypeError(f"MuData modality {name!r} is not an AnnData")
    return target


def resolve_roles(target: AnnData) -> tuple[ParseRule, ResolvedRoles]:
    """Resolve canonical scoring locations from the stored APB conversion rule."""
    rule = read_stored_rule(target)
    if rule is None:
        raise ValueError(
            "converted object has no string "
            "uns['anndata_proteomics']['rule_json']; rerun apb convert"
        )
    proteins = rule.column_roles.protein_assignment
    if proteins is None:
        raise ValueError(
            "ProteoBench scoring requires column_roles.protein_assignment in the stored APB rule"
        )
    if proteins not in target.var.columns:
        raise ValueError(
            "stored APB rule maps column_roles.protein_assignment to missing "
            f"var column {proteins!r}"
        )

    return rule, ResolvedRoles(proteins=proteins)
