"""Read the effective APB parsing rule stored on a converted container."""

from __future__ import annotations

from typing import Literal

from anndata import AnnData

from anndata_proteomics._containers import UnsHolder
from anndata_proteomics.rules.schema import ParseRule

_NAMESPACE = "anndata_proteomics"
_RULE_KEY = "rule_json"
type ColumnRoleName = Literal["protein_assignment", "fasta_accessions"]


def read_stored_rule(target: UnsHolder) -> ParseRule | None:
    """Return the effective parsing rule stored by conversion, if present."""
    namespace = target.uns.get(_NAMESPACE) or {}
    raw = namespace.get(_RULE_KEY)
    if raw is None:
        return None
    if not isinstance(raw, str):
        raise ValueError(
            "converted object has no string "
            f"uns[{_NAMESPACE!r}][{_RULE_KEY!r}]; got {type(raw).__name__}"
        )
    return ParseRule.model_validate_json(raw)


def read_stored_column_role(target: AnnData, role: ColumnRoleName) -> str | None:
    """Return one declared semantic ``var`` column without guessing its name."""
    rule = read_stored_rule(target)
    if rule is None:
        return None
    column = getattr(rule.column_roles, role)
    if column is None:
        return None
    if column not in target.var.columns:
        raise ValueError(
            f"stored APB rule maps column_roles.{role} to missing var column {column!r}"
        )
    return column
