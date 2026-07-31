"""AnnData persistence adapter for effective APB parsing rules."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from anndata import AnnData

from anndata_proteomics.adapters.anndata.namespace import (
    MissingNamespaceText,
    has_namespace_key,
    read_namespace_text,
)
from anndata_proteomics.rules.schema import ParseRule, QuantificationLevel

_RULE_KEY = "rule_json"
type ColumnRoleName = Literal["protein_assignment", "fasta_accessions"]


@dataclass(frozen=True, slots=True)
class MissingStoredRule:
    """Conversion stored no effective parsing rule on this container."""


MISSING_STORED_RULE = MissingStoredRule()


def has_stored_rule(target: AnnData) -> bool:
    """Return whether conversion stored an effective parsing rule."""
    return has_namespace_key(target, _RULE_KEY)


def read_stored_rule(target: AnnData) -> ParseRule | MissingStoredRule:
    """Read the effective parsing rule, or report its absence, with one validation."""
    payload = read_namespace_text(target, _RULE_KEY)
    if isinstance(payload, MissingNamespaceText):
        return MISSING_STORED_RULE
    return ParseRule.model_validate_json(payload)


def require_stored_rule(target: AnnData) -> ParseRule:
    """Read and validate the effective parsing rule stored by conversion."""
    rule = read_stored_rule(target)
    if isinstance(rule, MissingStoredRule):
        raise ValueError(
            "converted object has no string "
            "uns['anndata_proteomics']['rule_json']; rerun apb convert"
        )
    return rule


def require_quantification_level(target: AnnData) -> QuantificationLevel:
    """Return the validated quantification level stored in the effective rule."""
    return require_stored_rule(target).quantification_level


def has_stored_column_role(target: AnnData, role: ColumnRoleName) -> bool:
    """Return whether the stored rule declares one semantic ``var`` role."""
    rule = read_stored_rule(target)
    if isinstance(rule, MissingStoredRule):
        return False
    return getattr(rule.column_roles, role) is not None


def require_stored_column_role(target: AnnData, role: ColumnRoleName) -> str:
    """Return and validate the ``var`` column declared for one semantic role."""
    return column_role(require_stored_rule(target), target, role)


def column_role(rule: ParseRule, target: AnnData, role: ColumnRoleName) -> str:
    """Resolve one semantic ``var`` role against a rule the caller already holds."""
    column = getattr(rule.column_roles, role)
    if column is None:
        raise ValueError(f"stored APB rule does not declare column_roles.{role}")
    if column not in target.var.columns:
        raise ValueError(
            f"stored APB rule maps column_roles.{role} to missing var column {column!r}"
        )
    return column
