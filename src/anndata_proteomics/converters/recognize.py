"""Match a DataFrame's headers to one of the packaged ParseRules."""

from __future__ import annotations

import re
from collections.abc import Iterable

from anndata_proteomics.rules.loader import load_rule
from anndata_proteomics.rules.registry import iter_packaged_rules
from anndata_proteomics.rules.schema import ParseRule

_SAMPLE_PLACEHOLDER = "<sample>"


def _expected_long_columns(rule: ParseRule) -> set[str]:
    """Vendor columns a long rule expects to see in the input."""
    out = set(rule.columns.obs.values()) | set(rule.columns.var.values())
    out.update(layer.source_column for layer in rule.layers if layer.source_column)
    out.discard(_SAMPLE_PLACEHOLDER)
    return out


def _required_var_columns(rule: ParseRule) -> set[str]:
    """Vendor columns a wide rule expects on the var axis (per-feature, not per-sample)."""
    return {v for v in rule.columns.var.values() if v != _SAMPLE_PLACEHOLDER}


def matches(headers: Iterable[str], rule: ParseRule) -> bool:
    """Does the given header set plausibly match this rule?

    Long: every referenced vendor column must be present.
    Wide: every layer's `column_pattern` must match at least one header, and
    the var-side vendor columns must be present.
    """
    headers_set = set(headers)
    if rule.input_shape == "long":
        return _expected_long_columns(rule).issubset(headers_set)
    # wide
    for layer in rule.layers:
        if layer.column_pattern is None:
            return False
        pattern = re.compile(layer.column_pattern)
        if not any(pattern.match(h) for h in headers_set):
            return False
    return _required_var_columns(rule).issubset(headers_set)


def recognize(headers: Iterable[str]) -> ParseRule | None:
    """Find the unique packaged ParseRule that matches the headers.

    Returns None if zero rules match or multiple match (caller must specify
    the rule explicitly in that case).
    """
    headers_set = set(headers)
    candidates = [load_rule(p) for p in iter_packaged_rules()]
    hits = [r for r in candidates if matches(headers_set, r)]
    return hits[0] if len(hits) == 1 else None
