"""Shared axis-frame and key-index helpers for the long/wide converters."""

from __future__ import annotations

import pandas as pd

from anndata_proteomics.rules.schema import ParseRule

KEY_SEPARATOR = "_"
_SAMPLE_PLACEHOLDER = "<sample>"

# Added by `modifications.pipeline.apply_modifications` alongside the rule's own
# `output_column`; never present in the raw vendor input.
_MODIFICATION_OUTPUTS = frozenset({"stripped_sequence", "unknown_mod_tokens"})


def non_sample_columns(rule: ParseRule) -> frozenset[str]:
    """Every column a wide rule accounts for by name, so none can be a sample column.

    By the time the wide converter runs, ``apply_modifications`` and
    ``_materialize_columns`` have both been applied: the frame carries APB's derived
    columns *and* the rule's ``select`` outputs under their declared names, not the
    vendor's. Both spellings are excluded so a rule whose sample pattern cannot anchor
    on a suffix — AlphaDIA's run columns are bare run names — does not match them as
    extra samples.
    """
    out = set(_MODIFICATION_OUTPUTS)
    if rule.modifications is not None:
        out.add(rule.modifications.output_column)
        out.update(rule.modifications.source_columns)
    for group in (rule.columns.var, rule.columns.obs):
        out.update(group.select)
        out.update(group.optional_select)
        out.update(column.name for column in group.compute)
        out.update(group.select.values())
        out.update(group.optional_select.values())
    out.discard(_SAMPLE_PLACEHOLDER)
    return frozenset(out)


def join_keys(row: pd.Series) -> str:
    """Join a row of key values into a single string index token."""
    return KEY_SEPARATOR.join(str(v) for v in row)


def build_index(df: pd.DataFrame, keys: list[str]) -> pd.Series:
    """Build a string index from one or more key columns.

    Vectorised string concatenation (not a row-wise apply) so it stays cheap on the
    full, un-deduplicated frame the long converter scatters from.
    """
    if len(keys) == 1:
        return df[keys[0]].astype("string")
    joined = df[keys[0]].astype("string")
    for key in keys[1:]:
        joined = joined + KEY_SEPARATOR + df[key].astype("string")
    return joined


def build_axis_frame(df: pd.DataFrame, keys: list[str], output_columns: list[str]) -> pd.DataFrame:
    """Take first occurrence per key tuple for already-materialized output columns.

    ``output_columns`` is the rule's *declared* column set, which may name
    ``optional_select`` entries this export does not carry. Those are the only declared
    columns that can be absent here — ``_materialize_column_group`` raises on a missing
    required ``select`` source, and every compute assigns its column — so filtering to the
    present ones drops exactly the skipped optional columns and nothing else.
    """
    present = [column for column in output_columns if column in df.columns]
    needed_cols = list(dict.fromkeys(list(keys) + present))
    block = df[needed_cols].drop_duplicates(subset=keys).copy()
    out = block[present].copy()
    out.index = pd.Index(build_index(block, keys), name=KEY_SEPARATOR.join(keys))
    return out
