"""Convert vendor tables into backend-neutral quantitative pieces."""

from __future__ import annotations

import math

import pandas as pd

from anndata_proteomics.converters._fragments import explode_fragments
from anndata_proteomics.converters._pieces import ConversionPieces
from anndata_proteomics.converters.checks import check_layer_occupancy
from anndata_proteomics.converters.long import convert_long
from anndata_proteomics.converters.wide import convert_wide
from anndata_proteomics.modifications.pipeline import apply_modifications
from anndata_proteomics.rules.schema import ColumnCompute, ColumnGroup, ParseRule


def convert_table(
    df: pd.DataFrame,
    rule: ParseRule,
    *,
    strict: bool = False,
) -> ConversionPieces:
    """Normalize and convert one vendor table into backend-neutral pieces.

    Parameters
    ----------
    df
        Vendor-quant DataFrame (already loaded via ``readers``).
    rule
        Parsed JSON rule.
    strict
        Promote non-``X`` layer-contract warnings to errors.
    Notes
    -----
    If ``rule.modifications`` is set *and* a ``proforma_sequence`` /
    ``stripped_sequence`` compute consumes it, modification normalization runs
    before the long/wide dispatch so those computes can read the normalized
    columns. A rule that only *inherits* a ``[modifications]`` block from its
    vendor base without a consuming compute (e.g. protein levels) skips it — the
    output is identical and the per-row tokenization cost is avoided.
    """
    if rule.modifications is not None and _modifications_consumed(rule):
        # apply_modifications adds its output columns in place, so copy first: callers
        # converting several levels from one table must not see another level's columns.
        df = apply_modifications(df.copy(), rule.modifications)

    if rule.fragments is not None:
        # Fan packed per-precursor fragment lists out to one row per fragment before
        # the computed columns (ProForma_ion → ProForma_fragment) are materialized.
        # Drop columns the rule never reads first: explode multiplies the row count ~12x,
        # so carrying all ~60 vendor columns (mostly unused strings) through it is what
        # makes a full report cost many GB.
        df = df[_columns_needed_for_long(df, rule)]
        df = explode_fragments(df, rule.fragments)

    df = _materialize_columns(df, rule)

    if rule.input_shape == "long":
        pieces = convert_long(df, rule)
    else:
        pieces = convert_wide(df, rule)

    check_layer_occupancy(pieces.layers, x_layer=rule.axis.x_layer, strict=strict)
    return pieces


def _modifications_consumed(rule: ParseRule) -> bool:
    """True iff a var compute reads the modification output (proforma/stripped sequence).

    ``apply_modifications`` only adds the ``proforma_sequence`` / ``stripped_sequence`` columns;
    the sole consumers are ``how="proforma_sequence"`` / ``"stripped_sequence"`` computes. When a
    rule merely inherits a ``[modifications]`` block from its vendor base but has no such compute
    (protein levels), applying it would add unused columns — skip it.
    """
    return any(
        column.how in {"proforma_sequence", "stripped_sequence"}
        for column in rule.columns.var.compute
    )


def _columns_needed_for_long(df: pd.DataFrame, rule: ParseRule) -> list[str]:
    """Vendor/derived columns a long rule reads downstream of modifications.

    Used to trim the frame before the fragment explode. Keeps select sources, layer
    sources, the fragment packed columns, and the modification-derived columns that
    later computes (ProForma_*) consume — everything else is dead weight once exploded.
    """
    needed: set[str] = set(rule.columns.obs.select.values())
    needed |= set(rule.columns.obs.optional_select.values())
    needed |= set(rule.columns.var.select.values())
    needed |= set(rule.columns.var.optional_select.values())
    needed |= {layer.source for layer in rule.layers}
    if rule.modifications is not None:
        needed |= {*rule.modifications.source_columns, rule.modifications.output_column}
        needed.add("stripped_sequence")
    if rule.fragments is not None:
        if rule.fragments.label_strategy == "column":
            needed.add(rule.fragments.label_column)
        needed |= set(rule.fragments.value_columns)
    needed.discard("<sample>")
    # preserve original column order; keep only columns actually present
    return [column for column in df.columns if column in needed]


def _materialize_columns(df: pd.DataFrame, rule: ParseRule) -> pd.DataFrame:
    """Materialize declared selected and computed columns on a working DataFrame."""
    out = df.copy()
    _materialize_column_group(out, rule.columns.obs)
    _materialize_column_group(out, rule.columns.var)
    return out


def _materialize_column_group(df: pd.DataFrame, group: ColumnGroup) -> frozenset[str]:
    """Materialize one group's columns; return the optional names whose source was absent.

    A ``select`` source must be present. An ``optional_select`` source that this export does
    not carry is skipped, and its name is reported so ``_compute_column`` can tell a
    legitimately absent optional source from a broken rule.
    """
    for name, source in group.select.items():
        if source == "<sample>":
            continue
        if source not in df.columns:
            raise ValueError(f"cannot select column {name!r}; source {source!r} is missing")
        df[name] = df[source]
    skipped: set[str] = set()
    for name, source in group.optional_select.items():
        if source in df.columns:
            df[name] = df[source]
        else:
            skipped.add(name)
    for column in group.compute:
        df[column.name] = _compute_column(df, column, allow_missing=frozenset(skipped))
    return frozenset(skipped)


def _compute_generic_column(
    df: pd.DataFrame,
    column: ColumnCompute,
    allow_missing: frozenset[str],
) -> pd.Series:
    """Compute a ``coalesce`` / ``join_nonempty`` column over its present sources.

    Sources named in ``allow_missing`` are skipped ``optional_select`` columns and drop out
    of the chain. Anything else missing is a broken rule and still raises.
    """
    sources = [key for key in column.from_ if key not in allow_missing]
    missing = [key for key in sources if key not in df.columns]
    if missing:
        raise ValueError(
            f"cannot compute column {column.name!r}; source column(s) missing: {missing}"
        )
    if not sources:
        raise ValueError(
            f"cannot compute column {column.name!r}; every source column is an "
            f"optional_select absent from this input: {list(column.from_)}"
        )
    if column.how == "coalesce":
        return _coalesce_columns(df, sources)
    if column.separator is None:
        raise ValueError(
            f"cannot compute column {column.name!r}; how={column.how!r} requires a separator"
        )
    return _join_nonempty_columns(df, sources, column.separator)


def _compute_column(
    df: pd.DataFrame,
    column: ColumnCompute,
    *,
    allow_missing: frozenset[str] = frozenset(),
) -> pd.Series:
    if column.how in {"coalesce", "join_nonempty"}:
        return _compute_generic_column(df, column, allow_missing)
    if column.how in {"proforma_sequence", "stripped_sequence"}:
        source_key = column.how
        if source_key not in df.columns:
            raise ValueError(
                f"cannot compute column {column.name!r}; APB column {source_key!r} is missing"
            )
        return df[source_key]
    if column.how == "proforma_ion":
        sequence_key, charge_key = column.from_
        missing = [key for key in (sequence_key, charge_key) if key not in df.columns]
        if missing:
            raise ValueError(
                f"cannot compute column {column.name!r}; source column(s) missing: {missing}"
            )
        return pd.Series(
            [
                f"{sequence}/{_format_charge(charge)}"
                for sequence, charge in zip(df[sequence_key], df[charge_key], strict=True)
            ],
            index=df.index,
        )
    if column.how == "proforma_fragment":
        ion_key, label_key = column.from_
        missing = [key for key in (ion_key, label_key) if key not in df.columns]
        if missing:
            raise ValueError(
                f"cannot compute column {column.name!r}; source column(s) missing: {missing}"
            )
        return pd.Series(
            [f"{ion}/{label}" for ion, label in zip(df[ion_key], df[label_key], strict=True)],
            index=df.index,
        )
    raise ValueError(f"unsupported column compute mode: {column.how!r}")


def _coalesce_columns(df: pd.DataFrame, sources: list[str]) -> pd.Series:
    """Return the first non-null source value in declaration order."""
    result = df[sources[0]].astype(object).copy()
    for source in sources[1:]:
        result = result.where(result.notna(), df[source].astype(object))
    return result


def _join_nonempty_columns(
    df: pd.DataFrame,
    sources: list[str],
    separator: str,
) -> pd.Series:
    """Join non-null, non-empty source values in declaration order."""
    result = pd.Series("", index=df.index, dtype="string")
    has_value = pd.Series(False, index=df.index)
    for source in sources:
        values = df[source].astype("string")
        valid = values.notna() & values.ne("")
        append = has_value & valid
        result = result.mask(append, result + separator + values)
        first = ~has_value & valid
        result = result.mask(first, values)
        has_value |= valid
    return result.mask(~has_value, pd.NA)


def _format_charge(value: object) -> str:
    """Normalize charge values for ProForma ion identifiers."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        raise ValueError("cannot derive proforma_ion from missing charge")
    if isinstance(value, (int, float)):
        numeric = float(value)
    else:
        text = str(value).strip()
        if not text:
            raise ValueError("cannot derive proforma_ion from empty charge")
        try:
            numeric = float(text)
        except ValueError as exc:
            raise ValueError(f"charge must be numeric, got {value!r}") from exc
    if not numeric.is_integer():
        raise ValueError(f"charge must be an integer value, got {value!r}")
    charge = int(numeric)
    if charge <= 0:
        raise ValueError(f"charge must be positive, got {value!r}")
    return str(charge)
