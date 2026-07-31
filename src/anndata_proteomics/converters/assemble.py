"""Convert vendor tables into backend-neutral quantitative pieces."""

from __future__ import annotations

import pandas as pd

from anndata_proteomics.converters._fragments import explode_fragments
from anndata_proteomics.converters._pieces import ConversionPieces
from anndata_proteomics.converters.axis_types import AxisColumnContext, AxisName, coerce_axis_column
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
    _materialize_column_group(out, rule.columns.obs, "obs")
    _materialize_column_group(out, rule.columns.var, "var")
    return out


def _materialize_column_group(
    df: pd.DataFrame,
    group: ColumnGroup,
    axis: AxisName,
) -> frozenset[str]:
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
        df[name] = _coerce_selected_column(df[source], group, axis, name, source)
    skipped: set[str] = set()
    for name, source in group.optional_select.items():
        if source in df.columns:
            df[name] = _coerce_selected_column(df[source], group, axis, name, source)
        else:
            skipped.add(name)
    for column in group.compute:
        df[column.name] = _compute_column(
            df,
            column,
            allow_missing=frozenset(skipped),
        ).astype("string")
    return frozenset(skipped)


def _coerce_selected_column(
    values: pd.Series,
    group: ColumnGroup,
    axis: AxisName,
    output_name: str,
    source_name: str,
) -> pd.Series:
    """Coerce one selected source through its rule-owned logical contract."""
    return coerce_axis_column(
        values,
        AxisColumnContext(
            axis=axis,
            output_name=output_name,
            source_name=source_name,
            logical_type=group.type_for(output_name),
        ),
    )


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
        return _compute_proforma_ion(df, column)
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


def _compute_proforma_ion(df: pd.DataFrame, column: ColumnCompute) -> pd.Series:
    """Combine a string peptidoform and already-typed positive integer charge."""
    sequence_key, charge_key = column.from_
    missing_sources = [key for key in (sequence_key, charge_key) if key not in df.columns]
    if missing_sources:
        raise ValueError(
            f"cannot compute column {column.name!r}; source column(s) missing: {missing_sources}"
        )
    charges = df[charge_key]
    if charges.isna().any():
        raise ValueError("cannot derive proforma_ion from missing charge")
    nonpositive = charges.le(0)
    if nonpositive.any():
        examples = charges.loc[nonpositive].drop_duplicates().head(5).tolist()
        raise ValueError(f"charge must be positive; examples={examples}")
    sequences = df[sequence_key].astype("string")
    return sequences + "/" + charges.astype("string")
