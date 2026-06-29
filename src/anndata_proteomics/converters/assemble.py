"""Assemble ConversionPieces into an AnnData object with provenance metadata."""

from __future__ import annotations

import json
import math
from pathlib import Path

import anndata as ad
import pandas as pd

from anndata_proteomics.converters._pieces import ConversionPieces
from anndata_proteomics.modifications.pipeline import apply_modifications
from anndata_proteomics.params.anndata_io import write_search_parameters
from anndata_proteomics.params.registry import available_software, parse_params
from anndata_proteomics.rules.schema import ColumnCompute, ColumnGroup, ParseRule


def to_anndata(pieces: ConversionPieces, rule: ParseRule) -> ad.AnnData:
    """Build an AnnData from the converter pieces, recording the rule under `uns`.

    The full rule is stored as a JSON string under `uns['anndata_proteomics']['rule_json']`
    rather than a nested dict — h5py can't serialize the heterogeneous list-of-dicts
    structure of `layers`. Top-level fields are duplicated as plain strings for
    convenience (no need to parse JSON to read software_name, etc.).
    """
    adata = ad.AnnData(
        X=pieces.X,
        obs=pieces.obs,
        var=pieces.var,
        layers=pieces.layers,
    )
    adata.uns["anndata_proteomics"] = {
        "rule_json": json.dumps(rule.model_dump(mode="json", by_alias=True)),
        "schema_version": rule.schema_version,
        "software_name": rule.software_name,
        "input_shape": rule.input_shape,
        "quantification_level": rule.quantification_level,
    }
    if pieces.uns:
        for k, v in pieces.uns.items():
            adata.uns[k] = v
    return adata


def convert(
    df: pd.DataFrame,
    rule: ParseRule,
    *,
    params_path: str | Path | None = None,
) -> ad.AnnData:
    """One-shot: normalize modifications, dispatch long/wide, then assemble.

    Parameters
    ----------
    df
        Vendor-quant DataFrame (already loaded via ``readers``).
    rule
        Parsed TOML rule.
    params_path
        Optional vendor parameter file. When provided, the matching
        parameter parser (looked up by ``rule.software_name``) is invoked
        and its result is stored under
        ``uns['anndata_proteomics']['search_parameters']``.

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
        df = apply_modifications(df, rule.modifications)

    if rule.fragments is not None:
        # Fan packed per-precursor fragment lists out to one row per fragment before
        # the computed columns (ProForma_ion → ProForma_fragment) are materialized.
        # Drop columns the rule never reads first: explode multiplies the row count ~12x,
        # so carrying all ~60 vendor columns (mostly unused strings) through it is what
        # makes a full report cost many GB.
        from anndata_proteomics.converters._fragments import explode_fragments

        df = df[_columns_needed_for_long(df, rule)]
        df = explode_fragments(df, rule.fragments)

    df = _materialize_columns(df, rule)

    if rule.input_shape == "long":
        from anndata_proteomics.converters.long import convert_long

        pieces = convert_long(df, rule)
    else:
        from anndata_proteomics.converters.wide import convert_wide

        pieces = convert_wide(df, rule)

    adata = to_anndata(pieces, rule)
    if params_path is not None:
        _attach_search_parameters(adata, params_path, rule.software_name)
    return adata


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
    needed |= set(rule.columns.var.select.values())
    needed |= {layer.source for layer in rule.layers}
    if rule.modifications is not None:
        needed |= {rule.modifications.source_column, rule.modifications.output_column}
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


def _materialize_column_group(df: pd.DataFrame, group: ColumnGroup) -> None:
    for name, source in group.select.items():
        if source == "<sample>":
            continue
        if source not in df.columns:
            raise ValueError(f"cannot select column {name!r}; source {source!r} is missing")
        df[name] = df[source]
    for column in group.compute:
        df[column.name] = _compute_column(df, column)


def _compute_column(df: pd.DataFrame, column: ColumnCompute) -> pd.Series:
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


def _attach_search_parameters(adata: ad.AnnData, params_path: str | Path, software: str) -> None:
    """Parse ``params_path`` and store the result under ``uns``."""
    software_key = software.lower()
    available = {s.lower() for s in available_software()}
    if software_key not in available:
        # No parser yet for this vendor — keep the path as provenance, skip parsing.
        adata.uns["anndata_proteomics"]["search_parameters_path"] = str(params_path)
        return
    params = parse_params(params_path, software=software_key)
    write_search_parameters(adata, params, source_path=str(params_path))
