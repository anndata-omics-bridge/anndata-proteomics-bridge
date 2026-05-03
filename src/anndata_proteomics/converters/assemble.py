"""Assemble ConversionPieces into an AnnData object with provenance metadata."""

from __future__ import annotations

import anndata as ad

from anndata_proteomics.converters._pieces import ConversionPieces
from anndata_proteomics.rules.schema import ParseRule


def to_anndata(pieces: ConversionPieces, rule: ParseRule) -> ad.AnnData:
    """Build an AnnData from the converter pieces, recording the rule under `uns`."""
    adata = ad.AnnData(
        X=pieces.X,
        obs=pieces.obs,
        var=pieces.var,
        layers=pieces.layers,
    )
    adata.uns["anndata_proteomics"] = {
        "rule": rule.model_dump(mode="json"),
        "schema_version": rule.schema_version,
        "software_name": rule.software_name,
        "input_shape": rule.input_shape,
        "quantification_level": rule.quantification_level,
    }
    if pieces.uns:
        for k, v in pieces.uns.items():
            adata.uns[k] = v
    return adata


def convert(df, rule: ParseRule) -> ad.AnnData:
    """One-shot: dispatch to long/wide based on rule.input_shape, then assemble."""
    if rule.input_shape == "long":
        from anndata_proteomics.converters.long import convert_long

        pieces = convert_long(df, rule)
    else:
        from anndata_proteomics.converters.wide import convert_wide

        pieces = convert_wide(df, rule)
    return to_anndata(pieces, rule)
