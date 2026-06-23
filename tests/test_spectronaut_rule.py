"""Spectronaut TOML semantics: report-backed levels and correct layer placement."""

from __future__ import annotations

import pandas as pd
import pytest

from anndata_proteomics.converters.assemble import convert
from anndata_proteomics.converters.recognize import matches
from anndata_proteomics.readers.dispatch import read_table
from anndata_proteomics.rules.loader import load_packaged_rule
from anndata_proteomics.rules.registry import resolve_rule_path
from anndata_proteomics.rules.schema import ParseRule
from anndata_proteomics.scripts import _ui_support as ui


def _spectronaut_catalog() -> pd.DataFrame:
    catalog = ui.load_catalog()
    return catalog[catalog["slug"] == "spectronaut"].reset_index(drop=True)


def _resolve_var_feature_sources(rule: ParseRule) -> set[str]:
    """Vendor columns that define a feature's identity, by walking axis.var_keys back through
    columns.var.compute to the selected vendor columns they are built from."""
    select = rule.columns.var.select
    compute = {column.name: column.from_ for column in rule.columns.var.compute}

    def resolve(name: str) -> set[str]:
        if name in select:
            return {select[name]}
        if name in compute:
            return set().union(*(resolve(source) for source in compute[name]))
        return set()

    return set().union(*(resolve(key) for key in rule.axis.var_keys))


def _sample_sources(rule: ParseRule) -> set[str]:
    return {rule.columns.obs.select[key] for key in rule.axis.obs_keys}


def _feature_invariance(df: pd.DataFrame, rule: ParseRule) -> tuple[pd.Series, pd.Series]:
    """For every declared vendor source column, how feature-bound is it?

    Returns two Series indexed by vendor column:
      - frac_constant: among features seen in >=2 samples, the fraction whose value is identical
        across those samples (1.0 => the value is a property of the feature, not the measurement).
      - global_distinct: number of distinct values across the whole table (1 => degenerate here).
    """
    feature_cols = sorted(_resolve_var_feature_sources(rule))
    sample_cols = sorted(_sample_sources(rule))
    feature = df[feature_cols].astype(str).agg("\x1f".join, axis=1)
    sample = df[sample_cols].astype(str).agg("\x1f".join, axis=1)

    samples_per_feature = sample.groupby(feature).nunique()
    multi_sample = set(samples_per_feature.index[samples_per_feature >= 2])
    mask = feature.isin(multi_sample)

    cols = sorted(
        {layer.source for layer in rule.layers if layer.source in df.columns}
        | {src for src in rule.columns.var.select.values() if src in df.columns}
    )
    distinct_per_feature = df.loc[mask, cols].groupby(feature[mask]).nunique(dropna=False)
    frac_constant = (distinct_per_feature <= 1).mean()
    global_distinct = df[cols].nunique(dropna=False)
    return frac_constant, global_distinct


def test_spectronaut_has_report_backed_ion_protein_and_fragment_rules() -> None:
    assert resolve_rule_path("spectronaut", "ion") is not None
    assert resolve_rule_path("spectronaut", "protein") is not None
    assert resolve_rule_path("spectronaut", "fragment") is not None
    assert resolve_rule_path("spectronaut", "peptidoform") is None
    assert resolve_rule_path("spectronaut", "peptide") is None


def test_spectronaut_rule_matches_cached_common_headers() -> None:
    rules = [
        load_packaged_rule("spectronaut", "ion"),
        load_packaged_rule("spectronaut", "protein"),
    ]
    for _, row in _spectronaut_catalog().iterrows():
        headers = pd.read_csv(ui._dataset_path(row["input_file_path"]), sep="\t", nrows=0).columns
        for rule in rules:
            assert matches(headers, rule), f"{rule.quantification_level}: {row['input_file_path']}"
        assert not matches(headers, load_packaged_rule("spectronaut", "fragment"))


def test_spectronaut_catalog_offers_mudata() -> None:
    catalog = _spectronaut_catalog()
    if catalog.empty:
        return
    assert (
        catalog["targets"].apply(lambda targets: {"ion", "protein", "mudata"} <= set(targets)).all()
    )
    assert not catalog["targets"].apply(lambda targets: "fragment" in set(targets)).any()


@pytest.mark.parametrize("level", ["ion", "protein"])
def test_spectronaut_conversion_matches_declared_columns(level: str) -> None:
    """Conformance: the converted AnnData realizes exactly what the TOML declares — every
    selected/computed var and obs column is present, and the layer set is precisely the declared
    one (no stray, no missing). The expectations are read from the rule, not hardcoded."""
    catalog = _spectronaut_catalog()
    if catalog.empty:
        return
    rule = load_packaged_rule("spectronaut", level)
    df = read_table(ui._dataset_path(catalog.iloc[0]["input_file_path"]))
    run = df["R.FileName"].iloc[0]
    subset = df[df["R.FileName"] == run].head(2000).copy()

    adata = convert(subset, rule)

    assert set(rule.columns.var.names) <= set(adata.var.columns)
    assert set(rule.columns.obs.names) <= set(adata.obs.columns)
    assert {layer.name for layer in rule.layers} == set(adata.layers.keys())

    # when the feature key is a directly selected vendor column, conversion must collapse to one
    # var per distinct key value (computed keys such as ProForma_ion are exercised at ion level).
    if set(rule.axis.var_keys) <= set(rule.columns.var.select):
        key_source = rule.columns.var.select[rule.axis.var_keys[0]]
        assert adata.n_vars == subset[key_source].nunique()


@pytest.mark.parametrize("level", ["ion", "protein"])
def test_spectronaut_declared_layout_matches_data(level: str) -> None:
    """Placement correctness: a layer is an obs x var matrix, so a column belongs in [[layers]]
    only if it varies across samples for the same feature. This checks the TOML's split against
    the real multi-sample report (the dimensionality test), which is what conformance alone cannot
    do — deriving the truth from the same TOML would let a misplacement pass unnoticed."""
    catalog = _spectronaut_catalog()
    if catalog.empty:
        return
    rule = load_packaged_rule("spectronaut", level)
    df = read_table(ui._dataset_path(catalog.iloc[0]["input_file_path"]))
    if df[list(_sample_sources(rule))].drop_duplicates().shape[0] < 2:
        pytest.skip("single-sample export cannot distinguish layers from .var")

    frac_constant, global_distinct = _feature_invariance(df, rule)
    invariant = 0.999  # feature-bound columns are constant for ~every multi-sample feature

    # A selected .var column must be a property of the feature: identical across samples.
    for name, source in rule.columns.var.select.items():
        if source not in frac_constant.index:
            continue
        assert frac_constant[source] >= invariant, (
            f"{level}: var column {name!r} ({source}) differs across samples for the same feature "
            f"(constant for only {frac_constant[source]:.1%} of features) — it should be a layer"
        )

    # An informative layer must carry per-sample signal. A layer that is constant per feature yet
    # has more than one distinct value is really a .var attribute. Columns with a single global
    # value are degenerate in this export (e.g. all False / all NaN) and stay layers by design.
    for layer in rule.layers:
        source = layer.source
        if source not in frac_constant.index or global_distinct[source] <= 1:
            continue
        assert frac_constant[source] < invariant, (
            f"{level}: layer {layer.name!r} ({source}) is constant per feature "
            f"({frac_constant[source]:.1%}, {int(global_distinct[source])} distinct values) — "
            f"it should be a .var column, not a layer"
        )
