"""Param-driven conversion orchestration: vendor file + params → AnnData / MuData.

The non-UI core that backs the ``apb convert`` CLI. Given a vendor table plus the software
version (parsed from the parameter file), it selects the matching parsing-rule variant, converts
each quantification level, and (for the default target) wraps them into a MuData on a shared run
axis. No GUI / marimo / test-data-cache dependency — this is plain library code.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Iterable, Mapping
from pathlib import Path

import pandas as pd
from anndata import AnnData
from mudata import MuData

from anndata_proteomics.converters.recognize import matches
from anndata_proteomics.params.anndata_io import write_search_parameters
from anndata_proteomics.params.registry import parse_params
from anndata_proteomics.readers.summary import store_quantification_summary
from anndata_proteomics.rules.loader import load_rule, resolve_rule_for_version
from anndata_proteomics.rules.registry import iter_packaged_rules
from anndata_proteomics.rules.schema import ParseRule, QuantificationLevel

# Quantification levels, coarse to fine. Not every vendor exposes every level.
LEVELS: tuple[QuantificationLevel, ...] = (
    "ion",
    "peptidoform",
    "peptide",
    "protein",
    "fragment",
)
MUDATA = "mudata"
# Per-level var_names prefix so modalities don't collide on the global axis.
_PREFIX = {
    "fragment": "frg:",
    "ion": "ion:",
    "peptidoform": "pfm:",
    "peptide": "pep:",
    "protein": "prt:",
}


def software_slug(software_name: str) -> str:
    """Map a catalog ``software_name`` (e.g. "DIA-NN") to a parsing-rule vendor slug ("diann")."""
    return re.sub(r"[^a-z0-9]", "", software_name.lower())


def recognize_software(headers: Iterable[str]) -> str | None:
    """Infer the vendor slug from a file's column headers.

    Unlike ``recognize`` (which needs a *unique* rule match), a single vendor file can match
    several of its own level rules (e.g. a DIA-NN report matches ion/fragment/protein). This
    returns the vendor slug when every matching packaged rule belongs to the same vendor, else
    ``None`` (zero matches, or an ambiguous match spanning multiple vendors).
    """
    header_set = set(headers)
    slugs = {
        software_slug(rule.software_name)
        for path in iter_packaged_rules()
        if matches(header_set, rule := load_rule(path))
    }
    return next(iter(slugs)) if len(slugs) == 1 else None


def param_version(param_path: Path | None, slug: str) -> str | None:
    """Software version parsed from the param file, or None if absent/unparseable."""
    if param_path is None:
        return None
    try:
        return parse_params(param_path, software=slug).software_version
    except Exception:  # noqa: BLE001 — any parse failure means "version unknown"
        return None


def select_rule(
    slug: str,
    level: QuantificationLevel,
    version: str | None,
    headers: Iterable[str],
) -> ParseRule:
    """Resolve the rule for (slug, level) at ``version`` and validate it against ``headers``.

    Raises ValueError if no rule variant covers the version, or if the file's columns don't match
    the version-selected rule (no silent fallback — the caller fixes the version / param file).
    """
    rule = resolve_rule_for_version(slug, level, version)
    if rule is None:
        raise ValueError(f"{slug} {level}: no rule covers software version {version!r}")
    header_set = set(headers)
    label_col = (
        rule.fragments.label_column
        if rule.fragments is not None and rule.fragments.label_strategy == "column"
        else None
    )
    if not matches(header_set, rule) or (label_col is not None and label_col not in header_set):
        raise ValueError(
            f"{slug} {level}: file columns don't match the rule for software version "
            f"{version!r} — verify the version / provide the right param file"
        )
    return rule


def convertible_levels(
    slug: str,
    version: str | None,
    headers: Iterable[str],
) -> list[QuantificationLevel]:
    """Levels whose version-selected rule both exists and matches this file's columns."""
    header_set = set(headers)
    out = []
    for level in LEVELS:
        try:
            select_rule(slug, level, version, header_set)
        except (LookupError, ValueError):
            continue
        out.append(level)
    return out


def available_targets(slug: str, version: str | None, headers: Iterable[str]) -> list[str]:
    """Convertible levels plus the all-level MuData target when any level resolves."""
    levels = convertible_levels(slug, version, headers)
    targets = [str(level) for level in levels]
    if levels:
        targets.append(MUDATA)
    return targets


def matching_rules(
    rules: Mapping[QuantificationLevel, ParseRule],
    headers: Iterable[str],
) -> dict[QuantificationLevel, ParseRule]:
    """Return document levels whose required source columns match a vendor table."""
    header_set = set(headers)
    matched = {}
    for level, rule in rules.items():
        label_column = (
            rule.fragments.label_column
            if rule.fragments is not None and rule.fragments.label_strategy == "column"
            else None
        )
        if matches(header_set, rule) and (label_column is None or label_column in header_set):
            matched[level] = rule
    return matched


def _noop(_msg: str) -> None:
    pass


def convert_level(
    df: pd.DataFrame,
    slug: str,
    level: QuantificationLevel,
    version: str | None,
    *,
    params_path: Path | str | None = None,
    log: Callable[[str], None] = _noop,
) -> AnnData:
    from anndata_proteomics.converters.assemble import convert

    rule = select_rule(slug, level, version, df.columns)
    adata = convert(df, rule, params_path=params_path)
    log(f"  {level}: {adata.shape[0]} obs × {adata.shape[1]} var")
    return adata


def build_mudata(
    df: pd.DataFrame,
    slug: str,
    version: str | None,
    *,
    params_path: Path | str | None = None,
    log: Callable[[str], None] = _noop,
) -> MuData:
    """Build a MuData over the levels whose version-selected rule fits this file (shared run axis).

    Levels the version doesn't provide (e.g. fragment on DIA-NN 2.x) are skipped, not failed.
    """
    resolvable = set(convertible_levels(slug, version, df.columns))
    skipped = [level for level in LEVELS if level not in resolvable]
    if skipped:
        log(f"skipping levels not provided by software version {version!r}: {skipped}")
    if not resolvable:
        raise ValueError(f"{slug}: no level resolves for software version {version!r}")

    rules: dict[QuantificationLevel, ParseRule] = {
        level: select_rule(slug, level, version, df.columns)
        for level in LEVELS
        if level in resolvable
    }
    return build_mudata_from_rules(
        df,
        rules,
        params_path=params_path,
        software=slug,
        log=log,
    )


def build_mudata_from_rules(
    df: pd.DataFrame,
    rules: Mapping[QuantificationLevel, ParseRule],
    *,
    params_path: Path | str | None = None,
    software: str | None = None,
    log: Callable[[str], None] = _noop,
) -> MuData:
    """Build MuData from already selected effective rules."""
    import mudata

    if not rules:
        raise ValueError("no levels supplied")
    mods: dict[str, AnnData] = {}
    for level in LEVELS:
        if level not in rules:
            continue
        log(f"converting level: {level}")
        from anndata_proteomics.converters.assemble import convert

        adata = convert(df.copy(), rules[level], params_path=params_path)
        log(f"  {level}: {adata.shape[0]} obs × {adata.shape[1]} var")
        adata.var_names = [_PREFIX[level] + str(v) for v in adata.var_names]
        mods[level] = adata
    # Adopt the mudata 0.4 default now (no auto-pull of per-modality obs/var into the global
    # frames); each modality keeps its own obs/var. Silences the 0.3 FutureWarning.
    with mudata.set_options(pull_on_update=False):
        md = MuData(mods, axis=0)
    if params_path is not None and software is not None:
        params = parse_params(params_path, software=software)
        write_search_parameters(md, params, source_path=str(params_path))
    store_quantification_summary(md)
    log(f"  MuData: {md.n_obs} obs × {sum(a.n_vars for a in mods.values())} var, {len(mods)} mods")
    return md
