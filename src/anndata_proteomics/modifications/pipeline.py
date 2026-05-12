"""Apply a [modifications] rule to a DataFrame, adding normalized columns."""

from __future__ import annotations

import pandas as pd

from anndata_proteomics.modifications.apply_rules import (
    MapEntry,
    ModificationRule,
    apply_rule,
)
from anndata_proteomics.modifications.unimod_registry import resolve
from anndata_proteomics.rules.schema import Modifications


def _to_runtime_rule(mods: Modifications) -> ModificationRule:
    """Convert the validated TOML model into the runtime dataclass.

    Fills ``name``, ``target``, ``position``, ``mass_delta`` from the bundled
    Unimod registry; raises ``KeyError`` if any entry references an
    unknown accession.
    """
    runtime_entries: list[MapEntry] = []
    for e in mods.map:
        record = resolve(e.accession)
        runtime_entries.append(
            MapEntry(
                token=e.token,
                name=record.name,
                accession=record.accession,
                target=record.target,
                position=record.position,
                mass_delta=record.mass_delta,
            )
        )
    return ModificationRule(
        source_column=mods.source_column,
        token_pattern=mods.token_pattern or "",
        token_position=mods.token_position,
        case_sensitive=mods.case_sensitive,
        unknown_policy=mods.unknown_policy,
        sequence_column=mods.sequence_column,
        output_column=mods.output_column,
        entries=tuple(runtime_entries),
    )


def apply_modifications(df: pd.DataFrame, mods: Modifications) -> pd.DataFrame:
    """Add normalized modification columns to ``df`` based on ``mods``.

    Adds (and returns the same frame for convenience):
    - ``mods.output_column`` (default ``"proforma_sequence"``): ProForma string
    - ``"stripped_sequence"``: amino-acid-only sequence
    - ``"unknown_mod_tokens"``: list of unresolved vendor tokens per row

    The original ``mods.source_column`` is left untouched. If the parser
    mode is anything other than ``token_regex`` this is a no-op for now —
    ``already_proforma`` and ``separate_mod_column`` will be added when a
    rule that needs them lands.
    """
    if mods.parser != "token_regex":
        return df

    if mods.source_column not in df.columns:
        raise KeyError(
            f"[modifications].source_column={mods.source_column!r} not found "
            f"in DataFrame; available: {list(df.columns)[:10]}…"
        )

    runtime = _to_runtime_rule(mods)
    results = df[mods.source_column].astype(str).map(lambda s: apply_rule(s, runtime))
    df[mods.output_column] = [r.proforma_sequence for r in results]
    df["stripped_sequence"] = [r.stripped_sequence for r in results]
    df["unknown_mod_tokens"] = [r.unknown_tokens for r in results]
    return df
