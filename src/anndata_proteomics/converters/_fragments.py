"""Explode packed parallel-list fragment columns into one row per fragment.

DIA-NN-style reports pack per-fragment values as delimiter-joined lists inside each
precursor row (``Fragment.Info`` plus parallel ``Fragment.Quant.*`` lists, aligned by
index, often terminated by a trailing delimiter). A fragment-level AnnData needs one row
per fragment, so ``explode_fragments`` splits these columns and uses
``pandas.DataFrame.explode`` to fan the precursor row out — replicating the scalar
precursor columns and raising loudly if the parallel lists have mismatched lengths.
"""

from __future__ import annotations

import pandas as pd

from anndata_proteomics.rules.schema import Fragments


def _split_packed(value: object, delimiter: str) -> list[str]:
    """Split one packed cell into tokens, dropping a trailing empty terminator."""
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    text = text.rstrip(delimiter)  # drop DIA-NN's trailing-delimiter terminator
    if not text:
        return []
    return [token.strip() for token in text.split(delimiter)]


def explode_fragments(df: pd.DataFrame, fragments: Fragments) -> pd.DataFrame:
    """Return a long DataFrame with one row per fragment.

    Splits every ``fragments.value_columns`` entry on ``fragments.delimiter`` and explodes them
    together (index-aligned; mismatched list lengths raise ``ValueError``). The per-fragment
    ``fragments.label_output`` is either:

    - **column-labelled** (``label_strategy="column"``, e.g. ``Fragment.Info``): the token before
      ``/`` of ``fragments.label_column`` (``b4-unknown^1`` from ``b4-unknown^1/327.166``); or
    - **positional** (``label_strategy="positional"``, e.g. older DIA-NN with no ``Fragment.Info``):
      ``frag_0``, ``frag_1``, … by index within the precursor.

    Value columns are coerced to numeric; all other columns replicate per fragment.
    """
    value_columns = list(fragments.value_columns)
    if fragments.label_strategy == "column":
        packed_columns = [fragments.label_column, *value_columns]
    else:
        packed_columns = list(value_columns)
    missing = [c for c in packed_columns if c not in df.columns]
    if missing:
        raise KeyError(
            f"[fragments] references column(s) missing from the input: {missing}; "
            f"available: {list(df.columns)[:10]}…"
        )

    work = df.copy()
    for column in packed_columns:
        work[column] = work[column].map(lambda v: _split_packed(v, fragments.delimiter))

    if fragments.label_strategy == "positional":
        # Positional: a parallel index list per precursor, exploded alongside the values.
        work["_frag_pos"] = work[value_columns[0]].map(lambda tokens: list(range(len(tokens))))
        work = work.explode([*value_columns, "_frag_pos"], ignore_index=True)
        work = work.dropna(subset=[value_columns[0]]).reset_index(drop=True)
        work[fragments.label_output] = [f"frag_{int(pos)}" for pos in work["_frag_pos"]]
        work = work.drop(columns=["_frag_pos"])
    else:
        # Multi-column explode keeps the lists aligned and raises if their lengths differ.
        work = work.explode(packed_columns, ignore_index=True)
        # Precursors with no fragments explode to a NaN row; drop them.
        work = work.dropna(subset=[fragments.label_column]).reset_index(drop=True)
        work[fragments.label_output] = (
            work[fragments.label_column].astype(str).str.split("/").str[0]
        )
        # Drop the packed label column so the (now ~12x longer) frame doesn't carry a redundant
        # long string column.
        work = work.drop(columns=[fragments.label_column])

    # Coerce values to numeric so they ride the expanded frame as float64, not object strings.
    for column in value_columns:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    return work
