"""Dispatch by file extension to the right tabular reader."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from anndata_proteomics.readers.tabular import read_csv, read_parquet, read_tsv


class UnknownFormat(ValueError):
    """Raised when a file extension has no registered reader."""


# .txt is treated as tab-delimited because MaxQuant evidence.txt and similar
# proteomics exports are tab-delimited despite the generic extension.
EXTENSION_TO_READER = {
    ".csv": read_csv,
    ".tsv": read_tsv,
    ".txt": read_tsv,
    ".parquet": read_parquet,
}


def read_table(path: Path | str) -> pd.DataFrame:
    """Read a tabular file, dispatching by extension.

    Raises UnknownFormat if the extension is not registered.
    """
    p = Path(path)
    reader = EXTENSION_TO_READER.get(p.suffix.lower())
    if reader is None:
        raise UnknownFormat(
            f"unsupported extension {p.suffix!r} for {p}; "
            f"known: {sorted(EXTENSION_TO_READER)}"
        )
    return reader(p)
