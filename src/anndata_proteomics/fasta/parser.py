"""Minimal FASTA reader.

Yields ``FastaRecord(header, sequence)`` tuples from a path, an open
text stream, or a raw FASTA string. No biology semantics here — header
parsing and per-protein derivations live in :mod:`annotation`.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import IO, Union


@dataclass(frozen=True, slots=True)
class FastaRecord:
    header: str
    sequence: str


FastaSource = Union[str, Path, IO[str]]


def iter_fasta(source: FastaSource) -> Iterator[FastaRecord]:
    """Yield ``FastaRecord`` instances from a file path, stream, or raw FASTA text."""
    if isinstance(source, Path):
        with source.open("r", encoding="utf-8") as handle:
            yield from _iter_lines(handle)
        return
    if isinstance(source, str):
        if "\n" in source or source.lstrip().startswith(">"):
            yield from _iter_lines(StringIO(source))
        else:
            with Path(source).open("r", encoding="utf-8") as handle:
                yield from _iter_lines(handle)
        return
    yield from _iter_lines(source)


def _iter_lines(lines: Iterable[str]) -> Iterator[FastaRecord]:
    header: str | None = None
    seq_parts: list[str] = []
    for raw in lines:
        line = raw.rstrip("\r\n")
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield FastaRecord(header=header, sequence="".join(seq_parts))
            header = line[1:]
            seq_parts = []
        else:
            seq_parts.append(line.strip())
    if header is not None:
        yield FastaRecord(header=header, sequence="".join(seq_parts))
