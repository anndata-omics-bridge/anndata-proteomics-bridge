"""Shared parser helpers."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import IO, Union

Source = Union[str, Path, IO[bytes], IO[str]]


def read_text(source: Source, *, errors: str = "strict") -> str:
    """Read a path, text file-like, or bytes file-like into text.

    Rewinds seekable streams first (a no-op on fresh sources) and decodes bytes
    as UTF-8. Centralizes the source-acquisition logic that each vendor parser
    used to re-implement; only the per-vendor parse step legitimately varies.
    """
    if hasattr(source, "read"):
        try:
            source.seek(0)
        except Exception:
            pass
        raw = source.read()
        if isinstance(raw, bytes):
            return raw.decode("utf-8", errors=errors)
        return raw
    return Path(source).read_text(encoding="utf-8", errors=errors)


def read_lines(source: Source, *, strip: bool = False) -> list[str]:
    """Read *source* into a list of lines, optionally stripping each line."""
    lines = read_text(source).splitlines()
    return [line.strip() for line in lines] if strip else lines


def format_tolerance_range(tolerance: Mapping[str, Sequence[float | int]]) -> str:
    """Format a mass-tolerance dict ``{unit: [low, high]}`` as a bracketed string.

    Output format matches ProteoBench's `_format_tolerance_range` so that
    existing expected-output CSVs round-trip unchanged.

    Examples
    --------
    >>> format_tolerance_range({"ppm": [-20.0, 20.0]})
    '[-20.0 ppm, 20.0 ppm]'
    """
    if not isinstance(tolerance, dict):
        raise ValueError(f"tolerance must be dict, got {type(tolerance).__name__}")

    unit_lookup = {"ppm": "ppm", "da": "Da"}
    for key, values in tolerance.items():
        normalized = str(key).lower()
        if normalized in unit_lookup:
            unit = unit_lookup[normalized]
            return "[" + ", ".join(f"{v} {unit}" for v in values) + "]"
    raise KeyError(f"unsupported tolerance unit(s): {list(tolerance.keys())}")
