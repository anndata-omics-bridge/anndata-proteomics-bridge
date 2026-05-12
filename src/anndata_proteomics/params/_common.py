"""Shared parser helpers."""

from __future__ import annotations

from collections.abc import Mapping, Sequence


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
