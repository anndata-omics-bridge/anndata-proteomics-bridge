"""Linux-filename-style sanitisation for AnnData column names.

Implements the rule documented in
``anndata_omics_bridge/docs/conventions.md`` (applied to ``obs.columns`` /
``var.columns`` / layer names — never to row IDs or ``uns`` keys).
"""

from __future__ import annotations

import re
import unicodedata


def sanitize_name(name: str) -> str:
    """Sanitise a single column/layer name.

    - Case is preserved (Linux filesystems are case-sensitive).
    - Dots are kept (DIA-NN's ``Protein.Group`` stays ``Protein.Group``).
    - Whitespace and other special characters are replaced with ``_``.
    - Repeated underscores are collapsed; leading/trailing ``_`` and ``.`` are stripped.
    - Empty results fall back to ``col``.
    """
    name = unicodedata.normalize("NFKD", name).encode("ascii", "ignore").decode("ascii")
    name = re.sub(r"[^A-Za-z0-9_.]", "_", name)
    name = re.sub(r"_+", "_", name).strip("_.")
    return name or "col"


def sanitize_columns(names: list[str]) -> list[str]:
    """Sanitise a list of names, raising on a post-sanitisation collision.

    Per the conventions "Conflict policy": if two distinct originals sanitise to
    the same string, raise rather than silently de-duplicate. The right place to
    resolve a collision is the upstream name, not a generated suffix.
    """
    sanitised = [sanitize_name(n) for n in names]
    groups: dict[str, list[str]] = {}
    for original, clean in zip(names, sanitised):
        groups.setdefault(clean, []).append(original)
    collisions = {clean: origs for clean, origs in groups.items() if len(set(origs)) > 1}
    if collisions:
        lines = [
            f"  {original!r} -> {clean!r}"
            for clean, origs in collisions.items()
            for original in origs
        ]
        raise ValueError("column-name collision after sanitisation:\n" + "\n".join(lines))
    return sanitised
