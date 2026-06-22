"""Load and validate a parsing-rule TOML into a ParseRule."""

from __future__ import annotations

import tomllib
from pathlib import Path

from pydantic import ValidationError

from anndata_proteomics.rules.registry import find_rule, resolve_rule_path
from anndata_proteomics.rules.schema import ParseRule


def load_rule(path: Path | str) -> ParseRule:
    """Load a TOML file and validate it as a ParseRule.

    Raises FileNotFoundError if the path doesn't exist.
    On pydantic validation failure, the file path is attached as an exception note
    so it shows up in the traceback alongside pydantic's field-level message.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    data = tomllib.loads(p.read_text())
    try:
        return ParseRule.model_validate(data)
    except ValidationError as e:
        e.add_note(f"in {p}")
        raise


def load_packaged_rule(
    software: str, quantification_level: str, version: str | None = None
) -> ParseRule:
    """Load the packaged rule for (software, level) at ``version`` (None → version-agnostic root)."""
    return load_rule(find_rule(software, quantification_level, version))


def resolve_rule_for_version(
    software: str, quantification_level: str, version: str | None
) -> ParseRule | None:
    """The rule whose version subfolder covers ``version``, or ``None`` if no variant applies.

    DIA-NN report columns vary by version; this picks the right variant by the software version
    (parsed from the param file), via ``registry.resolve_rule_path``. ``None`` means the level
    is not available for that version (e.g. fragment on DIA-NN 2.x).
    """
    path = resolve_rule_path(software, quantification_level, version)
    return load_rule(path) if path is not None else None
