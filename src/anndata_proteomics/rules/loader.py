"""Load and validate a parsing-rule TOML into a ParseRule."""

from __future__ import annotations

import tomllib
from pathlib import Path

from pydantic import ValidationError

from anndata_proteomics.rules.registry import find_rule
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
    software: str, quantification_level: str, file_version: str = "1"
) -> ParseRule:
    """Convenience: registry.find_rule(...) + load_rule(...)."""
    return load_rule(find_rule(software, quantification_level, file_version))
