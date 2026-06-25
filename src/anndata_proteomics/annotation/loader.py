"""Load and validate an annotation TOML into an AnnotationSpec."""

from __future__ import annotations

import tomllib
from pathlib import Path

from pydantic import ValidationError

from anndata_proteomics.annotation.schema import AnnotationSpec


def load_annotation(path: Path | str) -> AnnotationSpec:
    """Load a TOML file and validate it as an AnnotationSpec.

    Raises FileNotFoundError if the path doesn't exist. On pydantic validation failure the
    file path is attached as an exception note (same pattern as ``rules.loader.load_rule``).
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    data = tomllib.loads(p.read_text())
    try:
        return AnnotationSpec.model_validate(data)
    except ValidationError as e:
        e.add_note(f"in {p}")
        raise
