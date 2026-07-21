"""Load and validate an annotation JSON document into an AnnotationSpec."""

from __future__ import annotations

import json
from pathlib import Path

from pydantic import ValidationError

from anndata_proteomics.annotation.schema import AnnotationSpec


def parse_annotation(text: str, *, path: Path | str = "<memory>") -> AnnotationSpec:
    """Parse and validate an annotation JSON document from text."""
    source = Path(path)
    try:
        data = json.loads(text)
    except json.JSONDecodeError as exc:
        exc.add_note(f"in {source}")
        raise
    try:
        return AnnotationSpec.model_validate(data)
    except ValidationError as exc:
        exc.add_note(f"in {source}")
        raise


def load_annotation(path: Path | str) -> AnnotationSpec:
    """Load a JSON file and validate it as an AnnotationSpec.

    Raises FileNotFoundError if the path doesn't exist. On pydantic validation failure the
    file path is attached as an exception note (same pattern as ``rules.loader.load_rule``).
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    return parse_annotation(p.read_text(encoding="utf-8"), path=p)
