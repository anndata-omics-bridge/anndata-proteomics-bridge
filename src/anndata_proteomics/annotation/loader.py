"""Load external observation-annotation tables without a schema model."""

from __future__ import annotations

import tomllib
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

ANNOTATION_SUFFIXES = frozenset({".csv", ".toml", ".tsv"})


@dataclass(slots=True)
class AnnotationTable:
    """A sample table plus the two fields needed to join it onto ``obs``."""

    samples: pd.DataFrame
    match_on: str = "index"
    key_field: str = "raw_file"
    source: Path | None = None


def load_annotation(path: Path | str) -> AnnotationTable:
    """Load a ProteoBench TOML or a delimited sample-annotation table.

    ProteoBench ``module_settings.toml`` files expose their table as top-level
    ``[[samples]]`` records. A sample may declare exact fallback identifiers as
    ``raw_file_alias`` or ``raw_file_aliases``. APB's earlier annotation-only
    TOMLs used ``[[obs.samples]]`` and may additionally set ``obs.match_on`` and
    ``obs.key_field``; both TOML shapes remain readable. CSV and TSV inputs use
    the conventional ``raw_file`` join column and match it against ``obs_names``.
    """
    source = Path(path).expanduser()
    if not source.exists():
        raise FileNotFoundError(source)
    if not source.is_file():
        raise ValueError(f"Annotation path is not a file: {source}")

    suffix = source.suffix.lower()
    if suffix == ".toml":
        annotation = _load_toml(source)
    elif suffix in {".csv", ".tsv"}:
        separator = "," if suffix == ".csv" else "\t"
        annotation = AnnotationTable(
            samples=pd.read_csv(source, sep=separator),
            source=source.resolve(),
        )
    else:
        supported = ", ".join(sorted(ANNOTATION_SUFFIXES))
        raise ValueError(
            f"Unsupported annotation format {suffix or '<none>'!r}; expected {supported}"
        )

    _validate_annotation_table(annotation)
    return annotation


def _load_toml(source: Path) -> AnnotationTable:
    document = tomllib.loads(source.read_text(encoding="utf-8"))
    obs = document.get("obs")
    if isinstance(obs, dict) and "samples" in obs:
        samples = obs["samples"]
        match_on = str(obs.get("match_on", "index"))
        key_field = str(obs.get("key_field", "raw_file"))
    else:
        samples = document.get("samples")
        match_on = "index"
        key_field = "raw_file"
    if not isinstance(samples, list) or not samples:
        raise ValueError(f"Annotation TOML has no [[samples]] or [[obs.samples]] records: {source}")
    if not all(isinstance(record, dict) for record in samples):
        raise ValueError(f"Annotation TOML samples must be tables: {source}")
    return AnnotationTable(
        samples=pd.DataFrame(samples),
        match_on=match_on,
        key_field=key_field,
        source=source.resolve(),
    )


def _validate_annotation_table(annotation: AnnotationTable) -> None:
    if annotation.samples.empty:
        raise ValueError("Annotation table must contain at least one sample row.")
    if annotation.key_field not in annotation.samples.columns:
        raise ValueError(
            f"Annotation table is missing key field {annotation.key_field!r}; "
            f"present columns: {list(annotation.samples.columns)}"
        )
