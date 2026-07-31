"""Build a protein-annotation DataFrame from one or more FASTA files.

Replicates the data-extraction half of prolfquapp's
``get_annot_from_fasta()`` (R6_ProteinAnnotation.R / get_annot_from_FASTA.R):

  fasta.id, fasta.header, proteinname, gene_name (optional),
  protein_length, nr_peptides, sequence (optional)

Decoy and contaminant records are retained and classified.  FASTA annotation
must never remove quantified features merely because a search engine reported
a decoy or contaminant identification.  The ``gene_name`` column is added only
when more than one record produces a UniProt-style match — matching
prolfquapp's gating rule so non-UniProt FASTAs don't end up with an all-empty
column.
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO, TypeGuard

import pandas as pd
from loguru import logger

from anndata_proteomics.fasta.config import (
    FastaConfig,
    ResolvedFastaConfig,
    matches_any,
    resolve_fasta_config,
)
from anndata_proteomics.fasta.parser import FastaSource, FastaSources, iter_fasta

_GN_RE = re.compile(r" GN=(\S+) PE=")
_UNIPROT_MIDDLE_RE = re.compile(r".+\|(.+)\|.*")


@dataclass(frozen=True, slots=True)
class CleavageRule:
    """A protease cleavage rule: a residue pattern plus which side it cuts.

    ``pattern`` matches the residue adjacent to a cut. ``after=True`` cuts
    C-terminal to the match (the common case: trypsin cuts after K/R);
    ``after=False`` cuts N-terminal (Asp-N cuts before D).
    """

    pattern: re.Pattern[str]
    after: bool = True


@dataclass(frozen=True, slots=True)
class FastaAnnotationConfig:
    """Configuration for building a protein annotation table from FASTA."""

    identifiers: FastaConfig = field(default_factory=FastaConfig)
    is_uniprot: bool = True
    cleavage: str | CleavageRule | None = None
    min_length: int = 7
    max_length: int = 30
    include_sequence: bool = False


DEFAULT_FASTA_ANNOTATION_CONFIG = FastaAnnotationConfig()


# Enzyme → cleavage rule, keyed by the canonical display names emitted by
# ``params.model.Parameters.enzyme`` (the ``_ENZYME_MAP`` values) so the two
# cannot drift. 99% of searches are trypsin, but Lys-C / Glu-C / etc. happen,
# so the peptide count uses the *actual* enzyme rather than assuming trypsin.
_CLEAVAGE_RULES: dict[str, CleavageRule] = {
    "Trypsin": CleavageRule(re.compile(r"[KR](?!P)")),
    "Trypsin/P": CleavageRule(re.compile(r"[KR]")),
    "Lys-C": CleavageRule(re.compile(r"K(?!P)")),
    "Arg-C": CleavageRule(re.compile(r"R(?!P)")),
    "Glu-C": CleavageRule(re.compile(r"[DE](?!P)")),
    "Chymotrypsin": CleavageRule(re.compile(r"[FYW](?!P)")),
    "Asp-N": CleavageRule(re.compile(r"D"), after=False),
}
_DEFAULT_ENZYME = "Trypsin"


def resolve_cleavage(cleavage: str | CleavageRule | None) -> tuple[CleavageRule, str]:
    """Resolve a cleavage spec to ``(rule, effective_enzyme_name)``.

    ``None`` is the documented trypsin default (no warning). An unknown enzyme
    name warns once and falls back to trypsin. A pre-built :class:`CleavageRule`
    is returned verbatim with the name ``"custom"``.
    """
    if isinstance(cleavage, CleavageRule):
        return cleavage, "custom"
    if cleavage is None:
        return _CLEAVAGE_RULES[_DEFAULT_ENZYME], _DEFAULT_ENZYME
    rule = _CLEAVAGE_RULES.get(cleavage)
    if rule is None:
        logger.warning(
            f"unknown enzyme {cleavage!r}; using {_DEFAULT_ENZYME} cleavage rule for peptide count"
        )
        return _CLEAVAGE_RULES[_DEFAULT_ENZYME], _DEFAULT_ENZYME
    return rule, cleavage


def extract_gene_name(header: str) -> str:
    """Return the UniProt ``GN=`` value from a FASTA header, or ``""`` if absent."""
    match = _GN_RE.search(header)
    return match.group(1) if match else ""


def _find_cleavage_sites(sequence: str, rule: CleavageRule) -> list[int]:
    """Return the cut positions for *rule* in *sequence* (0-based offsets).

    For an ``after`` rule the cut is C-terminal to the matched residue
    (``m.end()``); a cut at the very C-terminus is dropped (zero-length tail).
    For a ``before`` rule (Asp-N) the cut is N-terminal (``m.start()``); a cut
    at position 0 is dropped (zero-length head).
    """
    seq = sequence.upper()
    if rule.after:
        sites = [m.end() for m in rule.pattern.finditer(seq)]
        return [s for s in sites if s != len(seq)]
    sites = [m.start() for m in rule.pattern.finditer(seq)]
    return [s for s in sites if s != 0]


def count_peptides(
    sequence: str,
    *,
    cleavage: str | CleavageRule | None = None,
    min_length: int = 6,
    max_length: int = 30,
) -> int:
    """Count theoretical fully-cleaved peptides with ``min_length <= L < max_length``.

    The in-silico digest count behind the ``nr_peptides`` column. Mirrors the
    algorithm of prolfquapp's ``nr_tryptic_peptides`` (the upper bound is strict
    ``<``, not inclusive, even though the R docstring says "maximum length"), but
    the cleavage rule is configurable via *cleavage* — an enzyme name, a
    :class:`CleavageRule`, or ``None`` for trypsin — so it is not trypsin-specific.
    """
    rule, _ = resolve_cleavage(cleavage)
    cleavage_sites = _find_cleavage_sites(sequence, rule)
    starts = [0, *cleavage_sites]
    ends = [*cleavage_sites, len(sequence)]
    return sum(
        1
        for start, end in zip(starts, ends, strict=True)
        if min_length <= (end - start) < max_length
    )


def parse_header_id(header: str) -> tuple[str, str]:
    """Split a header on the first whitespace into (fasta.id, fasta.header)."""
    parts = header.split(maxsplit=1)
    fasta_id = parts[0].lstrip(">").rstrip(";")
    fasta_header = parts[1] if len(parts) > 1 else ""
    return fasta_id, fasta_header


def uniprot_proteinname(fasta_id: str) -> str:
    """Return the accession while preserving a prefix before ``sp|``/``tr|``.

    A decoy such as ``REV_sp|P12345|NAME`` must remain ``REV_P12345``;
    collapsing it to ``P12345`` would make it indistinguishable from its target.
    """
    match = re.match(r"^(.*?)(?:sp|tr)\|([^|]+)\|", fasta_id)
    if match:
        return f"{match.group(1)}{match.group(2)}"
    match = _UNIPROT_MIDDLE_RE.match(fasta_id)
    return match.group(1) if match else fasta_id


def fasta_to_dataframe(
    sources: FastaSources,
    config: FastaAnnotationConfig = DEFAULT_FASTA_ANNOTATION_CONFIG,
) -> pd.DataFrame:
    """Read one or more FASTA inputs into a protein-annotation DataFrame.

    Configuration controls identifier classification, accession parsing,
    digestion, and whether the raw sequence is retained.
    """
    frame, _ = fasta_to_dataframe_with_config(sources, config)
    return frame


def fasta_to_dataframe_with_config(
    sources: FastaSources,
    config: FastaAnnotationConfig = DEFAULT_FASTA_ANNOTATION_CONFIG,
) -> tuple[pd.DataFrame, ResolvedFastaConfig]:
    """Build the annotation frame and return its resolved ID configuration."""
    records: list[tuple[str, str, str]] = []
    for source in _iter_sources(sources):
        for record in iter_fasta(source):
            fasta_id, fasta_header = parse_header_id(record.header)
            records.append((fasta_id, fasta_header, record.sequence))

    resolved = resolve_fasta_config(
        (fasta_id for fasta_id, _, _ in records),
        config.identifiers,
    )

    if not records:
        return _empty_frame(include_sequence=config.include_sequence), resolved

    frame = pd.DataFrame(
        [
            {"fasta.id": fasta_id, "fasta.header": header, "sequence": sequence}
            for fasta_id, header, sequence in records
        ]
    )
    frame["is_decoy"] = frame["fasta.id"].map(
        lambda value: matches_any(str(value), resolved.decoy.patterns)
    )
    frame["is_contaminant"] = frame["fasta.id"].map(
        lambda value: matches_any(str(value), resolved.contaminant.patterns)
    )
    frame = _add_annotation_columns(
        frame,
        config,
    )

    if not config.include_sequence:
        frame = frame.drop(columns=["sequence"])

    return frame, resolved


def _add_annotation_columns(
    frame: pd.DataFrame,
    config: FastaAnnotationConfig,
) -> pd.DataFrame:
    """Add proteinname, optional gene_name, protein_length, and peptide counts."""
    if config.is_uniprot:
        frame["proteinname"] = frame["fasta.id"].map(uniprot_proteinname)
    else:
        frame["proteinname"] = frame["fasta.id"]

    gene_names = frame["fasta.header"].map(extract_gene_name)
    if (gene_names != "").sum() > 1:
        frame["gene_name"] = gene_names

    rule, _ = resolve_cleavage(config.cleavage)
    frame["protein_length"] = frame["sequence"].map(len)
    frame["nr_peptides"] = frame["sequence"].map(
        lambda seq: count_peptides(
            str(seq),
            cleavage=rule,
            min_length=config.min_length,
            max_length=config.max_length,
        )
    )
    return frame


def _is_text_stream(
    source: FastaSources,
) -> TypeGuard[IO[str]]:
    """Return whether *source* is one open text stream."""
    return hasattr(source, "read")


def _iter_sources(
    sources: FastaSources,
) -> Iterable[FastaSource]:
    if isinstance(sources, str | Path):
        yield sources
        return
    if _is_text_stream(sources):
        yield sources
        return
    yield from sources


def materialize_sources(
    sources: FastaSources,
) -> list[FastaSource]:
    """Materialize source iterables once so scanning cannot erase provenance."""
    return list(_iter_sources(sources))


def describe_sources(sources: FastaSources) -> list[str]:
    """Human-readable list of FASTA sources, for provenance.

    Paths are recorded as strings; inline FASTA text (the parser's string-content
    path) as ``<inline-fasta>``; open streams as ``<stream>``. Shared by both FASTA
    annotators so the same input is described identically everywhere.
    """
    items = materialize_sources(sources)
    out: list[str] = []
    for item in items:
        if isinstance(item, Path):
            out.append(str(item))
        elif isinstance(item, str):
            out.append("<inline-fasta>" if is_inline_fasta(item) else item)
        else:
            out.append("<stream>")
    return out


def is_inline_fasta(text: str) -> bool:
    """True when *text* is FASTA content rather than a path (matches the parser)."""
    return "\n" in text or text.lstrip().startswith(">")


def _empty_frame(*, include_sequence: bool) -> pd.DataFrame:
    columns = [
        "fasta.id",
        "fasta.header",
        "is_decoy",
        "is_contaminant",
        "proteinname",
        "protein_length",
        "nr_peptides",
    ]
    if include_sequence:
        columns.insert(2, "sequence")
    return pd.DataFrame(columns=columns)
