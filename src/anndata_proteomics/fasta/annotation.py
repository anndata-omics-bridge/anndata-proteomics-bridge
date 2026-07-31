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


@dataclass(frozen=True, slots=True)
class ResolvedCleavage:
    """A cleavage rule paired with the effective enzyme name."""

    rule: CleavageRule
    enzyme: str


@dataclass(frozen=True, slots=True)
class GeneName:
    """A UniProt header declares one concrete gene name."""

    value: str


@dataclass(frozen=True, slots=True)
class MissingGeneName:
    """A FASTA header contains no UniProt ``GN=`` field."""


type GeneNameResult = GeneName | MissingGeneName


@dataclass(frozen=True, slots=True)
class FastaHeaderDescription:
    """Text following the identifier in one FASTA header."""

    value: str


@dataclass(frozen=True, slots=True)
class MissingFastaHeaderDescription:
    """A FASTA header contains only its identifier."""


type FastaHeaderDescriptionResult = FastaHeaderDescription | MissingFastaHeaderDescription


@dataclass(frozen=True, slots=True)
class ParsedFastaHeader:
    """Explicit identifier and description result parsed from one FASTA header."""

    identifier: str
    description: FastaHeaderDescriptionResult


DEFAULT_CLEAVAGE = ResolvedCleavage(
    rule=_CLEAVAGE_RULES[_DEFAULT_ENZYME],
    enzyme=_DEFAULT_ENZYME,
)


@dataclass(frozen=True, slots=True)
class FastaAnnotationConfig:
    """Configuration for building a protein annotation table from FASTA."""

    identifiers: FastaConfig = field(default_factory=FastaConfig)
    is_uniprot: bool = True
    cleavage: ResolvedCleavage = DEFAULT_CLEAVAGE
    min_length: int = 7
    max_length: int = 30
    include_sequence: bool = False


DEFAULT_FASTA_ANNOTATION_CONFIG = FastaAnnotationConfig()


def resolve_cleavage_name(enzyme: str) -> ResolvedCleavage:
    """Resolve one enzyme name, falling back visibly to trypsin when unknown."""
    rule = _CLEAVAGE_RULES.get(enzyme)
    if rule is None:
        logger.warning(
            f"unknown enzyme {enzyme!r}; using {_DEFAULT_ENZYME} cleavage rule for peptide count"
        )
        return DEFAULT_CLEAVAGE
    return ResolvedCleavage(rule=rule, enzyme=enzyme)


def custom_cleavage(rule: CleavageRule) -> ResolvedCleavage:
    """Name an explicitly supplied cleavage rule for provenance."""
    return ResolvedCleavage(rule=rule, enzyme="custom")


def extract_gene_name(header: str) -> GeneNameResult:
    """Return the declared UniProt gene name or an explicit missing result."""
    match = _GN_RE.search(header)
    if match is None:
        return MissingGeneName()
    return GeneName(match.group(1))


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
    cleavage: CleavageRule = DEFAULT_CLEAVAGE.rule,
    min_length: int = 6,
    max_length: int = 30,
) -> int:
    """Count theoretical fully-cleaved peptides with ``min_length <= L < max_length``.

    The in-silico digest count behind the ``nr_peptides`` column. Mirrors the
    algorithm of prolfquapp's ``nr_tryptic_peptides`` (the upper bound is strict
    ``<``, not inclusive, even though the R docstring says "maximum length"), but
    the cleavage rule is configurable via *cleavage*, so it is not
    trypsin-specific.
    """
    cleavage_sites = _find_cleavage_sites(sequence, cleavage)
    starts = [0, *cleavage_sites]
    ends = [*cleavage_sites, len(sequence)]
    return sum(
        1
        for start, end in zip(starts, ends, strict=True)
        if min_length <= (end - start) < max_length
    )


def parse_header_id(header: str) -> ParsedFastaHeader:
    """Split a header into an identifier and explicit description result."""
    parts = header.split(maxsplit=1)
    fasta_id = parts[0].lstrip(">").rstrip(";")
    description: FastaHeaderDescriptionResult = (
        FastaHeaderDescription(parts[1]) if len(parts) > 1 else MissingFastaHeaderDescription()
    )
    return ParsedFastaHeader(identifier=fasta_id, description=description)


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
    records: list[tuple[str, FastaHeaderDescriptionResult, str]] = []
    for source in _iter_sources(sources):
        for record in iter_fasta(source):
            parsed = parse_header_id(record.header)
            records.append((parsed.identifier, parsed.description, record.sequence))

    resolved = resolve_fasta_config(
        (fasta_id for fasta_id, _, _ in records),
        config.identifiers,
    )

    if not records:
        return _empty_frame(include_sequence=config.include_sequence), resolved

    frame = pd.DataFrame(
        [
            {
                "fasta.id": fasta_id,
                "fasta.header": (
                    header.value if isinstance(header, FastaHeaderDescription) else pd.NA
                ),
                "sequence": sequence,
            }
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

    gene_names = frame["fasta.header"].map(
        lambda header: extract_gene_name(str(header)) if not pd.isna(header) else MissingGeneName()
    )
    present_gene_names = gene_names.map(lambda result: isinstance(result, GeneName))
    if int(present_gene_names.sum()) > 1:
        frame["gene_name"] = pd.Series(
            [result.value if isinstance(result, GeneName) else pd.NA for result in gene_names],
            index=frame.index,
            dtype="string",
        )

    frame["protein_length"] = frame["sequence"].map(len)
    frame["nr_peptides"] = frame["sequence"].map(
        lambda seq: count_peptides(
            str(seq),
            cleavage=config.cleavage.rule,
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
