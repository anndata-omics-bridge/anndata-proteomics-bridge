"""Build a protein-annotation DataFrame from one or more FASTA files.

Replicates the data-extraction half of prolfquapp's
``get_annot_from_fasta()`` (R6_ProteinAnnotation.R / get_annot_from_FASTA.R):

  fasta.id, fasta.header, proteinname, gene_name (optional),
  protein_length, nr_peptides, sequence (optional)

Decoy records matching ``decoy_pattern`` are filtered out before
downstream columns are computed; the ``gene_name`` column is added only
when more than one record produces a UniProt-style match — matching
prolfquapp's gating rule so non-UniProt fastas don't end up with an
all-empty column.
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from loguru import logger

from anndata_proteomics.fasta.parser import FastaSource, iter_fasta

_GN_RE = re.compile(r" GN=(\S+) PE=")
_UNIPROT_MIDDLE_RE = re.compile(r".+\|(.+)\|.*")
_DEFAULT_DECOY_PATTERN = r"^REV_|^rev_"


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


def _parse_header_id(header: str) -> tuple[str, str]:
    """Split a header on the first whitespace into (fasta.id, fasta.header)."""
    parts = header.split(maxsplit=1)
    fasta_id = parts[0].lstrip(">").rstrip(";")
    fasta_header = parts[1] if len(parts) > 1 else ""
    return fasta_id, fasta_header


def _uniprot_proteinname(fasta_id: str) -> str:
    match = _UNIPROT_MIDDLE_RE.match(fasta_id)
    return match.group(1) if match else fasta_id


def fasta_to_dataframe(
    sources: FastaSource | Iterable[FastaSource],
    *,
    decoy_pattern: str = _DEFAULT_DECOY_PATTERN,
    is_uniprot: bool = True,
    cleavage: str | CleavageRule | None = None,
    min_length: int = 7,
    max_length: int = 30,
    include_sequence: bool = False,
) -> pd.DataFrame:
    """Read one or more FASTA inputs into a protein-annotation DataFrame.

    *cleavage* selects the protease rule for ``nr_peptides`` (an enzyme
    name, a :class:`CleavageRule`, or ``None`` for trypsin).
    """
    records: dict[str, tuple[str, str]] = {}
    for source in _iter_sources(sources):
        for record in iter_fasta(source):
            fasta_id, fasta_header = _parse_header_id(record.header)
            if fasta_id not in records:
                records[fasta_id] = (fasta_header, record.sequence)

    if not records:
        return _empty_frame(include_sequence=include_sequence)

    frame = pd.DataFrame(
        [
            {"fasta.id": fid, "fasta.header": hdr, "sequence": seq}
            for fid, (hdr, seq) in records.items()
        ]
    )
    frame = _drop_decoys(frame, decoy_pattern)
    frame = _add_annotation_columns(
        frame,
        is_uniprot=is_uniprot,
        cleavage=cleavage,
        min_length=min_length,
        max_length=max_length,
    )

    if not include_sequence:
        frame = frame.drop(columns=["sequence"])

    return frame


def _drop_decoys(frame: pd.DataFrame, decoy_pattern: str) -> pd.DataFrame:
    """Drop decoy rows matching *decoy_pattern*, warning on a suspiciously low hit rate."""
    if not decoy_pattern:
        return frame
    decoy_match = frame["fasta.id"].str.contains(decoy_pattern, regex=True)
    if 0 < decoy_match.mean() < 0.1:
        logger.warning(
            f"decoy pattern {decoy_pattern!r} matched only "
            f"{100 * decoy_match.mean():.1f}% of records"
        )
    return frame.loc[~decoy_match].reset_index(drop=True)


def _add_annotation_columns(
    frame: pd.DataFrame,
    *,
    is_uniprot: bool,
    cleavage: str | CleavageRule | None,
    min_length: int,
    max_length: int,
) -> pd.DataFrame:
    """Add proteinname, optional gene_name, protein_length, and peptide counts."""
    if is_uniprot:
        frame["proteinname"] = frame["fasta.id"].map(_uniprot_proteinname)
    else:
        frame["proteinname"] = frame["fasta.id"]

    gene_names = frame["fasta.header"].map(extract_gene_name)
    if (gene_names != "").sum() > 1:
        frame["gene_name"] = gene_names

    rule, _ = resolve_cleavage(cleavage)
    frame["protein_length"] = frame["sequence"].map(len)
    frame["nr_peptides"] = frame["sequence"].map(
        lambda seq: count_peptides(seq, cleavage=rule, min_length=min_length, max_length=max_length)
    )
    return frame


def _iter_sources(
    sources: FastaSource | Iterable[FastaSource],
) -> Iterable[FastaSource]:
    if isinstance(sources, str | Path):
        yield sources
        return
    if hasattr(sources, "read"):
        yield sources  # type: ignore[misc]
        return
    yield from sources


def _empty_frame(*, include_sequence: bool) -> pd.DataFrame:
    columns = [
        "fasta.id",
        "fasta.header",
        "proteinname",
        "protein_length",
        "nr_peptides",
    ]
    if include_sequence:
        columns.insert(2, "sequence")
    return pd.DataFrame(columns=columns)
