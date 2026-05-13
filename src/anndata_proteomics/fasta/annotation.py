"""Build a protein-annotation DataFrame from one or more FASTA files.

Replicates the data-extraction half of prolfquapp's
``get_annot_from_fasta()`` (R6_ProteinAnnotation.R / get_annot_from_FASTA.R):

  fasta.id, fasta.header, proteinname, gene_name (optional),
  protein_length, nr_tryptic_peptides, sequence (optional)

Decoy records matching ``decoy_pattern`` are filtered out before
downstream columns are computed; the ``gene_name`` column is added only
when more than one record produces a UniProt-style match — matching
prolfquapp's gating rule so non-UniProt fastas don't end up with an
all-empty column.
"""

from __future__ import annotations

import logging
import re
from collections.abc import Iterable
from pathlib import Path

import pandas as pd

from anndata_proteomics.fasta.parser import FastaSource, iter_fasta

logger = logging.getLogger(__name__)

_GN_RE = re.compile(r" GN=(\S+) PE=")
_UNIPROT_MIDDLE_RE = re.compile(r".+\|(.+)\|.*")
_CLEAVAGE_RE = re.compile(r"[KR](?!P|$)")
_DEFAULT_DECOY_PATTERN = r"^REV_|^rev_"


def extract_gene_name(header: str) -> str:
    """Return the UniProt ``GN=`` value from a FASTA header, or ``""`` if absent."""
    match = _GN_RE.search(header)
    return match.group(1) if match else ""


def _find_cleavage_sites(sequence: str) -> list[int]:
    """Return 1-based end positions of tryptic cleavage sites (K|R, not before P)."""
    return [m.end() for m in _CLEAVAGE_RE.finditer(sequence.upper())]


def count_tryptic_peptides(
    sequence: str, *, min_length: int = 6, max_length: int = 30
) -> int:
    """Count fully-tryptic peptides with ``min_length <= L < max_length``.

    Mirrors prolfquapp's ``nr_tryptic_peptides``: the upper bound is
    strict (``<``), not inclusive, even though the R docstring says
    "maximum length".
    """
    cleavage_sites = _find_cleavage_sites(sequence)
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
    min_length: int = 7,
    max_length: int = 30,
    include_sequence: bool = False,
) -> pd.DataFrame:
    """Read one or more FASTA inputs into a protein-annotation DataFrame."""
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

    if decoy_pattern:
        decoy_match = frame["fasta.id"].str.contains(decoy_pattern, regex=True)
        if 0 < decoy_match.mean() < 0.1:
            logger.warning(
                "decoy pattern %r matched only %.1f%% of records",
                decoy_pattern,
                100 * decoy_match.mean(),
            )
        frame = frame.loc[~decoy_match].reset_index(drop=True)

    if is_uniprot:
        frame["proteinname"] = frame["fasta.id"].map(_uniprot_proteinname)
    else:
        frame["proteinname"] = frame["fasta.id"]

    gene_names = frame["fasta.header"].map(extract_gene_name)
    if (gene_names != "").sum() > 1:
        frame["gene_name"] = gene_names

    frame["protein_length"] = frame["sequence"].map(len)
    frame["nr_tryptic_peptides"] = frame["sequence"].map(
        lambda seq: count_tryptic_peptides(
            seq, min_length=min_length, max_length=max_length
        )
    )

    if not include_sequence:
        frame = frame.drop(columns=["sequence"])

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
        "nr_tryptic_peptides",
    ]
    if include_sequence:
        columns.insert(2, "sequence")
    return pd.DataFrame(columns=columns)
