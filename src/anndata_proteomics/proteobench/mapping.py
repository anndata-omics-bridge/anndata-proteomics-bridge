"""ProteoBench accession-description normalization used before species matching."""

from __future__ import annotations

import csv
import hashlib
import re
from dataclasses import dataclass
from functools import cache
from pathlib import Path

import pandas as pd

from anndata_proteomics.modifications.unimod_registry import load_registry

_PROTEIN_SEPARATOR = re.compile(r"[;,]")
_UNIMOD_TAG = re.compile(r"\[(UNIMOD:\d+)\]", flags=re.IGNORECASE)
_FINAL_RESIDUE_MODS = re.compile(
    r"(?<=[A-Z])-?(?:\[UNIMOD:\d+\])+(?=(?:/\d+)?$)",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class ProteinMappingResult:
    """Mapped protein strings and compact mapper-use provenance."""

    proteins: pd.Series
    mapper_sha256: str
    mapper_entries: int
    matched_token_occurrences: int
    unmatched_token_occurrences: int


def map_reported_proteins(proteins: pd.Series) -> ProteinMappingResult:
    """Apply ProteoBench's bundled accession-to-description mapping."""
    mapper, mapper_sha256 = _protein_mapper()
    matched = 0
    unmatched = 0

    def normalize(value: str) -> str:
        nonlocal matched, unmatched
        tokens = [token.strip() for token in _PROTEIN_SEPARATOR.split(value)]
        normalized = []
        for token in tokens:
            if not token:
                continue
            replacement = mapper.get(token)
            if replacement is None:
                unmatched += 1
                normalized.append(token)
            else:
                matched += 1
                normalized.append(replacement)
        return ";".join(normalized)

    # Preserve pandas' explicit missing scalar. Missing protein assignments are not
    # empty accessions and must never enter token-mapping control flow.
    mapped = proteins.astype("string").map(normalize, na_action="ignore")
    return ProteinMappingResult(
        proteins=mapped,
        mapper_sha256=mapper_sha256,
        mapper_entries=len(mapper),
        matched_token_occurrences=matched,
        unmatched_token_occurrences=unmatched,
    )


def render_proteobench_features(
    features: pd.Series,
    *,
    drop_final_residue_modifications: bool = False,
) -> pd.Series:
    """Render canonical APB ProForma tags with ProteoBench's legacy names.

    APB deliberately retains Unimod accessions in its canonical feature axis,
    whereas ProteoBench's intermediate CSV uses modification names. This
    compatibility rendering is used only for the reconstructed intermediate;
    it never changes the AnnData feature identifiers.
    """
    names = {accession.upper(): entry.name for accession, entry in load_registry().items()}

    def replace(match: re.Match[str]) -> str:
        accession = match.group(1).upper()
        name = names.get(accession)
        return f"[{name}]" if name is not None else match.group(0)

    rendered = features.astype("string")
    if rendered.isna().any():
        raise ValueError("ProteoBench feature identifiers must not be missing")
    if drop_final_residue_modifications:
        # ProteoBench 0.17's ``before_aa = false`` parser does not visit the
        # position after the final residue. Preserve that behavior in the
        # compatibility table without changing APB's canonical feature axis.
        rendered = rendered.str.replace(_FINAL_RESIDUE_MODS, "", regex=True)
    return rendered.str.replace(_UNIMOD_TAG, replace, regex=True)


@cache
def _protein_mapper() -> tuple[dict[str, str], str]:
    mapper_path = Path(__file__).with_name("mapper.csv")
    source = mapper_path.read_bytes()
    with mapper_path.open(encoding="utf-8", newline="") as handle:
        mapper = {
            row["gene_name"]: row["description"]
            for row in csv.DictReader(handle)
            if row.get("gene_name") and row.get("description")
        }
    return mapper, hashlib.sha256(source).hexdigest()
