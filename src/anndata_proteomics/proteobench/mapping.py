"""ProteoBench accession-description normalization used before species matching."""

from __future__ import annotations

import csv
import hashlib
import re
from dataclasses import dataclass
from functools import cache
from pathlib import Path
from typing import Any

import pandas as pd

from anndata_proteomics.modifications.unimod_registry import load_registry
from anndata_proteomics.params.anndata_io import read_search_parameters
from anndata_proteomics.proteobench.config import ModificationParserSettings

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

    mapped = proteins.map(normalize)
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

    rendered = features.astype("string").fillna("")
    if drop_final_residue_modifications:
        # ProteoBench 0.17's ``before_aa = false`` parser does not visit the
        # position after the final residue. Preserve that behavior in the
        # compatibility table without changing APB's canonical feature axis.
        rendered = rendered.str.replace(_FINAL_RESIDUE_MODS, "", regex=True)
    return rendered.str.replace(_UNIMOD_TAG, replace, regex=True)


def parse_proteobench_features(
    raw_features: pd.Series,
    canonical_features: pd.Series,
    settings: ModificationParserSettings,
    *,
    level: str,
    target: Any,
) -> pd.Series:
    """Reproduce ProteoBench's per-tool modification rendering.

    The source notation retained by APB is used because ProteoBench's historical
    output preserves tool-specific parsing details, including modification-name
    case and its ``before_aa = false`` terminal-residue behavior. Canonical APB
    ProForma identifiers remain unchanged.
    """
    pattern = re.compile(settings.pattern)

    def convert(value: Any) -> str:
        source = "" if pd.isna(value) else str(value)
        lowered = pattern.sub(lambda match: match.group(0).lower(), source)
        matches = list(pattern.finditer(lowered))
        positions = [
            _count_sequence_characters(
                lowered[: match.start()],
                isalpha=settings.isalpha,
                isupper=settings.isupper,
            )
            for match in matches
        ]
        modifications = [
            settings.modification_dict.get(match.group(0), match.group(0)) for match in matches
        ]
        by_position = dict(zip(positions, modifications, strict=True))
        stripped = _strip_sequence(
            lowered,
            isalpha=settings.isalpha,
            isupper=settings.isupper,
        )
        rendered = ""
        for index, amino_acid in enumerate(stripped):
            if settings.before_aa:
                rendered += amino_acid
            modification = by_position.get(index)
            if modification is not None:
                if index == 0:
                    rendered += f"[{modification}]-"
                elif index == len(stripped):
                    rendered += f"-[{modification}]"
                else:
                    rendered += f"[{modification}]"
            if not settings.before_aa:
                rendered += amino_acid
        return rendered

    rendered = raw_features.map(convert)
    rendered = _apply_unrepresented_fixed_modifications(rendered, settings, target)
    if level == "ion":
        charges = canonical_features.astype("string").str.rsplit("/", n=1).str[-1]
        rendered = rendered.astype("string") + "/" + charges
    return rendered.astype("string")


def _count_sequence_characters(value: str, *, isalpha: bool, isupper: bool) -> int:
    return len(_strip_sequence(value, isalpha=isalpha, isupper=isupper))


def _strip_sequence(value: str, *, isalpha: bool, isupper: bool) -> str:
    if isalpha and isupper:
        return "".join(
            character for character in value if character.isalpha() and character.isupper()
        )
    if isalpha:
        return "".join(character for character in value if character.isalpha())
    if isupper:
        return "".join(character for character in value if character.isupper())
    return value


def _apply_unrepresented_fixed_modifications(
    features: pd.Series,
    settings: ModificationParserSettings,
    target: Any,
) -> pd.Series:
    """Add fixed residue modifications absent from a tool's result notation."""
    parameters = read_search_parameters(target)
    if parameters is None:
        return features
    represented_names = set(settings.modification_dict.values())
    result = features
    for fixed_modification in parameters.fixed_mods:
        match = re.fullmatch(r"(?P<targets>[A-Z]+)\[(?P<name>[^]]+)\]", fixed_modification.name)
        if match is None or match.group("name") in represented_names:
            continue
        targets = set(match.group("targets"))
        name = match.group("name")
        result = result.map(
            lambda feature: _add_fixed_residue_modification(
                str(feature),
                targets=targets,
                name=name,
            )
        )
    return result


def _add_fixed_residue_modification(
    feature: str,
    *,
    targets: set[str],
    name: str,
) -> str:
    charge = ""
    sequence = feature
    if "/" in feature:
        sequence, suffix = feature.rsplit("/", 1)
        charge = f"/{suffix}"
    rendered = ""
    in_brackets = False
    for character in sequence:
        rendered += character
        if character == "[":
            in_brackets = True
        elif character == "]":
            in_brackets = False
        elif not in_brackets and character in targets:
            rendered += f"[{name}]"
    return rendered + charge


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
