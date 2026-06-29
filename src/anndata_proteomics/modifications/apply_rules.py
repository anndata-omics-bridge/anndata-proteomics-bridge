"""Token extraction and mapping from vendor modified sequences.

Takes a vendor-specific modified sequence (e.g. ``"PEPM[15.9949]TIDE"`` or
``"_(ac)PEPTIDEM(ox)_"``) plus a modification rule (regex + map entries)
and produces a :class:`ModifiedSequence` with localized
:class:`ModificationOccurrence`\\s and a ProForma rendering.

Map lookup uses the tuple ``(mass_delta, target, position)`` per the plan,
not mass alone — so e.g. Acetyl-Nterm and Acetyl-K with the same mass
remain distinguishable.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import Iterable

from anndata_proteomics.modifications.model import (
    ModificationOccurrence,
    ModifiedSequence,
)
from anndata_proteomics.modifications.proforma import render_proforma


@dataclass(frozen=True)
class MapEntry:
    """One ``[[modifications.map]]`` entry from a parsing-rule TOML."""

    token: str
    name: str
    accession: str | None = None
    target: list[str] | None = None  # allowed residues/termini (from the Unimod registry)
    position: str | None = "Anywhere"
    mass_delta: float | None = None


@dataclass(frozen=True)
class ModificationRule:
    """Parsed ``[modifications]`` section."""

    source_column: str
    token_pattern: str
    token_position: str  # "before_residue" | "after_residue"
    case_sensitive: bool = False
    unknown_policy: str = "preserve"  # "preserve" | "drop" | "error"
    sequence_column: str | None = None
    output_column: str = "proforma_sequence"
    entries: tuple[MapEntry, ...] = ()


_MASS_TOLERANCE = 1e-3
_TERM_MARKERS = {"_", "-", "."}
_TERMINUS_TARGETS = {
    "N-term",
    "C-term",
    "Protein N-term",
    "Protein C-term",
    "Peptide N-term",
    "Peptide C-term",
}


def _target_matches(
    entry_target: list[str] | None, adjacent_residue: str | None, position: str
) -> bool:
    """Decide whether an entry's allowed ``target``s are compatible with a token's context.

    ``entry_target`` is the list of residues/termini the modification may sit on.
    Terminal targets (e.g. ``"N-term"``) match when the token's position is a
    corresponding terminus; residue targets (``"M"``, ``"C"``, …) match by
    amino-acid identity. Matches if ANY listed target is compatible. An empty/None
    target matches anything.
    """
    if not entry_target:
        return True
    for target in entry_target:
        if target in _TERMINUS_TARGETS:
            if target.endswith(position) or target == position:
                return True
        elif adjacent_residue is not None and target == adjacent_residue:
            return True
    return False


def _parse_mass(raw: str) -> float | None:
    cleaned = raw.strip().lstrip("+")
    try:
        return float(cleaned)
    except ValueError:
        return None


def _match_entry(
    entries: Iterable[MapEntry],
    raw_token: str,
    adjacent_residue: str | None,
    position: str,
    case_sensitive: bool,
) -> MapEntry | None:
    """Pick the best map entry for a vendor token.

    For numeric tokens the lookup key is the tuple
    ``(mass_delta, target, position)``. For non-numeric tokens (e.g.
    ``"ox"``, ``"ac"``) the exact token string is used as the fallback.
    """
    parsed_mass = _parse_mass(raw_token)
    if parsed_mass is not None:
        for entry in entries:
            if entry.mass_delta is None:
                continue
            if not math.isclose(entry.mass_delta, parsed_mass, abs_tol=_MASS_TOLERANCE):
                continue
            if entry.position and entry.position != position:
                continue
            if not _target_matches(entry.target, adjacent_residue, position):
                continue
            return entry

    cmp_token = raw_token if case_sensitive else raw_token.lower()
    for entry in entries:
        entry_token = entry.token if case_sensitive else entry.token.lower()
        if entry_token != cmp_token:
            continue
        if entry.position and entry.position != position:
            continue
        if not _target_matches(entry.target, adjacent_residue, position):
            continue
        return entry
    return None


@dataclass
class _PendingToken:
    raw_token: str
    residue_index: int | None  # None for terminal
    position: str  # "Anywhere" | "N-term" | "C-term"
    adjacent_residue: str | None


def _strip_terminal_markers(seq: str) -> str:
    """Drop leading/trailing terminal markers (``_``, ``-``, ``.``)."""
    while seq and seq[0] in _TERM_MARKERS:
        seq = seq[1:]
    while seq and seq[-1] in _TERM_MARKERS:
        seq = seq[:-1]
    return seq


def _tokenize(
    seq: str, pattern: re.Pattern[str], token_position: str
) -> tuple[list[str], list[_PendingToken]]:
    """Walk regex matches, building the stripped residue sequence and the
    position-classified pending tokens (N-term / C-term / before-residue / Anywhere).
    """
    stripped_chars: list[str] = []
    pending: list[_PendingToken] = []
    cursor = 0

    for match in pattern.finditer(seq):
        for ch in seq[cursor : match.start()]:
            if ch.isalpha():
                stripped_chars.append(ch)

        groups = [g for g in match.groups() if g is not None]
        raw_token = groups[0] if groups else match.group(0)

        if not stripped_chars and match.start() == 0:
            position = "N-term"
            residue_idx = None
            adjacent = seq[match.end() : match.end() + 1] or None
        elif match.end() == len(seq) and token_position != "before_residue":
            position = "C-term"
            residue_idx = None
            adjacent = stripped_chars[-1] if stripped_chars else None
        elif token_position == "before_residue":
            next_residue = seq[match.end() : match.end() + 1]
            if next_residue.isalpha():
                stripped_chars.append(next_residue)
                cursor = match.end() + 1
                position = "Anywhere"
                residue_idx = len(stripped_chars) - 1
                adjacent = next_residue
                pending.append(_PendingToken(raw_token, residue_idx, position, adjacent))
                continue
            position = "Anywhere"
            residue_idx = len(stripped_chars) - 1 if stripped_chars else None
            adjacent = stripped_chars[-1] if stripped_chars else None
        else:
            position = "Anywhere"
            residue_idx = len(stripped_chars) - 1 if stripped_chars else None
            adjacent = stripped_chars[-1] if stripped_chars else None

        pending.append(_PendingToken(raw_token, residue_idx, position, adjacent))
        cursor = match.end()

    for ch in seq[cursor:]:
        if ch.isalpha():
            stripped_chars.append(ch)

    return stripped_chars, pending


def _resolve_tokens(
    pending: list[_PendingToken], rule: ModificationRule, stripped: str
) -> tuple[list[ModificationOccurrence], dict[int, str], list[str]]:
    """Match pending tokens to map entries; apply ``unknown_policy`` and record
    terminal/residue indices for unresolved tokens.
    """
    occurrences: list[ModificationOccurrence] = []
    unknown_tokens: dict[int, str] = {}
    unknown_list: list[str] = []

    for tok in pending:
        entry = _match_entry(
            rule.entries,
            tok.raw_token,
            tok.adjacent_residue,
            tok.position,
            rule.case_sensitive,
        )
        if entry is not None:
            occurrences.append(
                ModificationOccurrence(
                    name=entry.name,
                    accession=entry.accession,
                    target_residue=tok.adjacent_residue,
                    sequence_index=tok.residue_index,
                    position=tok.position,
                    mass_delta=entry.mass_delta,
                    source_token=tok.raw_token,
                )
            )
            continue
        if rule.unknown_policy == "error":
            raise ValueError(f"unknown modification token: {tok.raw_token!r}")
        if rule.unknown_policy == "drop":
            continue
        unknown_list.append(tok.raw_token)
        if tok.position == "N-term":
            unknown_tokens[-1] = tok.raw_token
        elif tok.position == "C-term":
            unknown_tokens[len(stripped)] = tok.raw_token
        elif tok.residue_index is not None:
            unknown_tokens[tok.residue_index] = tok.raw_token

    return occurrences, unknown_tokens, unknown_list


def apply_rule(modified_sequence: str, rule: ModificationRule) -> ModifiedSequence:
    """Normalize a vendor modified sequence via ``rule``: strip → tokenize → resolve → render."""
    pattern = re.compile(rule.token_pattern)
    seq = _strip_terminal_markers(modified_sequence)
    stripped_chars, pending = _tokenize(seq, pattern, rule.token_position)
    stripped = "".join(stripped_chars)
    occurrences, unknown_tokens, unknown_list = _resolve_tokens(pending, rule, stripped)
    proforma = render_proforma(stripped, occurrences, unknown_tokens=unknown_tokens)
    return ModifiedSequence(
        stripped_sequence=stripped,
        proforma_sequence=proforma,
        occurrences=occurrences,
        source_sequence=modified_sequence,
        unknown_tokens=unknown_list,
    )
