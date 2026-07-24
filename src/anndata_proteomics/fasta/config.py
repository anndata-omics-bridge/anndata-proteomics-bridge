"""Typed configuration for FASTA decoy and contaminant identifiers.

Patterns are inferred conservatively from raw FASTA identifier tokens unless
the caller supplies an explicit tuple.  An explicit empty tuple disables the
classification.  Classification is annotation only: it never filters FASTA
records or quantified features.
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from functools import cache
from typing import Literal

from pydantic import BaseModel, ConfigDict, field_validator

DEFAULT_DECOY_CANDIDATES = (
    r"^REV_",
    r"^rev_",
    r"^DECOY_",
    r"^decoy_",
    r"^XXX_",
    r"^reverse_",
)
DEFAULT_CONTAMINANT_CANDIDATES = (
    r"^zz(?:\||_)",
    r"^CON__",
    r"^CON_",
    r"^Cont_",
    r"^contam_",
)


class FastaConfig(BaseModel):
    """Input policy for resolving FASTA identifier classifications."""

    model_config = ConfigDict(frozen=True)

    decoy_patterns: tuple[str, ...] | None = None
    contaminant_patterns: tuple[str, ...] | None = None
    decoy_candidates: tuple[str, ...] = DEFAULT_DECOY_CANDIDATES
    contaminant_candidates: tuple[str, ...] = DEFAULT_CONTAMINANT_CANDIDATES

    @field_validator(
        "decoy_patterns",
        "contaminant_patterns",
        "decoy_candidates",
        "contaminant_candidates",
    )
    @classmethod
    def _valid_regexes(cls, patterns: tuple[str, ...] | None) -> tuple[str, ...] | None:
        if patterns is None:
            return None
        for pattern in patterns:
            try:
                re.compile(pattern)
            except re.error as error:
                raise ValueError(f"invalid FASTA identifier regex {pattern!r}: {error}") from error
        return patterns


class ResolvedPatternConfig(BaseModel):
    """Effective patterns for one identifier class and their observed counts."""

    model_config = ConfigDict(frozen=True)

    patterns: tuple[str, ...]
    source: Literal["explicit", "inferred", "none"]
    match_counts: dict[str, int]


class ResolvedFastaConfig(BaseModel):
    """Resolved decoy/contaminant policy for one supplied FASTA database."""

    model_config = ConfigDict(frozen=True)

    schema_version: Literal["0.1"] = "0.1"
    decoy: ResolvedPatternConfig
    contaminant: ResolvedPatternConfig
    n_fasta_ids: int


class FastaConfigAccumulator:
    """Resolve a :class:`FastaConfig` while FASTA identifiers stream past."""

    def __init__(self, config: FastaConfig) -> None:
        self._config = config
        self._decoy_tested = (
            config.decoy_candidates if config.decoy_patterns is None else config.decoy_patterns
        )
        self._contaminant_tested = (
            config.contaminant_candidates
            if config.contaminant_patterns is None
            else config.contaminant_patterns
        )
        self._decoy_counts = dict.fromkeys(self._decoy_tested, 0)
        self._contaminant_counts = dict.fromkeys(self._contaminant_tested, 0)
        self._n_fasta_ids = 0

    def observe(self, fasta_id: str) -> None:
        """Count all configured candidate matches for one raw FASTA ID."""
        self._n_fasta_ids += 1
        for pattern, regex in zip(
            self._decoy_tested,
            _compile_patterns(self._decoy_tested),
            strict=True,
        ):
            if regex.search(fasta_id):
                self._decoy_counts[pattern] += 1
        for pattern, regex in zip(
            self._contaminant_tested,
            _compile_patterns(self._contaminant_tested),
            strict=True,
        ):
            if regex.search(fasta_id):
                self._contaminant_counts[pattern] += 1

    def resolve(self) -> ResolvedFastaConfig:
        """Return immutable effective patterns and the complete count audit."""
        return ResolvedFastaConfig(
            decoy=_resolve_patterns(self._config.decoy_patterns, self._decoy_counts),
            contaminant=_resolve_patterns(
                self._config.contaminant_patterns,
                self._contaminant_counts,
            ),
            n_fasta_ids=self._n_fasta_ids,
        )


def resolve_fasta_config(
    fasta_ids: Iterable[str],
    config: FastaConfig,
) -> ResolvedFastaConfig:
    """Resolve *config* from a materialized collection of raw FASTA IDs."""
    accumulator = FastaConfigAccumulator(config)
    for fasta_id in fasta_ids:
        accumulator.observe(fasta_id)
    return accumulator.resolve()


def matches_any(value: str, patterns: tuple[str, ...]) -> bool:
    """Return whether *value* matches at least one compiled regex."""
    return any(regex.search(value) for regex in _compile_patterns(patterns))


def _resolve_patterns(
    explicit: tuple[str, ...] | None,
    counts: dict[str, int],
) -> ResolvedPatternConfig:
    if explicit is not None:
        return ResolvedPatternConfig(
            patterns=explicit,
            source="explicit",
            match_counts=counts,
        )
    inferred = tuple(pattern for pattern, count in counts.items() if count)
    return ResolvedPatternConfig(
        patterns=inferred,
        source="inferred" if inferred else "none",
        match_counts=counts,
    )


@cache
def _compile_patterns(patterns: tuple[str, ...]) -> tuple[re.Pattern[str], ...]:
    return tuple(re.compile(pattern) for pattern in patterns)
