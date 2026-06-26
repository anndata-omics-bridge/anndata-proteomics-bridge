"""Annotate the protein layer from FASTA file(s), stored under ``varm['fasta']``.

The feature-axis counterpart to :func:`annotation.apply.annotate_obs`. It builds
the prolfquapp-style protein-annotation frame (:func:`fasta.annotation.fasta_to_dataframe`)
and attaches it as a var-aligned DataFrame at ``varm['fasta']`` of the **protein**
layer only — a standalone protein-level AnnData, or the ``protein`` modality of a
MuData. The ``fasta`` namespace makes the source self-evident, so a FASTA-derived
quantity is understood to be the *theoretical* in-silico count, not an observed one
(e.g. ``nr_peptides`` is the in-silico digest count; the other columns keep their
bare prolfquapp names).

The join is keyed on the leading accession of the protein group (prolfquapp's
``cleanID``: the first ``;``-separated token, UniProt-middle-extracted), matched
against the FASTA ``proteinname``. The peptide count uses the *actual* digestion
enzyme: when the object was converted with a search-parameters file, the enzyme /
peptide-length bounds stored under ``uns['anndata_proteomics']['search_parameters']``
drive the cleavage rule, rather than assuming trypsin. A ``cleavage`` argument (or the
CLI ``--cleavage`` flag) overrides that for objects converted without parameters.
"""

from __future__ import annotations

import json
import re
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd
from loguru import logger

from anndata_proteomics.annotation._sanitize import sanitize_columns
from anndata_proteomics.fasta.annotation import (
    CleavageRule,
    fasta_to_dataframe,
    resolve_cleavage,
)
from anndata_proteomics.fasta.parser import FastaSource
from anndata_proteomics.params.anndata_io import read_search_parameters

_MAX_REPORTED = 5
_SCHEMA_VERSION = "0.1"
_DEFAULT_MIN_LENGTH = 7
_DEFAULT_MAX_LENGTH = 30
_VARM_KEY = "fasta"  # the varm slot the annotation DataFrame is stored under
_JOIN_KEY = "proteinname"  # the FASTA-frame column the var key is matched against
_PRT_PREFIX_RE = re.compile(r"^prt:")
_UNIPROT_MIDDLE_RE = re.compile(r".+\|(.+)\|.*")


def annotate_var_from_fasta(
    obj: Any,
    fasta_sources: FastaSource | Iterable[FastaSource],
    *,
    match_on: str = "Protein_Group",
    is_uniprot: bool = True,
    decoy_pattern: str = "^REV_|^rev_",
    cleavage: str | CleavageRule | None = None,
    min_length: int | None = None,
    max_length: int | None = None,
    include_sequence: bool = False,
    columns: Iterable[str] | None = None,
) -> Any:
    """Attach FASTA-derived annotation at the protein layer's ``varm['fasta']``, in place.

    *obj* is a protein-level AnnData or a MuData (whose ``protein`` modality is
    annotated). *match_on* names the ``var`` column holding the protein group
    (``"index"`` uses ``var_names``). *cleavage* / *min_length* / *max_length*
    override the enzyme and peptide-length bounds otherwise read from the stored
    search parameters. *columns*, if given, restricts which FASTA columns are
    stored. Returns *obj*.

    The annotation is written as a var-aligned DataFrame at ``varm['fasta']``
    (rows in ``var`` order; unmatched proteins are NaN). Raises ValueError when
    *obj* is not a protein layer, when ``varm['fasta']`` already exists, or when
    no var row matches any FASTA record. Partial mismatches are logged as
    warnings, not raised — mirroring :func:`annotate_obs`.
    """
    target = _resolve_protein_target(obj)
    _ensure_varm_free(target)

    rule, enzyme_name, min_len, max_len = _resolve_digestion(
        target, cleavage, min_length, max_length
    )
    ann = fasta_to_dataframe(
        fasta_sources,
        decoy_pattern=decoy_pattern,
        is_uniprot=is_uniprot,
        cleavage=rule,
        min_length=min_len,
        max_length=max_len,
        include_sequence=include_sequence,
    )
    ann = _index_by_join_key(ann)

    keys = _var_join_keys(target, match_on, is_uniprot=is_uniprot)
    in_table = keys.isin(ann.index)
    n_matched = int(in_table.sum())
    if n_matched == 0:
        raise ValueError(
            f"no var rows matched any FASTA record on match_on={match_on!r} "
            f"(leading accession of the protein group). "
            f"first var keys: {list(keys[:_MAX_REPORTED])}; "
            f"first FASTA proteinnames: {list(ann.index[:_MAX_REPORTED])}"
        )

    frame = _build_varm_frame(target, keys, ann, columns)
    target.varm[_VARM_KEY] = frame
    cols_added = list(frame.columns)

    _warn_on_mismatch(keys, in_table, ann)
    _record_provenance(
        target,
        fasta_sources=fasta_sources,
        match_on=match_on,
        cols_added=cols_added,
        n_matched=n_matched,
        enzyme=enzyme_name,
        min_length=min_len,
        max_length=max_len,
    )
    logger.info(
        f"stored protein annotation in varm[{_VARM_KEY!r}]: "
        f"{len(cols_added)} column(s) {cols_added}, "
        f"{n_matched}/{len(keys)} rows matched (enzyme={enzyme_name})"
    )
    return obj


def _resolve_protein_target(obj: Any) -> Any:
    """Return the protein-level AnnData to annotate, or raise if *obj* is not one."""
    if hasattr(obj, "mod"):  # MuData: annotate the protein modality only
        if "protein" not in obj.mod:
            raise ValueError(f"MuData has no 'protein' modality; modalities: {list(obj.mod)}")
        return obj.mod["protein"]
    level = (obj.uns.get("anndata_proteomics") or {}).get("quantification_level")
    if level != "protein":
        raise ValueError(
            "FASTA var-annotation applies to the protein layer only; got "
            f"quantification_level={level!r}"
        )
    return obj


def _resolve_digestion(
    target: Any,
    cleavage: str | CleavageRule | None,
    min_length: int | None,
    max_length: int | None,
) -> tuple[CleavageRule, str, int, int]:
    """Resolve the cleavage rule and peptide-length bounds for the count.

    Explicit arguments win; otherwise values come from the stored search
    parameters; otherwise trypsin / 7 / 30 with a warning (no enzyme on hand).
    """
    params = read_search_parameters(target)
    enzyme_from_params = params.enzyme if params else None

    if cleavage is None and enzyme_from_params is None:
        logger.warning(
            "no enzyme in search parameters and no cleavage override; "
            "using Trypsin for the peptide count"
        )
    spec = cleavage if cleavage is not None else enzyme_from_params
    rule, enzyme_name = resolve_cleavage(spec)

    min_len = _first_present(
        min_length, params.min_peptide_length if params else None, _DEFAULT_MIN_LENGTH
    )
    max_len = _first_present(
        max_length, params.max_peptide_length if params else None, _DEFAULT_MAX_LENGTH
    )
    return rule, enzyme_name, min_len, max_len


def _first_present(*values: int | None) -> int:
    """First non-None value (the last one is always a concrete default)."""
    for value in values:
        if value is not None:
            return value
    raise AssertionError("the final fallback must be non-None")  # pragma: no cover


def _index_by_join_key(ann: pd.DataFrame) -> pd.DataFrame:
    """Index the FASTA frame by ``proteinname``, dropping duplicate keys (first wins)."""
    if _JOIN_KEY not in ann.columns:  # pragma: no cover - fasta_to_dataframe always emits it
        raise ValueError(f"FASTA frame is missing the {_JOIN_KEY!r} join column")
    duplicated = ann[_JOIN_KEY].duplicated()
    if duplicated.any():
        dups = sorted(ann[_JOIN_KEY][duplicated].unique())[:_MAX_REPORTED]
        logger.warning(
            f"{int(duplicated.sum())} duplicate {_JOIN_KEY!r} value(s) in FASTA; "
            f"keeping first occurrence: {dups}"
        )
        ann = ann.loc[~duplicated]
    return ann.set_index(_JOIN_KEY)


def _var_join_keys(target: Any, match_on: str, *, is_uniprot: bool) -> pd.Index:
    """Leading-accession join keys for the protein ``var`` axis (prolfquapp cleanID)."""
    if match_on == "index":
        raw = pd.Index(target.var_names, dtype="object").astype(str)
    else:
        var = target.var
        if match_on not in var.columns:
            raise ValueError(
                f"match_on column {match_on!r} not found in var columns: {list(var.columns)}"
            )
        raw = pd.Index(var[match_on].astype(str))
    return pd.Index([_leading_accession(v, is_uniprot=is_uniprot) for v in raw])


def _leading_accession(group_value: str, *, is_uniprot: bool) -> str:
    """First accession of a protein group, matched to the FASTA ``proteinname`` form.

    Strips a ``prt:`` modality prefix, takes the first ``;``-separated token, and
    (for UniProt) extracts the middle ``db|ACC|NAME`` field — mirroring how
    ``fasta_to_dataframe`` derives ``proteinname``.
    """
    token = _PRT_PREFIX_RE.sub("", str(group_value)).strip().split(";")[0].strip()
    if is_uniprot:
        match = _UNIPROT_MIDDLE_RE.match(token)
        return match.group(1) if match else token
    return token


def _ensure_varm_free(target: Any) -> None:
    """Raise if ``varm['fasta']`` already exists (don't silently clobber a prior run)."""
    if _VARM_KEY in target.varm:
        raise ValueError(
            f"varm[{_VARM_KEY!r}] already present on the protein layer; "
            "delete it before re-annotating"
        )


def _build_varm_frame(
    target: Any,
    keys: pd.Index,
    ann: pd.DataFrame,
    columns: Iterable[str] | None,
) -> pd.DataFrame:
    """Var-aligned FASTA annotation DataFrame (rows in ``var`` order; unmatched → NaN)."""
    available = list(ann.columns)
    selected = [c for c in columns if c in available] if columns is not None else available
    aligned = ann.reindex(keys)[selected]  # rows in var order; unmatched keys → NaN
    aligned.columns = sanitize_columns(selected)
    aligned.index = pd.Index(target.var_names)  # varm is aligned to the var axis
    return aligned


def _warn_on_mismatch(keys: pd.Index, in_table: pd.Series, ann: pd.DataFrame) -> None:
    n_unmatched = int((~in_table).sum())
    if n_unmatched:
        logger.warning(f"{n_unmatched}/{len(keys)} var rows had no matching FASTA record")
    key_set = set(keys)
    records_unmatched = [k for k in ann.index if k not in key_set]
    if records_unmatched:
        shown = records_unmatched[:_MAX_REPORTED]
        tail = " …" if len(records_unmatched) > _MAX_REPORTED else ""
        logger.info(f"{len(records_unmatched)} FASTA record(s) matched no var row: {shown}{tail}")


def _record_provenance(
    target: Any,
    *,
    fasta_sources: FastaSource | Iterable[FastaSource],
    match_on: str,
    cols_added: list[str],
    n_matched: int,
    enzyme: str,
    min_length: int,
    max_length: int,
) -> None:
    """Append an entry under ``uns['anndata_proteomics']['var_annotations_json']``.

    Stored as a JSON string (mirroring ``obs_annotations_json``) so h5py can
    serialise it. Records the FASTA source(s) and the enzyme / length bounds the
    peptide count actually used.
    """
    entry = {
        "schema_version": _SCHEMA_VERSION,
        "source": "fasta",
        "destination": f"varm[{_VARM_KEY!r}]",
        "fasta_sources": _describe_sources(fasta_sources),
        "match_on": match_on,
        "columns": list(cols_added),
        "n_var_matched": n_matched,
        "cleavage_enzyme": enzyme,
        "min_peptide_length": min_length,
        "max_peptide_length": max_length,
    }
    namespace = dict(target.uns.get("anndata_proteomics", {}))
    existing = json.loads(namespace.get("var_annotations_json", "[]"))
    existing.append(entry)
    namespace["var_annotations_json"] = json.dumps(existing)
    target.uns["anndata_proteomics"] = namespace


def _describe_sources(sources: FastaSource | Iterable[FastaSource]) -> list[str]:
    """Human-readable list of FASTA sources.

    Paths are recorded as strings; inline FASTA text (the parser's
    string-content path) as ``<inline-fasta>``; open streams as ``<stream>``.
    """
    single = isinstance(sources, str | Path) or hasattr(sources, "read")
    items = [sources] if single else list(sources)
    out: list[str] = []
    for item in items:
        if isinstance(item, Path):
            out.append(str(item))
        elif isinstance(item, str):
            out.append("<inline-fasta>" if _is_inline_fasta(item) else item)
        else:
            out.append("<stream>")
    return out


def _is_inline_fasta(text: str) -> bool:
    """True when *text* is FASTA content rather than a path (matches the parser)."""
    return "\n" in text or text.lstrip().startswith(">")
