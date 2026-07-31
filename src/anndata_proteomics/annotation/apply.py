"""Join an external sample-annotation table onto the ``obs`` axis."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from anndata import AnnData
from loguru import logger
from mudata import MuData
from numpy.typing import NDArray

from anndata_proteomics.annotation._sanitize import sanitize_columns
from anndata_proteomics.annotation.loader import AnnotationTable

_MAX_REPORTED = 5


def annotate_obs(obj: AnnData | MuData, annotation: AnnotationTable) -> AnnData | MuData:
    """Join sample annotations onto ``obj``'s ``obs`` axis in place.

    Raises ``ValueError`` if no observation matches any annotation record.
    Partial mismatches are logged as warnings.
    """
    match_on = annotation.match_on

    holders: list[AnnData | MuData] = [obj]
    if isinstance(obj, MuData):
        holders += list(obj.mod.values())

    primary_keys = _obs_keys(obj, match_on)
    frame, matched_key_field = _matching_annotation_frame(annotation, primary_keys)
    in_table = primary_keys.isin(frame.index)
    n_matched = int(in_table.sum())
    if n_matched == 0:
        raise ValueError(
            f"no obs rows matched any annotation record on match_on={match_on!r} "
            f"(key_field={matched_key_field!r}). "
            f"first obs keys: {list(primary_keys[:_MAX_REPORTED])}; "
            f"first record keys: {list(frame.index[:_MAX_REPORTED])}"
        )

    cols_added: list[str] = []
    for holder in holders:
        cols_added = _join_obs_frame(
            _obs_frame(holder),
            _obs_keys(holder, match_on),
            frame,
        )

    _warn_on_mismatch(
        primary_keys,
        in_table,
        frame,
        matched_key_field=matched_key_field,
    )
    _record_provenance(
        obj,
        annotation,
        cols_added,
        n_matched,
        matched_key_field=matched_key_field,
    )
    logger.info(
        f"annotated obs: +{len(cols_added)} column(s) {cols_added}, "
        f"{n_matched}/{len(primary_keys)} rows matched"
    )
    return obj


def _matching_annotation_frame(
    annotation: AnnotationTable,
    obs_keys: pd.Index,
) -> tuple[pd.DataFrame, str]:
    """Try the primary sample identifier, then declared exact aliases."""
    key_fields = [annotation.key_field]
    alias_field = f"{annotation.key_field}_alias"
    if alias_field in annotation.samples:
        key_fields.append(alias_field)
    aliases_field = f"{annotation.key_field}_aliases"
    if aliases_field in annotation.samples:
        key_fields.append(aliases_field)

    selected = _build_annotation_frame(annotation, key_fields[0], key_fields)
    for key_field in key_fields:
        candidate = _build_annotation_frame(annotation, key_field, key_fields)
        if obs_keys.isin(candidate.index).any():
            return candidate, key_field
        selected = candidate
    return selected, key_fields[-1]


def _build_annotation_frame(
    annotation: AnnotationTable,
    join_field: str,
    identifier_fields: list[str],
) -> pd.DataFrame:
    """Index annotation rows by their string join value and sanitize columns."""
    frame = annotation.samples.copy()
    if join_field.endswith("_aliases"):
        frame = frame.explode(join_field)
    frame = frame.loc[frame[join_field].notna()].copy()
    frame[join_field] = frame[join_field].astype(str)
    if frame[join_field].duplicated().any():
        duplicates = sorted(frame[join_field][frame[join_field].duplicated()].unique())
        raise ValueError(f"duplicate {join_field!r} values in annotation table: {duplicates}")
    frame.index = pd.Index(frame[join_field], name=join_field)
    frame = frame.drop(columns=identifier_fields)
    frame.columns = sanitize_columns(list(frame.columns))
    return frame


def _obs_frame(holder: AnnData | MuData) -> pd.DataFrame:
    """Return an in-memory ``obs`` frame, rejecting a backed/lazy observation axis."""
    obs = holder.obs
    if not isinstance(obs, pd.DataFrame):
        raise TypeError("sample annotation requires an in-memory obs DataFrame")
    return obs


def _obs_keys(holder: AnnData | MuData, match_on: str) -> pd.Index:
    """Return string join keys for one object's observation axis."""
    if match_on == "index":
        return pd.Index(holder.obs_names, dtype="object").astype(str)
    obs = _obs_frame(holder)
    if match_on not in obs.columns:
        raise ValueError(
            f"match_on column {match_on!r} not found in obs columns: {list(obs.columns)}"
        )
    return pd.Index(obs[match_on].astype(str))


def _join_obs_frame(
    obs: pd.DataFrame,
    keys: pd.Index,
    annotation: pd.DataFrame,
) -> list[str]:
    """Assign annotation columns onto ``obs`` aligned by join keys."""
    overlap = [column for column in annotation.columns if column in obs.columns]
    if overlap:
        raise ValueError(f"annotation columns already present in obs: {overlap}")
    aligned = annotation.reindex(keys)
    for column in annotation.columns:
        obs[column] = aligned[column].to_numpy()
    return list(annotation.columns)


def _warn_on_mismatch(
    keys: pd.Index,
    in_table: NDArray[np.bool_],
    annotation: pd.DataFrame,
    *,
    matched_key_field: str,
) -> None:
    n_unmatched = int((~in_table).sum())
    if n_unmatched:
        logger.warning(f"{n_unmatched}/{len(keys)} obs rows had no matching annotation record")
    if matched_key_field.endswith("_aliases"):
        return
    key_set = set(keys)
    records_unmatched = [key for key in annotation.index if key not in key_set]
    if records_unmatched:
        shown = records_unmatched[:_MAX_REPORTED]
        tail = " …" if len(records_unmatched) > _MAX_REPORTED else ""
        logger.warning(
            f"{len(records_unmatched)} annotation record(s) matched no obs row: {shown}{tail}"
        )


def _record_provenance(
    obj: AnnData | MuData,
    annotation: AnnotationTable,
    cols_added: list[str],
    n_matched: int,
    *,
    matched_key_field: str,
) -> None:
    """Append lightweight annotation provenance under APB's ``uns`` namespace."""
    entry = {
        "source": str(annotation.source) if annotation.source else None,
        "source_format": annotation.source.suffix.lower().lstrip(".")
        if annotation.source
        else None,
        "match_on": annotation.match_on,
        "key_field": matched_key_field,
        "obs_columns_added": list(cols_added),
        "n_obs_matched": n_matched,
    }
    namespace = dict(obj.uns.get("anndata_proteomics", {}))
    existing = json.loads(namespace.get("obs_annotations_json", "[]"))
    existing.append(entry)
    namespace["obs_annotations_json"] = json.dumps(existing)
    obj.uns["anndata_proteomics"] = namespace
