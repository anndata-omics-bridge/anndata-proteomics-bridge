"""Apply an AnnotationSpec to an AnnData/MuData object (obs axis).

The join is keyed on the run/file identifier: each ``obs.samples`` record is matched against
``obs_names`` (default) or a named ``obs`` column, and the record's remaining fields are added
as ``obs`` columns. For a MuData the obs axis is shared, so every modality is annotated
(the global ``mdata.obs`` and each ``mdata.mod[m].obs``).
"""

from __future__ import annotations

import json
from typing import Any

import pandas as pd
from loguru import logger

from anndata_proteomics.annotation._sanitize import sanitize_columns
from anndata_proteomics.annotation.schema import AnnotationSpec, ObsAnnotation

_MAX_REPORTED = 5


def annotate_obs(obj: Any, spec: AnnotationSpec) -> Any:
    """Join ``spec.obs.samples`` onto ``obj``'s obs axis, in place. Returns ``obj``.

    Raises ValueError if no obs row matches any record (a wrong axis/key is almost always
    a misconfiguration). Partial mismatches are logged as warnings, not raised.
    """
    ann = _build_annotation_frame(spec.obs)
    match_on = spec.obs.match_on

    holders = [obj]
    if hasattr(obj, "mod"):  # MuData: shared obs axis — annotate global obs and every modality
        holders += [obj.mod[name] for name in obj.mod]

    primary_keys = _obs_keys(obj, match_on)
    in_table = primary_keys.isin(ann.index)
    n_matched = int(in_table.sum())
    if n_matched == 0:
        raise ValueError(
            f"no obs rows matched any annotation record on match_on={match_on!r} "
            f"(key_field={spec.obs.key_field!r}). "
            f"first obs keys: {list(primary_keys[:_MAX_REPORTED])}; "
            f"first record keys: {list(ann.index[:_MAX_REPORTED])}"
        )

    cols_added: list[str] = []
    for holder in holders:
        cols_added = _join_obs_frame(holder.obs, _obs_keys(holder, match_on), ann)

    _warn_on_mismatch(primary_keys, in_table, ann)
    _record_provenance(obj, spec, cols_added, n_matched)
    logger.info(
        f"annotated obs: +{len(cols_added)} column(s) {cols_added}, "
        f"{n_matched}/{len(primary_keys)} rows matched"
    )
    return obj


def _build_annotation_frame(obs_spec: ObsAnnotation) -> pd.DataFrame:
    """Records → DataFrame indexed by the (string) join value, with sanitised obs columns."""
    ann = pd.DataFrame(list(obs_spec.samples))
    key = obs_spec.key_field
    ann[key] = ann[key].astype(str)
    if ann[key].duplicated().any():
        dups = sorted(ann[key][ann[key].duplicated()].unique())
        raise ValueError(f"duplicate {key!r} values in obs.samples: {dups}")
    ann = ann.set_index(key)
    ann.columns = sanitize_columns(list(ann.columns))
    return ann


def _obs_keys(holder: Any, match_on: str) -> pd.Index:
    """The string join keys for ``holder``'s obs axis."""
    if match_on == "index":
        return pd.Index(holder.obs_names, dtype="object").astype(str)
    obs = holder.obs
    if match_on not in obs.columns:
        raise ValueError(
            f"match_on column {match_on!r} not found in obs columns: {list(obs.columns)}"
        )
    return pd.Index(obs[match_on].astype(str))


def _join_obs_frame(obs: pd.DataFrame, keys: pd.Index, ann: pd.DataFrame) -> list[str]:
    """Assign each annotation column onto ``obs`` (in place), aligned by ``keys``."""
    overlap = [c for c in ann.columns if c in obs.columns]
    if overlap:
        raise ValueError(f"annotation columns already present in obs: {overlap}")
    aligned = ann.reindex(keys)  # rows in obs order; unmatched keys → NaN
    for col in ann.columns:
        obs[col] = aligned[col].to_numpy()
    return list(ann.columns)


def _warn_on_mismatch(keys: pd.Index, in_table: pd.Series, ann: pd.DataFrame) -> None:
    n_unmatched = int((~in_table).sum())
    if n_unmatched:
        logger.warning(f"{n_unmatched}/{len(keys)} obs rows had no matching annotation record")
    key_set = set(keys)
    records_unmatched = [k for k in ann.index if k not in key_set]
    if records_unmatched:
        shown = records_unmatched[:_MAX_REPORTED]
        tail = " …" if len(records_unmatched) > _MAX_REPORTED else ""
        logger.warning(
            f"{len(records_unmatched)} annotation record(s) matched no obs row: {shown}{tail}"
        )


def _record_provenance(
    obj: Any, spec: AnnotationSpec, cols_added: list[str], n_matched: int
) -> None:
    """Append a provenance entry under ``uns['anndata_proteomics']['obs_annotations_json']``.

    Stored as a JSON string (mirroring ``rule_json`` in assemble.py) so h5py can serialise it.
    """
    entry = {
        "schema_version": spec.schema_version,
        "match_on": spec.obs.match_on,
        "key_field": spec.obs.key_field,
        "obs_columns_added": list(cols_added),
        "n_obs_matched": n_matched,
    }
    namespace = dict(obj.uns.get("anndata_proteomics", {}))
    existing = json.loads(namespace.get("obs_annotations_json", "[]"))
    existing.append(entry)
    namespace["obs_annotations_json"] = json.dumps(existing)
    obj.uns["anndata_proteomics"] = namespace
