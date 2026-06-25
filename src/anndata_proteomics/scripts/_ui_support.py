"""Framework-agnostic logic behind the test-data browser GUI (ui_test_tool.py).

Kept free of any GUI/marimo import so it can be unit-tested and reused. It turns the ProteoBench
test-data index into a catalog and converts a dataset to AnnData/MuData. Conversion is **param-
driven**: the co-located param file is parsed for the software version, the version selects the rule
variant (folder), and the data columns must match that rule or conversion errors. See
TODO/TODO_ui_test_tool.md and the marimo-background-jobs skill.
"""

from __future__ import annotations

import csv
import re
import shlex
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from anndata_proteomics.converters.recognize import matches
from anndata_proteomics.params.anndata_io import (
    get_search_parameters_path,
    read_search_parameters,
    write_search_parameters,
)
from anndata_proteomics.params.registry import parse_params
from anndata_proteomics.rules.loader import resolve_rule_for_version
from anndata_proteomics.test_data import DOWNLOADED_DB, TEST_DATA_DIR

# Background-job conversions write their result + console.log into a per-run subdir here.
CONVERTED_DIR = TEST_DATA_DIR.parent / "logs" / "ui_converted"
CONVERTED_COLUMNS = [
    "run_name",
    "timestamp",
    "software_name",
    "software_version",
    "slug",
    "target",
    "status",
    "result_type",
    "nr_prec",
    "size_mb",
    "input_file_path",
    "param_path",
    "output_dir",
    "result_path",
    "log_path",
]

# Quantification levels, coarse to fine. Not every vendor exposes every level.
LEVELS = ["ion", "peptidoform", "peptide", "protein", "fragment"]
MUDATA = "mudata"
# Per-level var_names prefix so modalities don't collide on the global axis.
_PREFIX = {
    "fragment": "frg:",
    "ion": "ion:",
    "peptidoform": "pfm:",
    "peptide": "pep:",
    "protein": "prt:",
}
# Converting these targets on a large input is memory-heavy (the fragment explode); warn first.
_HEAVY_TARGETS = {"fragment", MUDATA}
_HEAVY_SIZE_MB = 100.0
_RUN_DIR_RE = re.compile(r"^(?P<timestamp>\d{8}T\d{6})_(?P<slug>[a-z0-9]+)_(?P<target>[a-z0-9_]+)$")


def software_slug(software_name: str) -> str:
    """Map a catalog ``software_name`` (e.g. "DIA-NN") to a parsing-rule vendor slug ("diann")."""
    return re.sub(r"[^a-z0-9]", "", software_name.lower())


def _dataset_path(input_file_path: str) -> Path:
    return TEST_DATA_DIR / "json_dir" / input_file_path


def _headers(path: Path) -> set[str]:
    """Read just the column names of a cached input (cheap; for convertibility checks)."""
    path = Path(path)
    if path.suffix == ".parquet":
        return set(pq.read_schema(path).names)
    return set(pd.read_csv(path, sep="\t", nrows=0).columns)


def param_path_for(input_file_path: str) -> Path | None:
    """The co-located ProteoBench param file (``param_0..*``) next to a dataset's input file."""
    candidates = sorted(_dataset_path(input_file_path).parent.glob("param_0.*"))
    return candidates[0] if candidates else None


def _param_version(param_path: Path | None, slug: str) -> str | None:
    """Software version parsed from the param file, or None if absent/unparseable."""
    if param_path is None:
        return None
    try:
        return parse_params(param_path, software=slug).software_version
    except Exception:  # noqa: BLE001 — any parse failure means "version unknown"
        return None


def select_rule(slug: str, level: str, version: str | None, headers: Iterable[str]):
    """Resolve the rule for (slug, level) at ``version`` and validate it against ``headers``.

    Raises ValueError if no rule variant covers the version, or if the file's columns don't match
    the version-selected rule (no silent fallback — the caller fixes the version / param file).
    """
    rule = resolve_rule_for_version(slug, level, version)
    if rule is None:
        raise ValueError(f"{slug} {level}: no rule covers software version {version!r}")
    header_set = set(headers)
    label_col = (
        rule.fragments.label_column
        if rule.fragments is not None and rule.fragments.label_strategy == "column"
        else None
    )
    if not matches(header_set, rule) or (label_col is not None and label_col not in header_set):
        raise ValueError(
            f"{slug} {level}: file columns don't match the rule for software version "
            f"{version!r} — verify the version / provide the right param file"
        )
    return rule


def convertible_levels(slug: str, version: str | None, headers: Iterable[str]) -> list[str]:
    """Levels whose version-selected rule both exists and matches this file's columns."""
    header_set = set(headers)
    out = []
    for level in LEVELS:
        try:
            select_rule(slug, level, version, header_set)
        except (LookupError, ValueError):
            continue
        out.append(level)
    return out


def available_targets(slug: str, version: str | None, headers: Iterable[str]) -> list[str]:
    """Convertible levels for this file/version, plus ``"mudata"`` when at least two resolve."""
    levels = convertible_levels(slug, version, headers)
    targets = list(levels)
    if len(levels) >= 2:
        targets.append(MUDATA)
    return targets


def load_catalog() -> pd.DataFrame:
    """Read the ProteoBench test-data index into a catalog DataFrame.

    Per-row columns added: ``size_mb``, ``slug``, ``param_path`` (co-located param file or ""),
    ``targets`` (tuple — convertible only when a param file gives a version whose rule matches the
    data columns), ``targets_str``. Only ``status == "ok"`` rows are kept. Empty when the cache
    index is absent (gitignored — regenerate via ``test_data_download/Makefile``).
    """
    if not DOWNLOADED_DB.exists():
        return pd.DataFrame(
            columns=[
                "software_name",
                "software_version",
                "nr_prec",
                "size_mb",
                "slug",
                "param_path",
                "targets",
                "targets_str",
                "input_file_path",
            ]
        )
    rows = []
    by_path: dict[str, tuple[tuple[str, ...], str]] = {}  # rel -> (targets, param_path)
    with open(DOWNLOADED_DB) as f:
        for row in csv.DictReader(f):
            if row.get("status") != "ok":
                continue
            slug = software_slug(row["software_name"])
            rel = row["input_file_path"]
            if rel not in by_path:
                param = param_path_for(rel)
                try:
                    if param is None:
                        targets: tuple[str, ...] = ()
                    else:
                        version = _param_version(param, slug)
                        headers = _headers(_dataset_path(rel))
                        targets = tuple(available_targets(slug, version, headers))
                except (OSError, ValueError):  # missing/unreadable cached file
                    targets = ()
                by_path[rel] = (targets, str(param) if param else "")
            targets, param_str = by_path[rel]
            size_bytes = (
                int(row["input_file_size_bytes"])
                if row.get("input_file_size_bytes", "").isdigit()
                else 0
            )
            rows.append(
                {
                    "software_name": row["software_name"],
                    "software_version": row.get("software_version", ""),
                    "nr_prec": int(row["nr_prec"]) if row.get("nr_prec", "").isdigit() else 0,
                    "size_mb": round(size_bytes / 1e6, 1),
                    "slug": slug,
                    "param_path": param_str,
                    "targets": targets,
                    "targets_str": ", ".join(targets) if targets else "—",
                    "input_file_path": rel,
                }
            )
    return pd.DataFrame(rows)


def filter_catalog(
    catalog: pd.DataFrame,
    *,
    target: str | None = None,
    software: str | None = None,
    max_size_mb: float | None = None,
) -> pd.DataFrame:
    """Apply the GUI filters: by conversion target, by software, by size."""
    df = catalog
    if target:
        df = df[df["targets"].apply(lambda ts: target in ts)]
    if software and software != "All":
        df = df[df["software_name"] == software]
    if max_size_mb is not None:
        df = df[df["size_mb"] <= max_size_mb]
    return df.reset_index(drop=True)


def is_heavy(target: str, size_mb: float) -> bool:
    """Whether converting ``target`` on a ``size_mb`` input is likely memory-heavy."""
    return target in _HEAVY_TARGETS and size_mb >= _HEAVY_SIZE_MB


def _empty_converted_runs() -> pd.DataFrame:
    return pd.DataFrame(columns=CONVERTED_COLUMNS)


def _parse_run_dir_name(name: str) -> dict[str, str]:
    match = _RUN_DIR_RE.match(name)
    if match is None:
        return {"timestamp": "", "slug": "", "target": ""}
    return match.groupdict()


def _result_file(run_dir: Path) -> Path | None:
    for name in ("result.h5mu", "result.h5ad"):
        path = run_dir / name
        if path.is_file():
            return path
    return None


def _command_metadata(log_path: Path) -> dict[str, str]:
    if not log_path.exists():
        return {}
    first_line = log_path.read_text(encoding="utf-8", errors="replace").splitlines()[0:1]
    if not first_line or not first_line[0].startswith("$ "):
        return {}
    try:
        parts = shlex.split(first_line[0][2:])
    except ValueError:
        return {}
    out: dict[str, str] = {}
    for flag, key in {
        "--input": "input_file_path",
        "--slug": "slug",
        "--target": "target",
        "--params": "param_path",
        "--outdir": "output_dir",
    }.items():
        if flag in parts:
            index = parts.index(flag)
            if index + 1 < len(parts):
                out[key] = parts[index + 1]
    return out


def _catalog_lookup(catalog: pd.DataFrame) -> dict[str, dict[str, Any]]:
    if catalog.empty:
        return {}
    return {
        str(row["input_file_path"]): {
            "software_name": row.get("software_name", ""),
            "software_version": row.get("software_version", ""),
            "nr_prec": row.get("nr_prec", ""),
            "size_mb": row.get("size_mb", ""),
            "param_path": row.get("param_path", ""),
        }
        for _, row in catalog.iterrows()
    }


def _stored_search_parameters(obj) -> dict[str, Any] | None:
    params = read_search_parameters(obj)
    if params is None:
        return None
    out = params.model_dump(mode="json", exclude_none=True)
    path = get_search_parameters_path(obj)
    if path:
        out["search_parameters_path"] = path
    return out


def _artifact_search_parameters(result_path: Path) -> dict[str, Any] | None:
    try:
        obj = load_converted_result(result_path)
    except Exception:  # noqa: BLE001 - table metadata falls back to log/catalog.
        return None
    if hasattr(obj, "mod"):
        params = _stored_search_parameters(obj)
        if params is not None:
            return params
        for modality in obj.mod.values():
            params = _stored_search_parameters(modality)
            if params is not None:
                return params
        return None
    return _stored_search_parameters(obj)


def list_converted_runs(converted_dir: Path | str = CONVERTED_DIR) -> pd.DataFrame:
    """Scan converted-output folders into a stable DataFrame for the GUI table."""
    root = Path(converted_dir).expanduser()
    if not root.exists():
        return _empty_converted_runs()

    catalog_by_input = _catalog_lookup(load_catalog())
    rows: list[dict[str, str]] = []
    for run_dir in sorted((p for p in root.iterdir() if p.is_dir()), reverse=True):
        parsed = _parse_run_dir_name(run_dir.name)
        result = _result_file(run_dir)
        log_path = run_dir / "console.log"
        command = _command_metadata(log_path)
        input_file_path = command.get("input_file_path", "")
        catalog_row = catalog_by_input.get(input_file_path, {})
        artifact_params = _artifact_search_parameters(result) if result is not None else None
        if result is not None:
            status = "finished"
            result_type = result.suffix.removeprefix(".")
            result_path = str(result)
        elif log_path.exists():
            status = "incomplete"
            result_type = ""
            result_path = ""
        else:
            status = "empty"
            result_type = ""
            result_path = ""
        rows.append(
            {
                "run_name": run_dir.name,
                "timestamp": parsed["timestamp"],
                "software_name": str(
                    (artifact_params or {}).get(
                        "software_name", catalog_row.get("software_name", "")
                    )
                ),
                "software_version": str(
                    (artifact_params or {}).get(
                        "software_version", catalog_row.get("software_version", "")
                    )
                ),
                "slug": command.get("slug", parsed["slug"]),
                "target": command.get("target", parsed["target"]),
                "status": status,
                "result_type": result_type,
                "nr_prec": str(catalog_row.get("nr_prec", "")),
                "size_mb": str(catalog_row.get("size_mb", "")),
                "input_file_path": input_file_path,
                "param_path": str(
                    (artifact_params or {}).get(
                        "search_parameters_path",
                        command.get("param_path", str(catalog_row.get("param_path", ""))),
                    )
                ),
                "output_dir": str(run_dir),
                "result_path": result_path,
                "log_path": str(log_path) if log_path.exists() else "",
            }
        )
    if not rows:
        return _empty_converted_runs()
    return pd.DataFrame(rows, columns=CONVERTED_COLUMNS)


def converted_runs_table(runs: pd.DataFrame) -> pd.DataFrame:
    """User-facing converted-runs table with no internal filesystem-path columns."""
    columns = [
        "run_name",
        "software_name",
        "software_version",
        "target",
        "status",
        "result_type",
        "nr_prec",
        "size_mb",
    ]
    if runs.empty:
        return pd.DataFrame(columns=columns)
    return runs[columns].copy()


def load_converted_result(result_path: Path | str) -> Any:
    """Load a converted ``result.h5ad`` or ``result.h5mu`` file."""
    path = Path(result_path).expanduser()
    if path.suffix == ".h5ad":
        import anndata as ad

        return ad.read_h5ad(path)
    if path.suffix == ".h5mu":
        import mudata

        # Adopt the mudata 0.4 default now (no auto-pull of per-modality obs/var into the
        # global frames); modalities keep their own obs/var. Silences the 0.3 FutureWarning.
        with mudata.set_options(pull_on_update=False):
            return mudata.read_h5mu(path)
    raise ValueError(f"unsupported converted result type: {path}")


def _noop(_msg: str) -> None:
    pass


def convert_target(
    input_file_path: str,
    slug: str,
    target: str,
    *,
    param_path: Path | str,
    log: Callable[[str], None] = _noop,
):
    """Convert one cached dataset to an AnnData (a level) or a MuData (``"mudata"``).

    The param file is mandatory: it provides the software version that selects the rule variant.
    """
    from anndata_proteomics.readers.dispatch import read_table

    version = _param_version(Path(param_path), slug)
    log(f"param: {Path(param_path).name} → software_version={version!r}")
    df = read_table(_dataset_path(input_file_path))
    log(f"  rows={len(df)} cols={len(df.columns)}")
    if target == MUDATA:
        return _build_mudata(df, slug, version, params_path=param_path, log=log)
    return _convert_level(df, slug, target, version, params_path=param_path, log=log)


def _convert_level(
    df: pd.DataFrame,
    slug: str,
    level: str,
    version: str | None,
    *,
    params_path: Path | str | None = None,
    log: Callable[[str], None] = _noop,
):
    from anndata_proteomics.converters.assemble import convert

    rule = select_rule(slug, level, version, df.columns)
    adata = convert(df, rule, params_path=params_path)
    log(f"  {level}: {adata.shape[0]} obs × {adata.shape[1]} var")
    return adata


def _build_mudata(
    df: pd.DataFrame,
    slug: str,
    version: str | None,
    *,
    params_path: Path | str | None = None,
    log: Callable[[str], None] = _noop,
):
    """Build a MuData over the levels whose version-selected rule fits this file (shared run axis).

    Levels the version doesn't provide (e.g. fragment on DIA-NN 2.x) are skipped, not failed.
    """
    import mudata
    from mudata import MuData

    resolvable = set(convertible_levels(slug, version, df.columns))
    skipped = [level for level in LEVELS if level not in resolvable]
    if skipped:
        log(f"skipping levels not provided by software version {version!r}: {skipped}")
    if len(resolvable) < 2:
        raise ValueError(
            f"{slug}: fewer than two levels resolve for software version {version!r} "
            f"(resolved: {sorted(resolvable)}) — nothing to wrap in a MuData"
        )

    mods = {}
    for level in LEVELS:
        if level not in resolvable:
            continue
        log(f"converting level: {level}")
        adata = _convert_level(df.copy(), slug, level, version, params_path=params_path, log=log)
        adata.var_names = [_PREFIX[level] + str(v) for v in adata.var_names]
        mods[level] = adata
    # Adopt the mudata 0.4 default now (no auto-pull of per-modality obs/var into the global
    # frames); each modality keeps its own obs/var. Silences the 0.3 FutureWarning.
    with mudata.set_options(pull_on_update=False):
        md = MuData(mods, axis=0)
    if params_path is not None:
        params = parse_params(params_path, software=slug)
        write_search_parameters(md, params, source_path=str(params_path))
    log(f"  MuData: {md.n_obs} obs × {sum(a.n_vars for a in mods.values())} var, {len(mods)} mods")
    return md


def _matrix_stats(matrix: np.ndarray) -> dict[str, float]:
    flat = np.asarray(matrix, dtype="float64").ravel()
    valid = flat[~np.isnan(flat)]
    if valid.size == 0:
        return {"min": float("nan"), "max": float("nan"), "mean": float("nan"), "nan_pct": 100.0}
    return {
        "min": float(np.min(valid)),
        "max": float(np.max(valid)),
        "mean": float(np.mean(valid)),
        "nan_pct": round(100 * (flat.size - valid.size) / flat.size, 1),
    }


def _search_parameters_summary(obj) -> dict[str, Any] | None:
    params = _stored_search_parameters(obj)
    if params is None:
        return None
    headline = [
        "software_name",
        "software_version",
        "search_engine_version",
        "quantification_method",
        "ident_fdr_psm",
        "ident_fdr_peptide",
        "ident_fdr_protein",
        "enable_match_between_runs",
    ]
    ordered = {key: params[key] for key in headline if key in params}
    ordered.update({key: value for key, value in params.items() if key not in ordered})
    return ordered


def _mudata_search_parameters_summary(obj) -> dict[str, Any] | None:
    params = _search_parameters_summary(obj)
    if params is not None:
        return params
    summaries = {name: _search_parameters_summary(modality) for name, modality in obj.mod.items()}
    present = {name: summary for name, summary in summaries.items() if summary is not None}
    if not present:
        return None
    first = next(iter(present.values()))
    mismatches = [name for name, summary in present.items() if summary != first]
    out = dict(first)
    out["source"] = "modalities"
    out["modalities"] = list(present)
    if mismatches:
        out["mismatched_modalities"] = mismatches
    return out


def _summarize_anndata(obj, *, include_search_parameters: bool = True) -> dict:
    summary = {
        "kind": "AnnData",
        "shape": (int(obj.n_obs), int(obj.n_vars)),
        "obs_columns": list(obj.obs.columns),
        "var_columns": list(obj.var.columns),
        "layers": list(obj.layers.keys()),
        "uns_keys": list(obj.uns.keys()),
        "x_stats": _matrix_stats(obj.X),
    }
    if include_search_parameters:
        summary["search_parameters"] = _search_parameters_summary(obj)
    return summary


def summarize(obj) -> dict:
    """Summary dict for an AnnData, or per-modality for a MuData. GUI-renderable."""
    if hasattr(obj, "mod"):  # MuData
        search_parameters = _mudata_search_parameters_summary(obj)
        return {
            "kind": "MuData",
            "n_obs": int(obj.n_obs),
            "uns_keys": list(obj.uns.keys()),
            "search_parameters": search_parameters,
            "modalities": {
                name: _summarize_anndata(ad, include_search_parameters=search_parameters is None)
                for name, ad in obj.mod.items()
            },
        }
    return _summarize_anndata(obj)
