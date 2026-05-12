"""Iterate every packaged ParseRule, convert one canonical input per converter, render reports.

Per packaged rule:
  * resolve a canonical input via `anndata_proteomics.test_data.find_test_data(software_name)`
  * if none → emit a `skipped` row (still gets a .log + .meta.json)
  * else read → recognize/load → convert → write `<stem>.h5ad` + `<stem>.meta.json`
  * shell out to annProtSum's `render_report.R` to render `<stem>.html`
  * any failure → `failed` row with traceback in the per-rule log

`<stem>` = `<software_token>_<sha8(input_path or software_name)>`. Per-run
output goes under `--output-dir` (default `examples/results/`); the run also
produces `<output-dir>/index.html` with one row per packaged rule and links
to: input path + size, .h5ad + size, dim+layer summary, .html report, .log file.

Usage:
    python tools/generate_report.py
    python tools/generate_report.py --output-dir examples/results
    python tools/generate_report.py --rule DIA-NN --rule WOMBAT
    python tools/generate_report.py --log-level DEBUG
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import subprocess
import sys
import traceback
from dataclasses import dataclass
from pathlib import Path

from loguru import logger

from anndata_proteomics._logging import configure_default_sink
from anndata_proteomics.converters.assemble import convert as run_convert
from anndata_proteomics.readers.dispatch import read_table
from anndata_proteomics.rules.loader import load_rule
from anndata_proteomics.rules.registry import iter_packaged_rules
from anndata_proteomics.rules.schema import ParseRule
from anndata_proteomics.test_data import find_param_file, find_test_data

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT_DIR = REPO_ROOT / "examples" / "results"
WORKSPACE_ROOT = REPO_ROOT.parent
DEV_RENDER_SCRIPT = WORKSPACE_ROOT / "annProtSum" / "inst" / "bin" / "render_report.R"

LOG_FORMAT_FILE = "{time:YYYY-MM-DD HH:mm:ss} | {level: <7} | {message}"


@dataclass(frozen=True)
class Outcome:
    status: str  # "ok" | "skipped" | "failed"
    software: str
    stem: str
    input_path: Path | None
    h5ad_path: Path | None
    html_path: Path | None
    log_path: Path
    meta_path: Path
    input_size_bytes: int | None
    h5ad_size_bytes: int | None
    n_obs: int | None
    n_var: int | None
    layers: list[dict[str, int | str]]
    error: str | None


def _sha8(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()[:8]


def _software_token(rule: ParseRule) -> str:
    return rule.software_name.lower().replace(" ", "_").replace("-", "")


def _stem_for(rule: ParseRule, input_path: Path | None) -> str:
    seed = str(input_path.resolve()) if input_path is not None else rule.software_name
    return f"{_software_token(rule)}_{_sha8(seed)}"


def _path_size(path: Path | None) -> int | None:
    """Return file size in bytes when `path` points at an existing file."""
    if path is None or not path.exists():
        return None
    return path.stat().st_size


def _format_bytes(size: int | None) -> str:
    """Format byte counts for compact display in the HTML index."""
    if size is None:
        return "(none)"
    value = float(size)
    for unit in ("B", "KiB", "MiB", "GiB"):
        if value < 1024 or unit == "GiB":
            if unit == "B":
                return f"{int(value)} {unit}"
            return f"{value:.1f} {unit}"
        value /= 1024
    return f"{value:.1f} GiB"


def _resolve_render_script() -> Path:
    """Find render_report.R via installed annProtSum, fall back to dev path."""
    try:
        out = subprocess.run(
            [
                "R",
                "-q",
                "-s",
                "-e",
                'cat(system.file("bin/render_report.R", package = "annProtSum"))',
            ],
            capture_output=True,
            text=True,
            check=True,
            timeout=30,
        ).stdout.strip()
        if out and Path(out).exists():
            return Path(out)
    except (FileNotFoundError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
        pass
    if DEV_RENDER_SCRIPT.exists():
        return DEV_RENDER_SCRIPT
    raise RuntimeError(
        "could not locate render_report.R. Install annProtSum "
        f"(devtools::install('{WORKSPACE_ROOT / 'annProtSum'}', dependencies = TRUE)) "
        f"or ensure the dev path exists at {DEV_RENDER_SCRIPT}."
    )


def _clear_output_dir(output_dir: Path) -> None:
    """Remove prior per-converter outputs and the index. Idempotent."""
    if not output_dir.exists():
        return
    for pattern in ("*.meta.json", "*.h5ad", "*.html", "*.log", "index.html"):
        for p in output_dir.glob(pattern):
            p.unlink()


def _write_meta(outcome: Outcome) -> None:
    meta = {
        "status": outcome.status,
        "software": outcome.software,
        "stem": outcome.stem,
        "input_path": str(outcome.input_path) if outcome.input_path else None,
        "h5ad_path": outcome.h5ad_path.name if outcome.h5ad_path else None,
        "html_path": outcome.html_path.name if outcome.html_path else None,
        "log_path": outcome.log_path.name,
        "input_size_bytes": outcome.input_size_bytes,
        "h5ad_size_bytes": outcome.h5ad_size_bytes,
        "n_obs": outcome.n_obs,
        "n_var": outcome.n_var,
        "layers": outcome.layers,
        "error": outcome.error,
    }
    outcome.meta_path.write_text(json.dumps(meta, indent=2) + "\n")


def _run_one(rule: ParseRule, output_dir: Path, log_level: str) -> Outcome:
    """Convert one rule's canonical input + render report. Always returns an Outcome."""
    output_dir.mkdir(parents=True, exist_ok=True)
    input_path = find_test_data(rule.software_name)
    stem = _stem_for(rule, input_path)
    log_path = output_dir / f"{stem}.log"
    meta_path = output_dir / f"{stem}.meta.json"

    handler_id = logger.add(log_path, format=LOG_FORMAT_FILE, level=log_level)
    try:
        logger.info(f"=== {rule.software_name} ({rule.quantification_level}) ===")
        if input_path is None or not input_path.exists():
            reason = (
                f"no row with status=ok for software_name={rule.software_name!r} "
                f"in test_data_download/raw_file_db_downloaded.csv "
                f"(regenerate via test_data_download/Makefile)"
            )
            logger.warning(reason)
            outcome = Outcome(
                status="skipped",
                software=rule.software_name,
                stem=stem,
                input_path=input_path,
                h5ad_path=None,
                html_path=None,
                log_path=log_path,
                meta_path=meta_path,
                input_size_bytes=None,
                h5ad_size_bytes=None,
                n_obs=None,
                n_var=None,
                layers=[],
                error=reason,
            )
            _write_meta(outcome)
            return outcome

        h5ad_path = output_dir / f"{stem}.h5ad"
        html_path = output_dir / f"{stem}.html"
        try:
            logger.info(f"reading {input_path}")
            df = read_table(input_path)
            logger.info(f"input shape: {df.shape}")

            params_path = find_param_file(rule.software_name)
            if params_path is not None:
                logger.info(f"parameter file: {params_path}")
            else:
                logger.info("no parameter file available for this tool")
            logger.info("converting to AnnData")
            adata = run_convert(df, rule, params_path=params_path)
            adata.write_h5ad(h5ad_path)
            logger.info(f"wrote {h5ad_path}  shape={adata.shape}  layers={list(adata.layers)}")

            layer_descs = [
                {
                    "name": name,
                    "n_obs": int(adata.layers[name].shape[0]),
                    "n_var": int(adata.layers[name].shape[1]),
                }
                for name in adata.layers
            ]

            render_script = _resolve_render_script()
            logger.info(f"rendering report via {render_script}")
            with open(log_path, "a") as fh:
                fh.write(f"\n--- Rscript {render_script} ---\n")
                fh.flush()
                subprocess.run(
                    ["Rscript", str(render_script), str(h5ad_path), str(html_path)],
                    check=True,
                    stdout=fh,
                    stderr=fh,
                )
            logger.info(f"wrote {html_path}")

            outcome = Outcome(
                status="ok",
                software=rule.software_name,
                stem=stem,
                input_path=input_path,
                h5ad_path=h5ad_path,
                html_path=html_path,
                log_path=log_path,
                meta_path=meta_path,
                input_size_bytes=_path_size(input_path),
                h5ad_size_bytes=_path_size(h5ad_path),
                n_obs=int(adata.n_obs),
                n_var=int(adata.n_vars),
                layers=layer_descs,
                error=None,
            )
            _write_meta(outcome)
            return outcome
        except Exception as exc:
            tb = traceback.format_exc()
            logger.error(f"{rule.software_name} failed: {exc}")
            logger.debug(tb)
            outcome = Outcome(
                status="failed",
                software=rule.software_name,
                stem=stem,
                input_path=input_path,
                h5ad_path=h5ad_path if h5ad_path.exists() else None,
                html_path=None,
                log_path=log_path,
                meta_path=meta_path,
                input_size_bytes=_path_size(input_path),
                h5ad_size_bytes=_path_size(h5ad_path),
                n_obs=None,
                n_var=None,
                layers=[],
                error=f"{type(exc).__name__}: {exc}",
            )
            _write_meta(outcome)
            return outcome
    finally:
        logger.remove(handler_id)


def _path_link(path_str: str) -> str:
    """Render an anchor whose visible text is the basename and whose title is the full path."""
    escaped = html.escape(path_str)
    base = html.escape(Path(path_str).name)
    return f'<a href="{escaped}" title="{escaped}"><code>{base}</code></a>'


def _meta_to_row(meta: dict) -> str:
    """Render one <tr> for the index."""
    software = html.escape(meta["software"])
    log_cell = f'<a href="{html.escape(meta["log_path"])}">log</a>' if meta.get("log_path") else ""
    input_cell = _path_link(meta["input_path"]) if meta.get("input_path") else "(none)"
    input_size_cell = html.escape(_format_bytes(meta.get("input_size_bytes")))
    h5ad_size_cell = html.escape(_format_bytes(meta.get("h5ad_size_bytes")))

    if meta["status"] == "ok":
        layer_lis = "".join(
            f"<li><code>{html.escape(l['name'])}</code>: ({l['n_obs']}, {l['n_var']})</li>"
            for l in meta["layers"]
        )
        dim_cell = (
            f"({meta['n_obs']}, {meta['n_var']})"
            f'<ul style="margin:0;padding-left:1em">{layer_lis}</ul>'
        )
        output_cell = _path_link(meta["h5ad_path"])
        report_cell = f'<a href="{html.escape(meta["html_path"])}">report</a>'
    else:
        status_label = meta["status"].upper()
        err = html.escape(meta.get("error") or "")
        dim_cell = f"<em>{status_label}</em>: {err}"
        output_cell = "(none)"
        report_cell = "(none)"

    return (
        "<tr>"
        f"<td>{software}</td>"
        f"<td>{input_cell}</td>"
        f"<td>{input_size_cell}</td>"
        f"<td>{output_cell}</td>"
        f"<td>{h5ad_size_cell}</td>"
        f"<td>{dim_cell}</td>"
        f"<td>{report_cell}</td>"
        f"<td>{log_cell}</td>"
        "</tr>"
    )


def rebuild_index(output_dir: Path) -> Path:
    """Glob `<output-dir>/*.meta.json` and write `<output-dir>/index.html`."""
    metas = sorted(output_dir.glob("*.meta.json"))
    rows = [_meta_to_row(json.loads(p.read_text())) for p in metas]
    body = (
        '<table border="1" cellpadding="6" cellspacing="0" '
        'style="border-collapse:collapse">'
        "<thead><tr>"
        "<th>software</th><th>input</th><th>input size</th>"
        "<th>output (.h5ad)</th><th>.h5ad size</th>"
        "<th>dim / layers</th><th>report</th><th>log</th>"
        "</tr></thead>"
        f"<tbody>{''.join(rows) if rows else '<tr><td colspan=8>(no conversions yet)</td></tr>'}</tbody>"
        "</table>"
    )
    out = output_dir / "index.html"
    out.write_text(
        '<!doctype html><html><head><meta charset="utf-8">'
        "<title>anndata_proteomics conversions</title>"
        "<style>"
        "body{font-family:system-ui,sans-serif;margin:1.5em;}"
        "table{border-collapse:collapse;}"
        "td,th{vertical-align:top;word-break:break-word;}"
        "td code{font-size:0.9em;}"
        "</style></head>"
        f"<body><h1>anndata_proteomics conversions</h1>{body}</body></html>\n"
    )
    return out


def _select_rules(rule_filter: list[str] | None) -> list[ParseRule]:
    rules = [load_rule(p) for p in iter_packaged_rules()]
    if not rule_filter:
        return rules
    wanted = {s.lower() for s in rule_filter}
    selected = [r for r in rules if r.software_name.lower() in wanted]
    if not selected:
        available = sorted(r.software_name for r in rules)
        raise SystemExit(
            f"--rule {rule_filter} matched none of the packaged rules. "
            f"Available software_name values: {available}"
        )
    return selected


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Convert one canonical input per packaged ParseRule and render a report index."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"where to write artifacts (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--rule",
        action="append",
        dest="rule_filter",
        help="restrict to specific software_name(s), case-insensitive; repeatable. "
        "When passed, the output dir is NOT cleared.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="loguru level for the stderr sink and per-rule log files (default: INFO).",
    )
    args = parser.parse_args(argv)

    configure_default_sink(level=args.log_level)

    rules = _select_rules(args.rule_filter)
    output_dir = args.output_dir.resolve()

    if not args.rule_filter:
        logger.info(f"clearing prior outputs in {output_dir}")
        _clear_output_dir(output_dir)

    logger.info(f"running {len(rules)} converter(s); output dir: {output_dir}")
    outcomes = [_run_one(rule, output_dir, args.log_level) for rule in rules]
    index = rebuild_index(output_dir)

    summary_counts = {"ok": 0, "skipped": 0, "failed": 0}
    for o in outcomes:
        summary_counts[o.status] += 1
    logger.info(
        f"summary: ok={summary_counts['ok']}  "
        f"skipped={summary_counts['skipped']}  failed={summary_counts['failed']}"
    )
    logger.info(f"index: {index}")
    return 0 if summary_counts["failed"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
