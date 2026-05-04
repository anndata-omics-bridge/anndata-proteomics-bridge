"""Convert one vendor file to AnnData and render an annProtSum HTML report.

Per call: read the vendor file → recognize / use given rule → convert to AnnData
→ write `<software>_<sha8(input_path)>.h5ad` and a sidecar `.meta.json` under
`--output-dir` (default `examples/results/`) → shell out to the annProtSum R
package's `render_report.R` to render an HTML report → rebuild
`<output-dir>/index.html` with one row per conversion (input, .h5ad,
per-layer shape, report).

Usage:
    python tools/generate_report.py path/to/data.tsv
    python tools/generate_report.py path/to/data.tsv --rule-toml my.toml
    python tools/generate_report.py path/to/data.tsv --output-dir examples/results
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from anndata_proteomics.converters.assemble import convert as run_convert
from anndata_proteomics.converters.recognize import recognize
from anndata_proteomics.readers.dispatch import read_table
from anndata_proteomics.rules.loader import load_rule
from anndata_proteomics.rules.schema import ParseRule


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT_DIR = REPO_ROOT / "examples" / "results"
WORKSPACE_ROOT = REPO_ROOT.parent
DEV_RENDER_SCRIPT = (
    WORKSPACE_ROOT / "annProtSum" / "inst" / "bin" / "render_report.R"
)


@dataclass(frozen=True)
class Conversion:
    input_path: Path
    h5ad_path: Path
    html_path: Path
    meta_path: Path
    software: str
    n_obs: int
    n_var: int
    layers: list[str]


def _sha8(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()[:8]


def _stem_for(input_path: Path, rule: ParseRule) -> str:
    software = rule.software_name.lower().replace(" ", "_").replace("-", "")
    return f"{software}_{_sha8(str(input_path.resolve()))}"


def _resolve_render_script() -> Path:
    """Find render_report.R via installed annProtSum, fall back to dev path."""
    try:
        out = subprocess.run(
            ["R", "-q", "-s", "-e",
             'cat(system.file("bin/render_report.R", package = "annProtSum"))'],
            capture_output=True, text=True, check=True, timeout=30,
        ).stdout.strip()
        if out and Path(out).exists():
            return Path(out)
    except (FileNotFoundError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
        pass
    if DEV_RENDER_SCRIPT.exists():
        return DEV_RENDER_SCRIPT
    raise RuntimeError(
        "could not locate render_report.R. Install annProtSum "
        "(devtools::install('/Users/wolski/projects/anndata_bridge/annProtSum', "
        "dependencies = TRUE)) or ensure the dev path exists at "
        f"{DEV_RENDER_SCRIPT}"
    )


def convert_one(
    input_path: Path,
    rule_toml: Path | None,
    output_dir: Path,
) -> Conversion:
    """Run the read → recognize/load → convert → render pipeline for one file."""
    output_dir.mkdir(parents=True, exist_ok=True)
    df = read_table(input_path)
    if rule_toml is not None:
        rule = load_rule(rule_toml)
    else:
        rule = recognize(list(df.columns))
        if rule is None:
            raise SystemExit(
                f"could not auto-recognize a rule for {input_path}; "
                f"pass --rule-toml PATH"
            )

    adata = run_convert(df, rule)
    stem = _stem_for(input_path, rule)
    h5ad_path = output_dir / f"{stem}.h5ad"
    html_path = output_dir / f"{stem}.html"
    meta_path = output_dir / f"{stem}.meta.json"

    adata.write_h5ad(h5ad_path)

    layers = list(adata.layers)
    meta = {
        "input_path": str(input_path),
        "h5ad_path": h5ad_path.name,
        "html_path": html_path.name,
        "software": rule.software_name,
        "input_shape": rule.input_shape,
        "quantification_level": rule.quantification_level,
        "n_obs": int(adata.n_obs),
        "n_var": int(adata.n_vars),
        "layers": [
            {"name": name, "n_obs": int(adata.layers[name].shape[0]),
             "n_var": int(adata.layers[name].shape[1])}
            for name in layers
        ],
    }
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")

    render_script = _resolve_render_script()
    subprocess.run(
        ["Rscript", str(render_script), str(h5ad_path), str(html_path)],
        check=True,
    )

    return Conversion(
        input_path=input_path,
        h5ad_path=h5ad_path,
        html_path=html_path,
        meta_path=meta_path,
        software=rule.software_name,
        n_obs=int(adata.n_obs),
        n_var=int(adata.n_vars),
        layers=layers,
    )


def rebuild_index(output_dir: Path) -> Path:
    """Globs `<output-dir>/*.meta.json` and writes `<output-dir>/index.html`."""
    metas = sorted(output_dir.glob("*.meta.json"))
    rows = []
    for p in metas:
        m = json.loads(p.read_text())
        layer_lis = "".join(
            f"<li><code>{html.escape(l['name'])}</code>: "
            f"({l['n_obs']}, {l['n_var']})</li>"
            for l in m["layers"]
        )
        rows.append(
            f"<tr>"
            f"<td>{html.escape(m['software'])}</td>"
            f"<td><code>{html.escape(m['input_path'])}</code></td>"
            f"<td><a href=\"{html.escape(m['h5ad_path'])}\">"
            f"{html.escape(m['h5ad_path'])}</a></td>"
            f"<td><ul style=\"margin:0;padding-left:1em\">{layer_lis}</ul></td>"
            f"<td><a href=\"{html.escape(m['html_path'])}\">report</a></td>"
            f"</tr>"
        )
    body = (
        "<table border=\"1\" cellpadding=\"6\" cellspacing=\"0\" "
        "style=\"border-collapse:collapse\">"
        "<thead><tr>"
        "<th>software</th><th>input</th><th>output (.h5ad)</th>"
        "<th>layers (shape)</th><th>report</th>"
        "</tr></thead>"
        f"<tbody>{''.join(rows) if rows else '<tr><td colspan=5>(no conversions yet)</td></tr>'}</tbody>"
        "</table>"
    )
    out = output_dir / "index.html"
    out.write_text(
        f"<!doctype html><html><head><meta charset=\"utf-8\">"
        f"<title>anndata_proteomics conversions</title></head>"
        f"<body><h1>anndata_proteomics conversions</h1>{body}</body></html>\n"
    )
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path, help="vendor data file (.csv / .tsv / .txt / .parquet)")
    parser.add_argument("--rule-toml", type=Path, default=None,
                        help="explicit ParseRule TOML; default = auto-recognize")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR,
                        help=f"where to write artifacts (default: {DEFAULT_OUTPUT_DIR})")
    args = parser.parse_args(argv)

    conv = convert_one(args.data, args.rule_toml, args.output_dir)
    index = rebuild_index(args.output_dir)
    print(f"h5ad:   {conv.h5ad_path}")
    print(f"html:   {conv.html_path}")
    print(f"meta:   {conv.meta_path}")
    print(f"index:  {index}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
