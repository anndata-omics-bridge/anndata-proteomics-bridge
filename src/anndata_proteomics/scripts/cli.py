"""anndata-proteomics CLI dispatcher.

Subcommands:
- validate [path ...]     validate one or more TOML rules; defaults to all packaged
- list                    list packaged rules
- export-schema           regenerate parse_rule.schema.json
- convert <data> <toml>   STUB until readers/ + converters/ land
"""

from __future__ import annotations

import sys
from pathlib import Path

from cyclopts import App

from anndata_proteomics.converters.assemble import convert as _run_convert
from anndata_proteomics.converters.recognize import recognize
from anndata_proteomics.readers.dispatch import read_table
from anndata_proteomics.rules import _export_schema
from anndata_proteomics.rules.loader import load_rule
from anndata_proteomics.rules.registry import iter_packaged_rules
from anndata_proteomics.rules.validate import (
    _print_and_exit_code,
    validate_all_packaged,
    validate_file,
)

app = App(name="anndata-proteomics", help="anndata_proteomics CLI")


@app.command
def validate(*paths: Path) -> int:
    """Validate one or more TOML rule files.

    With no paths, walks all packaged rules (same as `validate-rules`).
    """
    if not paths:
        results = validate_all_packaged()
    else:
        results = [validate_file(p) for p in paths]
    return _print_and_exit_code(results)


@app.command(name="list")
def list_rules() -> int:
    """List packaged parsing rules: software, level, file_version, path."""
    for p in iter_packaged_rules():
        parts = p.stem.split("_")
        # filename: parse_<software_tokens>_<level>_<version>.toml
        software = "_".join(parts[1:-2])
        level = parts[-2]
        version = parts[-1]
        print(f"{software:14}  {level:12}  v{version:<3}  {p}")
    return 0


@app.command(name="export-schema")
def export_schema_cmd() -> int:
    """Regenerate parse_rule.schema.json from the pydantic models."""
    _export_schema.main()
    return 0


@app.command
def convert(
    data: Path,
    rule_toml: Path | None = None,
    output: Path | None = None,
) -> int:
    """Convert a vendor file to AnnData and write a .h5ad.

    If --rule-toml is omitted, the rule is auto-recognized from the data's
    column headers. Use --rule-toml to override or when recognition is ambiguous.
    --output defaults to <data>.h5ad next to the input.
    """
    df = read_table(data)
    if rule_toml is None:
        rule = recognize(list(df.columns))
        if rule is None:
            print(
                f"error: could not auto-recognize a rule for {data}; "
                f"pass --rule-toml PATH",
                file=sys.stderr,
            )
            return 1
    else:
        rule = load_rule(rule_toml)
    adata = _run_convert(df, rule)
    out = output or data.with_suffix(".h5ad")
    adata.write_h5ad(out)
    print(f"wrote {out}  shape={adata.shape}  layers={list(adata.layers)}")
    return 0


def main() -> int:
    """Console-script entry point."""
    rc = app()
    return int(rc) if rc is not None else 0


if __name__ == "__main__":
    sys.exit(main())
