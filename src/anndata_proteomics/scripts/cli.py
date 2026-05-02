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

from anndata_proteomics.rules import _export_schema
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
def convert(data: Path, rule_toml: Path) -> int:
    """STUB. Not yet implemented; see docs/RESTART_PLAN.md steps 5-10."""
    print(
        "error: convert is not yet implemented; "
        "see docs/RESTART_PLAN.md steps 5-10",
        file=sys.stderr,
    )
    return 2


def main() -> int:
    """Console-script entry point."""
    rc = app()
    return int(rc) if rc is not None else 0


if __name__ == "__main__":
    sys.exit(main())
