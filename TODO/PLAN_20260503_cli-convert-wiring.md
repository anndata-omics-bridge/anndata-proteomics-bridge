# PLAN 2026-05-03 — Wire `anndata-proteomics convert` to the real conversion pipeline

## Context

The conversion pipeline (`read_table → recognize → convert(df, rule) → AnnData`) is fully implemented and tested. The CLI subcommand `anndata-proteomics convert` was left as a stub (returns exit 2, "not yet implemented") in [src/anndata_proteomics/scripts/cli.py](../src/anndata_proteomics/scripts/cli.py). This plan replaces the stub with a real implementation so end users can convert vendor files from the command line.

## Subcommand surface

```
anndata-proteomics convert <data> [--rule-toml PATH] [--output PATH]
```

- `<data>` — input vendor file (.csv / .tsv / .txt / .parquet).
- `--rule-toml PATH` — optional. If omitted, auto-recognize from headers via `converters.recognize.recognize(...)`. Fails (exit 1) if zero or multiple packaged rules match.
- `--output PATH` — optional. Defaults to `<data>.h5ad` next to the input.

Exit codes:
- `0` — wrote the .h5ad successfully.
- `1` — recognition failed, file unreadable, or conversion error.

## Files

- **`src/anndata_proteomics/scripts/cli.py`** — replace the `convert` stub.
- **`tests/test_cli.py`** — replace `test_convert_is_stub` with positive + recognition-failure cases.
- **`tests/test_cli_integration.py`** — replace `test_cli_convert_stub_returns_two` with a happy-path subprocess test.
- **`README.md`** — bump the Quickstart to mention `anndata-proteomics convert <data>`.

## Implementation

```python
from anndata_proteomics.converters.assemble import convert as _run_convert
from anndata_proteomics.converters.recognize import recognize
from anndata_proteomics.readers.dispatch import read_table
from anndata_proteomics.rules.loader import load_rule

@app.command
def convert(
    data: Path,
    rule_toml: Path | None = None,
    output: Path | None = None,
) -> int:
    """Convert a vendor file to AnnData and write a .h5ad.

    If --rule-toml is omitted, the rule is auto-recognized from the data's
    column headers. Use --rule-toml when recognition is ambiguous or to
    apply a custom rule.
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
```

## Tests

- Unit (test_cli.py):
  - convert with explicit rule_toml → writes h5ad, returns 0.
  - convert with auto-recognition → writes h5ad.
  - convert with unrecognizable headers (synthetic file) → returns 1.
- Integration (test_cli_integration.py):
  - One subprocess test: pick a small packaged vendor's data file, invoke `anndata-proteomics convert <data> --output <tmp>`, assert exit 0 and the file exists.

## Verification

```bash
.venv/bin/python -m pytest tests/ -q
anndata-proteomics convert test_data_download/json_dir/Results_quant_peptidoform_DDA/.../input_file.csv \
  --output /tmp/wombat.h5ad
ls -la /tmp/wombat.h5ad
.venv/bin/python -c "import anndata; a = anndata.read_h5ad('/tmp/wombat.h5ad'); print(a)"
```

## Out of scope

- A `--validate-only` flag (use `anndata-proteomics validate` for that).
- Multi-file batch conversion. Single file in / single file out for now.
- Output formats other than .h5ad. AnnData has zarr support too; add later if needed.
- Progress bars / verbose mode.
