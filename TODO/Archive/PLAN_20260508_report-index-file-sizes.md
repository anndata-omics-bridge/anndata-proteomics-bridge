# Plan: Add file sizes to report index

## Goal

Show both source input file size and generated AnnData `.h5ad` file size in the
`tools/generate_report.py` `index.html` table.

## Current behavior

`tools/generate_report.py` writes one `<stem>.meta.json` per packaged rule and then
`rebuild_index(output_dir)` builds `<output-dir>/index.html` from those metadata files.
The table currently has:

- software
- input
- output (`.h5ad`)
- dim / layers
- report
- log

The metadata does not store file sizes, and the index does not display them.

## Proposed implementation

1. Add two optional fields to the `Outcome` dataclass:
   - `input_size_bytes`
   - `h5ad_size_bytes`
2. Populate `input_size_bytes` from `input_path.stat().st_size` when an input file exists.
3. Populate `h5ad_size_bytes` after `adata.write_h5ad(h5ad_path)` when the `.h5ad` exists.
4. Write both values into `<stem>.meta.json`.
5. Add a small formatter such as `_format_bytes(size)` for human-readable values:
   - bytes below 1024 as `B`
   - then `KiB`, `MiB`, `GiB`
6. Add two index columns:
   - `input size`
   - `.h5ad size`
7. For skipped/failed rows where a file is unavailable, show `(none)` or an empty value
   consistently.
8. Update `tests/test_generate_report.py` to assert that:
   - generated metadata contains file-size keys
   - successful rows include size text in `index.html`
   - skipped rows keep `None` sizes
9. Update README/report docs if needed to mention that the index includes input and output sizes.

## Scope

Only `tools/generate_report.py`, `tests/test_generate_report.py`, and possibly README/docs
should change. No parsing-rule schema changes and no report-rendering changes in `annProtSum`.

## Verification

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_generate_report.py
```

If that passes quickly, optionally run the full suite:

```bash
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider
```
