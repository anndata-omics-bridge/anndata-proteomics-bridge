# PLAN 2026-05-03 — `converters/recognize.py` (RESTART_PLAN step 6)

## Context

Step 6: given the column headers of a vendor file, identify which packaged `ParseRule` matches. This is the bridge between "I have a DataFrame" and "I know how to convert it" — it lets the future `convert` CLI work without the user having to specify the vendor explicitly.

The signal differs by shape:
- **Long rules**: every vendor column referenced in `columns.obs` (RHS), `columns.var` (RHS), and `layers[*].source_column` must be present in the headers.
- **Wide rules**: every layer's `column_pattern` regex must match at least one header, plus `columns.var` (RHS) values must be present.

If exactly one packaged rule matches → return it. If zero or multiple match → return None (ambiguous; caller must specify).

## Files

```
src/anndata_proteomics/converters/__init__.py   empty
src/anndata_proteomics/converters/recognize.py  matches(headers, rule) + recognize(headers) -> ParseRule | None
tests/test_recognize.py                         per-vendor recognition + ambiguity + nothing-matches
```

## API

```python
def matches(headers: list[str] | set[str], rule: ParseRule) -> bool:
    """Does the given header set plausibly match this rule?"""

def recognize(headers: list[str] | set[str]) -> ParseRule | None:
    """Find the unique packaged ParseRule that matches; None if zero or multiple match."""
```

## Tests

For each packaged TOML, read the matching test_data_download file's header (via the same lookup we use in test_readers_integration.py), call `recognize(headers)`, assert it returns the correct rule. Plus negative cases:
- Empty headers → None
- Random unrelated headers → None
- Headers that match only some columns → None (long: exact subset required)

DRY: parametrize over `iter_packaged_rules()`, look up data file via `raw_file_db_downloaded.csv`. Skip when cache absent.

## Out of scope

- Heuristics for picking among multiple matches (use explicit `--rule` flag instead).
- Performance optimisation — 6 rules × 100 columns is microseconds.

## Verification

```bash
.venv/bin/python -m pytest tests/test_recognize.py -v
```
