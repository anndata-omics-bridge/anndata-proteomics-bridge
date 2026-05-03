# PLAN 2026-05-03 — `converters/long.py / wide.py / factors.py / assemble.py` (RESTART_PLAN steps 7–10)

## Context

Steps 7–10: actually convert a `(DataFrame, ParseRule)` into AnnData pieces. After this lands, the package can read a vendor file end-to-end and produce an AnnData object — completing the core "vendor file + parsing TOML → AnnData" goal of the restart.

Splitting the work across four files matches the RESTART_PLAN target and keeps each concern small:
- `long.py` — pivot a long DataFrame into per-layer (obs × var) matrices.
- `wide.py` — extract sample tokens from column headers via regex, build per-layer matrices.
- `factors.py` — encode string-valued layers to integer codes per `categories`.
- `assemble.py` — turn the pieces into an `AnnData` object, plus minimal `uns` provenance.

## Files

```
src/anndata_proteomics/converters/long.py       convert_long(df, rule) -> ConversionPieces
src/anndata_proteomics/converters/wide.py       convert_wide(df, rule) -> ConversionPieces
src/anndata_proteomics/converters/factors.py    encode_factor(series, categories) -> Series[int]
src/anndata_proteomics/converters/assemble.py   to_anndata(pieces, rule) -> AnnData
tests/test_converters_long.py                   long conversion happy path + dup handling
tests/test_converters_wide.py                   wide conversion happy path + factor layers
tests/test_converters_factors.py                encode_factor unit tests
tests/test_converters_assemble.py               assemble pieces → AnnData smoke
tests/test_converters_e2e.py                    one parametrized end-to-end test per packaged TOML
```

## Shared types

```python
@dataclass
class ConversionPieces:
    X: np.ndarray
    obs: pd.DataFrame      # indexed by obs_names (sample identifiers)
    var: pd.DataFrame      # indexed by var_names (joined feature keys)
    layers: dict[str, np.ndarray]
    uns: dict
```

Stored in `converters/__init__.py` for import by long.py, wide.py, assemble.py.

## Long converter

Input: DataFrame where each row is `(obs_key_values, var_key_values, layer_value_for_each_layer.source_column)`.

Algorithm:
1. Build the obs index — `df[obs_keys].drop_duplicates()` joined as string. obs DataFrame: pick the first row per obs-key tuple, take `columns.obs` columns and rename vendor → internal.
2. Build the var index — same but for `var_keys` and `columns.var`.
3. For each layer: `df.pivot_table(index=obs_keys, columns=var_keys, values=source_column, aggfunc=...)`. Reindex rows / columns to match obs / var order. Resulting `(n_obs × n_var)` matrix.
4. `X = layers[axis.x_layer]`.

Duplicate policy (`duplicates.mode`):
- `"error"` → use `aggfunc='first'` and validate uniqueness afterwards (raise if dup with different values).
- `"keep_first"` → `aggfunc='first'`, no validation.
- `"aggregate"` → `aggfunc='sum'` with `numeric_only=True`.
- `"keep_all_as_raw_table"` → not implemented in this pass; raise NotImplementedError.

For factor-encoded layers: pre-map `df[source_column]` through `encode_factor(...)` before pivoting so the resulting matrix is numeric.

## Wide converter

Input: DataFrame where columns include both var metadata and per-sample value columns matching `layer.column_pattern` regexes.

Algorithm:
1. For each layer, compile the regex and find matching columns. Extract `<sample>` token from each column via the named capture group.
2. Sample tokens (across all layers) form the `obs` axis. Sanity check: every layer must share the same set of sample tokens (or a subset; warn / error on mismatch).
3. Var axis: `df[var_keys].drop_duplicates()` joined as string. var DataFrame from `columns.var` (skipping the `"<sample>"` placeholder).
4. For each layer: gather the matching columns into an `(n_var × n_samples)` matrix from the DataFrame, transpose to `(n_obs × n_var)`, reorder rows / columns. For factor-encoded layers, apply `encode_factor` to each column's values before stacking.
5. obs DataFrame: just sample tokens for now (no per-sample metadata in wide rules — the data has no place for it; future plan can add cleanup-rule-derived obs columns).
6. Apply `sample_name_cleanup.pattern` to obs_names if present.

## Factors

```python
def encode_factor(series: pd.Series, categories: dict[str, int], default: int = -1) -> pd.Series:
    """Map string values to integer codes; unknown / NaN values become `default`."""
    out = series.map(categories)
    out = out.fillna(default).astype("int64")
    return out
```

## Assemble

```python
def to_anndata(pieces: ConversionPieces, rule: ParseRule) -> ad.AnnData:
    adata = ad.AnnData(X=pieces.X, obs=pieces.obs, var=pieces.var, layers=pieces.layers)
    adata.uns["anndata_proteomics"] = {
        "rule": rule.model_dump(mode="json"),
        "schema_version": rule.schema_version,
    }
    return adata
```

Per the ADR, future per-tool `column_roles` namespaces (e.g. `uns['exploreDE']`, `uns['bottom_up_proteomics']`) are out of scope for this plan — they get added when those consumers are wired up. For now we record the source rule under `uns['anndata_proteomics']` so downstream code can introspect.

## Tests

### Per-module unit tests

- `test_converters_long.py` — synthetic 3-sample × 2-feature long DataFrame; assert obs/var/X/layer shapes and values; dup handling.
- `test_converters_wide.py` — synthetic with `^(?P<sample>S\d+) Intensity$` columns; assert sample extraction and matrix shape; one factor layer.
- `test_converters_factors.py` — encode_factor with known categories, NaN handling, unknown-value default.
- `test_converters_assemble.py` — feed minimal pieces, assert AnnData has correct shape + `uns['anndata_proteomics']['rule']` round-trips.

### End-to-end DRY test

`test_converters_e2e.py` — one parametrized test that:

1. Iterates every packaged TOML.
2. Looks up matching test_data_download file via `raw_file_db_downloaded.csv` (skip if absent).
3. `read_table → recognize → convert_long/wide → to_anndata`.
4. Asserts non-empty `n_obs > 0`, `n_var > 0`, X has expected shape.

Single piece of test code covers all 6 vendors. Add a 7th TOML and the test extends automatically.

## Dependencies

- `anndata` already in deps.
- `pandas`, `numpy` already in deps.

## Verification

```bash
.venv/bin/python -m pytest tests/test_converters_*.py -v
.venv/bin/python -m pytest tests/ -q

# Smoke from Python
.venv/bin/python - <<'PY'
from anndata_proteomics.readers.dispatch import read_table
from anndata_proteomics.converters.recognize import recognize
from anndata_proteomics.converters.long import convert_long
from anndata_proteomics.converters.wide import convert_wide
from anndata_proteomics.converters.assemble import to_anndata
from pathlib import Path

# Pick e.g. WOMBAT (smallest)
import csv
with open("test_data_download/raw_file_db_downloaded.csv") as f:
    for row in csv.DictReader(f):
        if row["software_name"] == "WOMBAT" and row.get("status") == "ok":
            data_file = Path("test_data_download/json_dir") / row["input_file_path"]
            break
df = read_table(data_file)
rule = recognize(list(df.columns))
print("recognised:", rule.software_name, rule.input_shape)
pieces = convert_wide(df, rule) if rule.input_shape == "wide" else convert_long(df, rule)
adata = to_anndata(pieces, rule)
print("AnnData:", adata)
PY
```

## Out of scope

- `keep_all_as_raw_table` duplicate mode (raise NotImplementedError; no current TOML uses it).
- Per-tool `uns['<app_name>']['column_roles']` writeback (only `uns['anndata_proteomics']`).
- `obs` enrichment from `sample_name_cleanup.pattern` extracted groups (we apply the rename but don't promote groups to obs columns).
- Performance — pivot_table on multi-million-row DIA-NN reports is expected to take seconds; if it becomes a blocker, switch to manual unstack.
- Anything that hits `proteome_discoverer` (no rule, no test data).
