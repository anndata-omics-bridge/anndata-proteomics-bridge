# PLAN 2026-05-02 — `readers/` (RESTART_PLAN step 5)

## Context

Step 5 of [docs/RESTART_PLAN.md](../docs/RESTART_PLAN.md): generic vendor-file reading. The reader's only job is **`Path -> pandas.DataFrame`** — no vendor semantics, no rule application, no column renaming. Vendor knowledge lives in the TOMLs we already shipped; the reader is deliberately stupid.

User constraint (this turn): keep the reader code DRY across the 6 vendors. **The plan exploits the fact that the readers are already DRY by design** — there is exactly one dispatch function and one tabular reader per file extension. No per-vendor reader code anywhere. Vendor differences are absorbed by:

- **Extension dispatch** for the file format (`.csv` / `.tsv` / `.txt` / `.parquet`).
- **The TOML rule** (already shipped) for what the columns mean.

The 6 packaged vendors cover: `.parquet` (DIA-NN), `.tsv` (Spectronaut, FragPipe), `.txt` (MaxQuant), `.csv` (PEAKS, WOMBAT). Four extensions, one reader per extension, four lines of dispatch.

## Files to create

```
src/anndata_proteomics/readers/__init__.py        empty per Coding Rules
src/anndata_proteomics/readers/tabular.py         read_csv / read_tsv / read_parquet
src/anndata_proteomics/readers/dispatch.py        read_table(path) + UnknownFormat
tests/test_readers_tabular.py                     happy + edge cases per format
tests/test_readers_dispatch.py                    extension dispatch + UnknownFormat
tests/test_readers_integration.py                 every packaged rule reads its test_data_download file (skipif when cache absent)
```

No schema changes. No `[reader]` block on the TOML — current vendor files are all auto-readable by extension. If a future vendor file needs delimiter / encoding / skiprows overrides, that's a follow-up plan that adds an optional `[reader]` section to `ParseRule` then.

## Public API

### `readers/tabular.py`

```python
def read_csv(path: Path | str) -> pd.DataFrame:
    """Read a CSV file. UTF-8 (with BOM tolerance), pandas defaults otherwise."""

def read_tsv(path: Path | str) -> pd.DataFrame:
    """Read a tab-delimited file. UTF-8, pandas defaults."""

def read_parquet(path: Path | str) -> pd.DataFrame:
    """Read a parquet file via pyarrow (already a dep)."""
```

Each one is a 2–4 line wrapper around `pd.read_csv` / `pd.read_parquet`. The wrappers exist so the test suite has stable function entry points and so future overrides land in one place.

### `readers/dispatch.py`

```python
class UnknownFormat(ValueError):
    """Raised when a file extension is not recognised."""

EXTENSION_TO_READER = {
    ".csv":     read_csv,
    ".tsv":     read_tsv,
    ".txt":     read_tsv,    # MaxQuant evidence.txt and friends are tab-delimited
    ".parquet": read_parquet,
}

def read_table(path: Path | str) -> pd.DataFrame:
    """Dispatch to the right reader based on file extension. UnknownFormat on miss."""
    p = Path(path)
    reader = EXTENSION_TO_READER.get(p.suffix.lower())
    if reader is None:
        raise UnknownFormat(
            f"unsupported extension {p.suffix!r} for {p}; "
            f"known: {sorted(EXTENSION_TO_READER)}"
        )
    return reader(p)
```

That's the entire reader subsystem. ~25 LOC of public surface.

## Tests

### `tests/test_readers_tabular.py` — synthetic, no external data

For each of `read_csv`, `read_tsv`, `read_parquet`: write a tmp file with 2-3 columns and 3 rows, call the reader, assert column names and row count.

Plus targeted edge cases:

- **CSV with embedded commas in quoted fields** — pandas handles, but assert it.
- **TSV with UTF-8 BOM** in the first column header — many MaxQuant-adjacent exports have this; assert the BOM is stripped or absorbed.
- **Empty file** → either an empty DataFrame or a clear pandas error; pin whichever pandas does so the behaviour is documented.

### `tests/test_readers_dispatch.py` — extension routing

- `.csv` → calls `read_csv`
- `.tsv` → calls `read_tsv`
- `.txt` → calls `read_tsv` (MaxQuant convention)
- `.parquet` → calls `read_parquet`
- `.xyz` → `UnknownFormat`
- Case-insensitive extension (`.CSV`) — accepted

Verify by writing a tmp file per extension and checking the returned DataFrame matches.

### `tests/test_readers_integration.py` — DRY win

One parametrized test, looping over every packaged TOML, that:

1. Loads the TOML with `load_rule(toml_path)`.
2. Looks up the matching test data file by `rule.software_name` in `test_data_download/raw_file_db_downloaded.csv`.
3. If the file exists locally, calls `read_table(data_file)`.
4. Asserts the DataFrame is non-empty and has at least as many columns as `len(rule.columns.var) + len(rule.columns.obs) + len(rule.layers)` (a coarse "did we read something with the expected breadth" check; deeper column-presence checks belong to the converters tests in step 11).

Skip the test (don't fail) when `test_data_download/` is absent — it's a 52 GB cache, gitignored, regenerable via `test_data_download/Makefile`. The skip message points at the makefile.

This is the **single piece of test code** that exercises the reader against all 6 vendors. Adding a 7th vendor TOML automatically extends this test with no additional code.

```python
import csv
from pathlib import Path
import pytest
from anndata_proteomics.readers.dispatch import read_table
from anndata_proteomics.rules.loader import load_rule
from anndata_proteomics.rules.registry import iter_packaged_rules

PROJECT_ROOT = Path(__file__).resolve().parent.parent
TEST_DATA_DIR = PROJECT_ROOT / "test_data_download"
DOWNLOADED_DB = TEST_DATA_DIR / "raw_file_db_downloaded.csv"


def _find_test_data(software_name: str) -> Path | None:
    if not DOWNLOADED_DB.exists():
        return None
    with open(DOWNLOADED_DB) as f:
        for row in csv.DictReader(f):
            if row["software_name"] == software_name and row.get("status") == "ok":
                return TEST_DATA_DIR / "json_dir" / row["input_file_path"]
    return None


@pytest.mark.parametrize(
    "toml_path",
    list(iter_packaged_rules()),
    ids=lambda p: f"{p.parent.name}/{p.name}",
)
def test_reader_loads_test_data_for_packaged_rule(toml_path: Path) -> None:
    rule = load_rule(toml_path)
    data_file = _find_test_data(rule.software_name)
    if data_file is None or not data_file.exists():
        pytest.skip(
            f"no downloaded test data for {rule.software_name!r}; "
            f"regenerate via test_data_download/Makefile"
        )
    df = read_table(data_file)
    expected_min_cols = (
        len(rule.columns.var) + len(rule.columns.obs) + len(rule.layers)
    )
    assert not df.empty, f"{data_file} produced an empty DataFrame"
    assert len(df.columns) >= expected_min_cols, (
        f"{data_file}: got {len(df.columns)} columns, "
        f"rule expects at least {expected_min_cols}"
    )
```

## What this plan deliberately does NOT do

- **No vendor recognition.** That's `converters/recognize.py` (step 6) — given a header, pick the rule.
- **No rule application.** That's `converters/long.py` / `wide.py` (steps 7-8).
- **No column rename / sanitisation.** Readers return columns as-is. Sanitisation per `anndata_omics_bridge/docs/conventions.md` happens in converters.
- **No `[reader]` section on `ParseRule`.** Defer until a real vendor file forces it; YAGNI today.
- **No streaming / chunked reads.** Files in `test_data_download` are mostly < 100 MB; the largest (DIA-NN parquet) is ~22 MB. `pd.read_*` handles them in one shot. Premature optimisation otherwise.
- **No GZ / compressed handling.** None of the 6 packaged vendors emit compressed files. Add later if needed.

## Verification

```bash
cd /Users/wolski/projects/anndata_bridge/anndata_proteomics_bridge
source .venv/bin/activate
uv pip install -e '.[dev]'

# Modules import cleanly
.venv/bin/python -c "from anndata_proteomics.readers.dispatch import read_table, UnknownFormat"
.venv/bin/python -c "from anndata_proteomics.readers.tabular import read_csv, read_tsv, read_parquet"

# Unit tests
.venv/bin/python -m pytest tests/test_readers_tabular.py tests/test_readers_dispatch.py -v

# Integration test against the cache (skips if cache absent)
.venv/bin/python -m pytest tests/test_readers_integration.py -v

# Full suite stays green
.venv/bin/python -m pytest tests/ -q

# Smoke test from Python: read each packaged vendor's first downloaded file
.venv/bin/python - <<'PY'
import csv
from pathlib import Path
from anndata_proteomics.readers.dispatch import read_table
TD = Path("test_data_download")
if (TD / "raw_file_db_downloaded.csv").exists():
    seen = set()
    with open(TD / "raw_file_db_downloaded.csv") as f:
        for row in csv.DictReader(f):
            sw = row["software_name"]
            if sw in seen or row.get("status") != "ok":
                continue
            seen.add(sw)
            p = TD / "json_dir" / row["input_file_path"]
            if p.exists():
                df = read_table(p)
                print(f"{sw:14}  {p.suffix:8}  shape={df.shape}")
PY
```

Expected smoke output (one row per software in the catalog with downloaded test data):

```
DIA-NN          .parquet  shape=(N, 60)
Spectronaut     .tsv      shape=(N, 67)
MaxQuant        .txt      shape=(N, 90+)
FragPipe        .tsv      shape=(N, 48)
PEAKS           .csv      shape=(N, 40)
WOMBAT          .csv      shape=(N, 14)
```

## Commit plan

Single commit: `feat(readers): generic file → DataFrame dispatch (csv / tsv / txt / parquet)`.

Followed by a one-line `docs/RESTART_PLAN.md` update flipping step 5 to ✅, with `09e57bf...` style commit ref. Not a separate commit — fold into the same one.

## After this plan lands

Step 6: `converters/recognize.py` — given a DataFrame's header columns, pick the matching `ParseRule` from the packaged set. This is where header-based vendor recognition lives, distinct from the reader.
