# Integration test data — ProteoBench vendor submissions

APB's integration tests run the real converter over **vendor output that people
uploaded to ProteoBench** (DIA-NN, Spectronaut, MaxQuant, FragPipe, PEAKS,
WOMBAT, and more) spread across ProteoBench's **8 quant modules** (`dda_qexactive`,
`dda_astral`, `dda_peptidoform`, `dia_aif`, `dia_astral`, `dia_diapasef`,
`dia_zenotof`, `dia_singlecell`). "8" counts modules, not vendors — many vendors
appear across those modules, and the download catalog is broader than APB's
packaged parsing rules: only rows whose vendor maps to a packaged rule (DIA-NN,
Spectronaut, MaxQuant, FragPipe, PEAKS, WOMBAT) are fed to the converter; the
rest are indexed but skipped by the fixtures.

That data is large and lives on remote servers, so it is **not stored in the
repo** — it is pulled on demand into a gitignored cache and can be regenerated at
any time.

This document describes how the downloaded test data gets here, how to run the
tests, and how the tests consume the cache.

## Where it lives

Everything sits under `../test_data_download/` (gitignored — see
[`.gitignore`](../.gitignore), all regenerable):

```
test_data_download/
├── raw_file_db_full.csv       # every ProteoBench submission, one row each (from `catalog`)
├── raw_file_db_selected.csv   # one representative row per module+software+version (from `select`)
├── raw_file_db_downloaded.csv # the selected rows plus download results (from `download`) ← the index tests read
├── json_dir/
│   └── <repo_name>/                              # e.g. Results_quant_ion_DIA_AIF
│       ├── <repo_name>-main/<hash>.json          # per-submission metadata (from `catalog`)
│       └── <hash>/                               # one folder per downloaded submission (from `download`)
│           ├── input_file.{tsv,txt,csv,parquet}  # raw vendor output
│           └── param_0..<ext>                     # co-located search parameters (note the double dot)
└── fasta/                     # ProteoBench HYE / HY reference FASTAs (from `fasta`)
```

The three CSVs share the same 8 catalog columns — `module`, `repo_name`,
`intermediate_hash`, `software_name`, `software_version`, `nr_feature`,
`is_temporary`, `old_new`. `selected` differs from `full` only in row count;
`downloaded` adds `input_file_path`, `input_file_size_bytes`, and `status` (one of
`ok`, `not_on_server`, `input_file_missing`).

A few layout notes:

- Inside one `json_dir/<repo_name>/`, the `-main` metadata folder (from `catalog`)
  and the per-`<hash>` download folders (from `download`) live side by side.
  Re-running `catalog` replaces only the `<repo_name>-main` metadata snapshot and
  preserves already-downloaded hash folders.
- Param files are literally named `param_0.` **plus** the extension, e.g.
  `param_0..txt`, `param_0..workflow`, `param_0..xml`, `param_0..yml` — the
  basename between the two dots is empty. [`conftest.py`](conftest.py) pairs each
  input with its params via a `param_0.*` glob, which matches the double-dot names.
- Each `<hash>/` folder is the extracted ProteoBench submission, so it may also
  contain `comment.txt`, `result_performance.csv`, and occasionally
  `input_file_secondary.*` / `param_1..*`. The tool and tests use only
  `input_file.*` and `param_0.*`.

Each submission's `param_0..*` file is downloaded **with** its `input_file.*` and
supplies the **software version**, which selects the packaged parsing-rule variant
— so conversion is param-driven, exactly the way `apb convert --params` picks the
rule from a co-located parameter file. There is no separate parameter download.

## How to (re)generate it

**Prerequisites:** the package installed in your environment (`uv venv && source
.venv/bin/activate && uv pip install -e ".[dev]"`) — this provides the
`apb-testdata` command. Regeneration needs **network access** to GitHub and
`proteobench.cubimed.rub.de`.

Run the commands from any directory; each step feeds the next:

```bash
apb-testdata catalog   # 1. metadata JSONs + raw_file_db_full.csv
apb-testdata select    # 2. representative rows → raw_file_db_selected.csv
apb-testdata download  # 3. vendor files + raw_file_db_downloaded.csv
apb-testdata fasta     # 4. HYE + HY reference FASTAs → fasta/
apb-testdata clean     # remove json_dir/, fasta/, and the three CSVs
```

**How `select` chooses its one row:** it drops submissions with a missing
`nr_feature` (quantified feature count), then for each `(module, software_name,
software_version)` group keeps the row with the **smallest `nr_feature`**, breaking
ties by the lexicographically smallest `intermediate_hash`. Fewest precursors
means the smallest, fastest file to download and convert, so the cache stays small
while still exercising the parser for every module/software/version combination
that has a usable submission. (Dropping missing-`nr_feature` rows can silently
exclude a version that has no precursor count.)

The three data steps are the subcommands of the `apb-testdata` command (source:
[`extract_raw_file_db.py`](../src/anndata_proteomics/scripts/extract_raw_file_db.py)),
which you can run directly. All paths default under `test_data_download/`
regardless of your working directory:

```bash
uv run apb-testdata catalog
uv run apb-testdata select --selection-csv test_data_download/my_selection.csv
uv run apb-testdata download --selection-csv test_data_download/my_selection.csv
uv run apb-testdata fasta
```

`catalog` always refreshes all configured ProteoBench modules. `select` prints the
available module summary and writes `--selection-csv`; use `--module` to restrict
it to one module and `--strategy` to keep all rows, the smallest row per module,
the smallest per module/software, or the smallest per module/software/version.
`download` downloads every row in `--selection-csv`. Use `--catalog-csv`,
`--cache-dir`, and `--manifest-csv` to override their standard paths.
Use `clean --data-dir <folder>` to remove only the known generated artifacts
under a custom test-data root.

`download` is **idempotent**: hashes already present under `json_dir/` are skipped,
so a re-run fetches only what is missing. `catalog` refreshes repository metadata
without removing downloaded fixture directories.

## How to run the tests

```bash
pytest tests/            # runs everything; cache-backed tests skip cleanly if the cache is absent
pytest -m integration    # only the network-marked integration tests
```

The `integration` marker is declared in
[`pyproject.toml`](../pyproject.toml). A full run of the ProteoBench-backed tests
needs the regenerated cache from the steps above.
[`test_cli_integration.py`](test_cli_integration.py) drives the `apb convert` CLI
over the cached inputs and skips per tool when a tool's cached data is missing.

## How the tests consume it

The consumer side lives in the package, not in a catalog service: the test suite
reads `raw_file_db_downloaded.csv` directly and runs the real converter over each
cached vendor file. Everything degrades gracefully — when the gitignored cache is
absent, lookups return `None`/`[]` and the fixtures `pytest.skip(...)`, so the
suite is green on a fresh checkout.

- [`test_data.py`](../src/anndata_proteomics/test_data.py)
  exposes the cache paths (`TEST_DATA_DIR`, `DOWNLOADED_DB`, `FASTA_DIR`) and the
  lookups:
  - `find_test_data(software_name)` — the first successfully-downloaded
    (`status == "ok"`) cached input for a tool, or `None` if the cache index is
    absent or the tool has no cached data. Matches the **exact** catalog
    `software_name` (e.g. `"DIA-NN"`).
  - `find_fasta(module=... | dataset_dir=..., test_data_dir=...)` — the reference
    FASTA for a module, optionally resolved from an explicit cache root.
  - `find_param_file(software_name)` — a representative parameter file for a tool,
    read from the committed in-repo `tests/params/` fixtures (no external checkout).
- [`conftest.py`](conftest.py) turns the index into fixtures. `cached_datasets`
  hands each test the cached inputs for a vendor — the resolved input path, its
  paired param file, the detected version, the column headers, and the
  quantification levels that input can be converted to — so tests run the converter
  without knowing cache paths. Its argument is a **lowercased vendor slug** (e.g.
  `"diann"`), unlike `find_test_data`'s exact catalog name. `spectronaut_datasets`
  is the list of cached Spectronaut datasets; `diann_full_subset` is a single
  prepared DIA-NN dataset (`{df, version, slug}`) requiring ion+protein+fragment
  levels, which skips if none is cached.
