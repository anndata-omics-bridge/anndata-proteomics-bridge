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
[`.gitignore`](../.gitignore), which pegs it at ~52 GB, all regenerable):

```
test_data_download/
├── Makefile                   # thin wrapper: `make catalog/select/download` → `apb-testdata` (+ a `fasta` step)
├── raw_file_db_full.csv       # every ProteoBench submission, one row each (from `catalog`)
├── raw_file_db_selected.csv   # one representative row per module+software+version (from `select`)
├── raw_file_db_downloaded.csv # the selected rows plus download results (from `download`) ← the index tests read
├── json_dir/
│   └── <repo_name>/                              # e.g. Results_quant_ion_DIA_AIF
│       ├── <repo_name>-main/<hash>.json          # per-submission metadata (from `catalog`)
│       └── <hash>/                               # one folder per downloaded submission (from `download`)
│           ├── input_file.{tsv,txt,csv,parquet}  # raw vendor output
│           └── param_0..<ext>                     # co-located search parameters (note the double dot)
└── fasta/                     # ProteoBench HYE / HY reference FASTAs (from `make fasta`)
```

The three CSVs share the same 8 catalog columns — `module`, `repo_name`,
`intermediate_hash`, `software_name`, `software_version`, `nr_prec`,
`is_temporary`, `old_new`. `selected` differs from `full` only in row count;
`downloaded` adds `input_file_path`, `input_file_size_bytes`, and `status` (one of
`ok`, `not_on_server`, `input_file_missing`).

A few layout notes:

- Inside one `json_dir/<repo_name>/`, the `-main` metadata folder (from `catalog`)
  and the per-`<hash>` download folders (from `download`) live side by side.
  Re-running `catalog` for a repo **wipes that whole `<repo_name>/` folder** first
  (including already-downloaded hash folders), so re-catalog then re-download.
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
`apb-testdata` command that the `make` targets and the direct calls below both
use. Regeneration needs **network access** to GitHub and
`proteobench.cubimed.rub.de`.

Run the pipeline from the cache directory; each step feeds the next:

```bash
cd test_data_download
make catalog     # 1. pull datapoint JSONs from the 8 ProteoBench Results_quant_* GitHub
                 #    repos → json_dir/<repo>/<repo>-main/  +  raw_file_db_full.csv
make select      # 2. reduce to one representative row per (module, software, version)
                 #    → raw_file_db_selected.csv
make download    # 3. download the raw input_file.* for the selected rows
                 #    → json_dir/<repo>/<hash>/  +  raw_file_db_downloaded.csv
make fasta       # 4. download the HYE + HY reference FASTAs → fasta/
make clean       # remove json_dir/, fasta/, and the three CSVs
```

`make help` lists the targets.

**How `select` chooses its one row:** it drops submissions with a missing
`nr_prec` (precursor count), then for each `(module, software_name,
software_version)` group keeps the row with the **smallest `nr_prec`**, breaking
ties by the lexicographically smallest `intermediate_hash`. Fewest precursors
means the smallest, fastest file to download and convert, so the cache stays small
while still exercising the parser for every module/software/version combination
that has a usable submission. (Dropping missing-`nr_prec` rows can silently
exclude a version that has no precursor count.)

The three data steps are the subcommands of the `apb-testdata` command (source:
[`extract_raw_file_db.py`](../src/anndata_proteomics/scripts/extract_raw_file_db.py)),
which you can run directly for finer control (e.g. a subset of modules). Outputs
default under `test_data_download/` regardless of your working directory; pass
`--json-dir test_data_download/json_dir` to `catalog` so it doesn't fall back to
its own default of `test_data_download/temp_results`, which the rest of the
pipeline never reads (and which `make clean` does not remove):

```bash
uv run apb-testdata catalog --json-dir test_data_download/json_dir --modules dia_aif dda_astral
uv run apb-testdata select
uv run apb-testdata download
```

There is **no `fasta` subcommand** — the FASTAs are fetched by `make fasta` only
(curl + unzip in the [`Makefile`](../test_data_download/Makefile)).

`download` is **idempotent**: hashes already present under `json_dir/` are skipped,
so a re-run fetches only what is missing. `catalog` is the one step that wipes and
overwrites (per repo).

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
(and the report generator, [`generate_report.py`](../tools/generate_report.py))
read `raw_file_db_downloaded.csv` directly and run the real converter over each
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
  - `find_fasta(module=... | dataset_dir=...)` — the reference FASTA for a module.
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
