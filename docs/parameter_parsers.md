# Vendor parameter parsers

How the ~10 vendor parameter-file parsers under `src/anndata_proteomics/params/` work and how
they share the common machinery. For the package-level picture and diagrams see
[parsing_architecture.md](parsing_architecture.md) (the params class diagram and flow **A**).

## The shared contract

Every vendor module exposes the **same entry point** and produces the **same type**:

```python
extract_params(source) -> Parameters     # source = path | text file-like | bytes file-like
```

`registry.py` dispatches by software name (case-insensitive):

```python
parse_params(path, software)   # = get_parser(software)(path)
get_parser("DIA-NN")           # -> diann.extract_params   (keys: "dia-nn", "diann")
available_software()           # -> sorted list of keys
```

Registered keys → module:

| key(s) | module |
|---|---|
| `alphapept` | `alphapept.py` |
| `dia-nn`, `diann` | `diann.py` |
| `fragpipe` | `fragpipe.py` |
| `maxquant` | `maxquant.py` |
| `metamorpheus` | `metamorpheus.py` |
| `msaid` | `msaid.py` |
| `peaks` | `peaks.py` |
| `sage` | `sage.py` |
| `spectronaut` | `spectronaut.py` |
| `wombat` | `wombat.py` |

Every parser follows the same three steps: **read the source → vendor-specific parse → build
`Parameters(**fields)`**. Source reading is centralized in `_common.read_text` / `read_lines`
for text and structured-text formats (XML and CSV use `ElementTree` / pandas instead). The
final `Parameters(...)` construction runs the model's validators, so vendors only need to
produce *raw-ish* field values — normalization is shared (see "Shared normalization" below).

## Per-vendor overview

| Vendor | Input file (example fixture) | Format | How it's read | Parse technique | Mod family |
|---|---|---|---|---|---|
| **AlphaPept** | `alphapept_0.4.9.yaml` | YAML | `yaml.safe_load(read_text)` | dict key access (`summary`/`fasta`/`search`/`features`) | flat |
| **DIA-NN** | `DIANN_…report.log.txt`, `DIANN_cfg_*.txt` | log / cfg text | `read_lines` | command-line settings + in-log regex + optional cfg block | flat |
| **FragPipe** | `fragpipe.workflow` | key=value workflow | `read_text` → `_read_workflow` → pandas `Series` | field derivations on the Series | mass-lookup |
| **MaxQuant** | `mqpar_MQ2.1.3.0_noMBR.xml` | XML | `_read_xml` (`ElementTree`, flattened to dict) | version-gated MultiIndex paths | text `Name (Residue)` |
| **MetaMorpheus** | `…_search_task_config.toml` **+** `…_version_result.txt` | TOML **+** version text (two files) | `read_text(errors="replace")` + `tomllib.loads` per file | format-trial picks TOML vs version line | text `Name on Residue` |
| **MSAID** | `MSAID_default_params.*` | tabular | `pd.read_csv` | row dict | flat |
| **PEAKS** | `PEAKS_parameters*.txt` | text report | `read_lines(strip=True)` | label / regex line scan | flat |
| **Sage** | `sage_parameterfile.json` | JSON | `json.loads(read_text)` | nested dict (`database.enzyme`, …) | mass-lookup |
| **Spectronaut** | `Spectronaut_*.txt` | settings-export text | `read_lines(strip=True)` | label / regex line scan | text `Name (Residue)` |
| **WOMBAT** | `wombat_params.yaml` | YAML | `yaml.safe_load(read_text)` | dict key access (`params`) | text `Name of Residue` |

## Input-format families

- **Structured documents** — read whole, then indexed by key: YAML (AlphaPept, WOMBAT),
  JSON (Sage), XML (MaxQuant), TOML (MetaMorpheus), tabular/CSV (MSAID). Acquisition is a
  one-liner (`read_text` then a library `load`, or pandas/ElementTree); the work is field mapping.
- **Line/regex text logs** — DIA-NN, PEAKS, Spectronaut. Read with `read_lines` and scanned by
  label lookups (`_value`) or regex. PEAKS/Spectronaut strip each line.
- **key=value workflow** — FragPipe `.workflow` is flattened into a pandas `Series` indexed by
  property name, then individual settings are derived from it.

## Modification handling — three families

All parsers ultimately emit `fixed_mods` / `variable_mods` as comma-joined **ProForma-style**
tokens (e.g. `C[Carbamidomethyl]`, `M[Oxidation]`). They get there three different ways:

1. **Flat token map** — `MAP.get(token, token)` over a per-vendor dict.
   *AlphaPept* (`cC`→`C[Carbamidomethyl]`), *DIA-NN* (UniMod + `Carbamidomethyl (C)` forms),
   *MSAID*, *PEAKS* (`Carbamidomethyl (+57.02)` forms).
2. **Numeric mass lookup** — `_lookup_mod_name(mass)` matches a mass shift within `1e-3` against
   a `*MASS_TO_MOD*` table. *FragPipe* (`57.02146`→`Carbamidomethyl`), *Sage*.
3. **Text `Name <sep> Residue` homogenizer** — `_homogenize_mod*` splits a descriptive spec into
   name + residue/terminus and renders `Residue[Name]`. *MaxQuant* (`Name (Residue)`),
   *MetaMorpheus* (`Name on Residue`), *Spectronaut* (`Name (Residue)`), *WOMBAT* (`Name of Residue`).

> These three mechanics are currently re-implemented per vendor (the mapping *data* is genuinely
> per-vendor, but the *parsing mechanics* are duplicated). Consolidating them into `_common` is
> tracked as **Action 15** in [TODO/TODO_code_review_june.md](../TODO/TODO_code_review_june.md).

## Shared normalization (applies to every vendor)

Because all parsers return a `Parameters` object, the model's validators normalize uniformly,
so individual parsers don't repeat this:

- **enzyme** → canonical name via `_ENZYME_MAP` (e.g. `kr` → `Trypsin/P`).
- **FDR** (`ident_fdr_*`) → `Probability` in `[0,1]`; a numeric `>= 1` is treated as a percentage
  and divided by 100.
- **tolerances** → `MassTolerance.parse` (handles `20 ppm`, `[-20 ppm, 20 ppm]`, and
  auto-calibration sentinels like `dynamic` / `0 ppm`).
- **missing values** → sentinel strings (`""`, `none`, `n/a`, `placeholder`, `-`, …) become `None`.
- **fixed/variable mods** → split into `list[SearchedModification]`.
- **ranges** → `min_* <= max_*` enforced for charge, peptide length, and m/z.

`Parameters.to_series()` / `from_series()` round-trip to the ProteoBench CSV layout used by the
equivalence tests.

## Per-vendor specifics (the non-obvious bits)

- **DIA-NN** — layered precedence: `_DIANN_IMPLICIT_DEFAULTS` ← command-line `--` settings ←
  in-log regex fallbacks ← optional cfg block. Enzyme strings like `K*,R*` → `Trypsin/P`,
  `K*,R*,!P` → `Trypsin`.
- **FragPipe** — FDR/MBR depends on whether the run used DIA-NN (`diann.run-dia-nn`) or the
  Philosopher report; precursor/fragment tolerance units from `msfragger.*_mass_units`; m/z from
  digest mass + charge.
- **MaxQuant** — `mqpar.xml` flattened to a MultiIndex; the fixed/variable-mods path is
  version-gated (different XML location for MaxQuant > 1.6.0.0).
- **MetaMorpheus** — the only parser taking **two** inputs (a search-task TOML and a version-text
  file); `_load_pair` tries to parse each as TOML and treats the non-TOML one as the version line.
- **Sage** — enzyme from `database.enzyme` (`cleave_at` + `restrict`: with `restrict="P"` →
  `Trypsin`); static/variable mods are `residue → mass(es)` maps.
- **PEAKS / Spectronaut** — plain text exports scanned by labels; PEAKS strips a trailing `%`
  from FDR values; Spectronaut "dynamic"/"0 ppm" tolerances become "Automatic calibration".
- **AlphaPept / WOMBAT** — YAML configs read by key; tolerances assembled as bracketed
  `[-x unit, x unit]` strings from a value + `ppm`/`Da` flag.
- **MSAID** — parameters supplied as a small table, read with pandas.

## Tests & fixtures

`tests/params/` holds, per vendor, an **input file** plus an **expected CSV** (the ProteoBench
golden output). The `tests/test_params_*.py` suites call `extract_params(input).to_series()` and
compare it to `Parameters.from_series(expected_csv).to_series()` — so both sides pass through the
same model normalization. (56 param tests; see also `test_params_model.py` for the model itself.)

## Adding a vendor

1. Add `params/<vendor>.py` with `extract_params(source) -> Parameters` — read via
   `_common.read_text`/`read_lines` (or pandas/ElementTree), then build `Parameters(...)`.
2. Register it in `registry._REGISTRY`.
3. Add an input + expected-CSV fixture under `tests/params/` and a `test_params_<vendor>.py`.

(Keeping `registry`, packaged `parsing_rules/`, and fixtures in sync is tracked as **Action 8**
in the code review.)
