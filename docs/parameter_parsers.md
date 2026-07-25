# Vendor Parameter Parsers

Vendor search-parameter parsers live under
`src/anndata_proteomics/params/parsers/`. They read whole-experiment search
settings such as enzyme, tolerances, FDR thresholds, precursor ranges, and
fixed/variable modifications.

Generic parameter code remains one level higher:

```text
params/
  anndata_io.py
  model.py
  registry.py
  parsers/
```

## Contract

Every vendor parser exposes:

```python
extract_params(source) -> Parameters
```

`params.registry` dispatches by software name:

```python
parse_params(path, software)
get_parser("DIA-NN")
available_software()
```

Registered parsers:

| Key(s) | Implementation |
|---|---|
| `alphapept` | `params.parsers.alphapept` |
| `dia-nn`, `diann` | `params.parsers.diann` |
| `fragpipe` | `params.parsers.fragpipe` |
| `maxquant` | `params.parsers.maxquant` |
| `metamorpheus` | `params.parsers.metamorpheus` |
| `msaid` | `params.parsers.msaid` |
| `peaks` | `params.parsers.peaks` |
| `sage` | `params.parsers.sage` |
| `spectronaut` | `params.parsers.spectronaut` |
| `wombat` | `params.parsers.wombat` |

## Parser Pattern

Each parser does three things:

1. Read the source file.
2. Extract vendor-specific fields.
3. Build `Parameters(**fields)`.

Shared file-reading helpers live in `params.parsers._common`:

- `read_text`
- `read_lines`
- `format_tolerance_range`
- `homogenize_paren_mods`
- `lookup_mass_mod`

The `Parameters` model then applies shared normalization: enzyme names, FDR
probabilities, mass tolerances, missing values, fixed/variable modifications,
and min/max ranges.

`Parameters.acquisition_method` is a non-nullable
`"DDA" | "DIA" | "unknown"` value. It defaults to `"unknown"` so parameter
payloads written before the field existed remain readable. DIA-NN reports
`"DDA"` when its command line contains `--dda` or its log states
`All runs will be analysed as DDA runs`; otherwise a valid DIA-NN parameter
file reports `"DIA"`, the software default. Other parsers retain `"unknown"`
until their acquisition markers are explicitly supported.

Rule documents may consume this and other typed search parameters while
materializing a level. Parameter parsing selects the effective conversion
rule; converters still receive an ordinary flat `ParseRule`.

## Input Families

| Vendor | Input | Parser style | Modification style |
|---|---|---|---|
| AlphaPept | YAML | dict access | flat token map |
| DIA-NN | log / cfg text | command line + regex + cfg block | flat token map |
| FragPipe | `.workflow` | key/value series | mass lookup |
| MaxQuant | `mqpar.xml` | XML flattening | `Name (Residue)` text |
| MetaMorpheus | TOML + version text | paired-file loader | `Name on Residue` text |
| MSAID | table | row dict | flat token map |
| PEAKS | text report | label scan | flat token map |
| Sage | JSON | nested dict | mass lookup |
| Spectronaut | text export | label/regex scan | `Name (Residue)` text |
| WOMBAT | YAML | dict access | `Name of Residue` text |

## Modification Handling

Parsers emit `fixed_mods` and `variable_mods` in ProForma-like token form
(`C[Carbamidomethyl]`, `M[Oxidation]`). Three mechanics are used:

- flat token maps for vendor names or UniMod-like tokens;
- numeric mass lookup for tools that report mass deltas;
- text homogenizers for `Name <separator> Residue` formats.

The mapping data remains vendor-specific. Shared mechanics belong in
`params.parsers._common`.

## Tests

`tests/params/` holds parser inputs and ProteoBench-style expected CSV files.
The parser suites compare `extract_params(...).to_series()` with
`Parameters.from_series(expected).to_series()`, so both sides pass through the
same model normalization. The literal acquisition value `"unknown"` is
field-specific data and is preserved even though that token means “missing”
for older free-text parameter fields.

Current focused coverage: `57` parameter tests.

## Adding A Parser

1. Add `params/parsers/<vendor>.py` with
   `extract_params(source) -> Parameters`.
2. Register it in `params/registry.py`.
3. Add input and expected CSV fixtures under `tests/params/`.
4. Add or extend `tests/test_params_<vendor>.py`.
