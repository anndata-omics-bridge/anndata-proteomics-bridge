# Plan: Modification Handling And ProForma/SDRF Normalization

## Goal

Make APB the upstream owner of reusable proteomics modification handling.

Modification handling is needed in two related places:

1. **Parameter files**: searched fixed/variable modifications from vendor
   parameter files should become structured metadata that can be exported to
   SDRF-Proteomics `comment[modification parameters]`.
2. **Quantification result columns**: modified peptide / peptidoform columns in
   vendor result files should be normalized so APB can build stable `var`
   identities and optionally emit ProForma strings when localization is known.

ProteoBench currently carries per-tool modification parsing rules in parsing
TOMLs. Those rules should migrate into APB parsing TOMLs/schema rather than
staying downstream.

## Standards Context

### SDRF-Proteomics

SDRF-Proteomics uses repeated `comment[modification parameters]` columns for
searched modifications. The value format is an SDRF key=value object, not
ProForma.

Example:

```text
NT=Oxidation;AC=UNIMOD:35;TA=M;MT=variable;PP=Anywhere
```

Fields:

- `NT`: name term, e.g. `Oxidation`
- `AC`: accession, e.g. `UNIMOD:35`
- `TA`: target amino acid or site, e.g. `M`
- `MT`: modification type, e.g. `fixed` or `variable`
- `PP`: position qualifier, e.g. `Anywhere`

Parse order-insensitively, but export in canonical order:

```text
NT=<name>;AC=<accession>;MT=<fixed|variable>;TA=<target>;PP=<position>
```

### ProForma

ProForma is appropriate for modified peptide/proteoform sequences. It is the
right output format for localized peptidoforms from quantification result
columns, but it should not be the only internal representation because searched
parameter-file modifications are often rules, not concrete sequence
occurrences.

## Current APB Gap

APB parsing-rule schema currently has no `modifications_parser` section.
Current APB TOMLs preserve vendor modified sequence columns as `var` columns,
for example:

- `parsing_rules/diann/parse_diann_ion_1.toml`
- `parsing_rules/fragpipe/parse_fragpipe_ion_1.toml`
- `parsing_rules/maxquant/parse_maxquant_ion_1.toml`
- `parsing_rules/peaks/parse_peaks_ion_1.toml`
- `parsing_rules/spectronaut/parse_spectronaut_ion_1.toml`
- `parsing_rules/wombat/parse_wombat_peptidoform_1.toml`

These files need a first-class rule section for modification parsing and
normalization.

## Current ProteoBench Source Rules

ProteoBench TOMLs with `[modifications_parser]` live under:

```text
ProteoBench/proteobench/io/parsing/io_parse_settings/
```

There are 60 TOML files with modification parser configuration. Grouped by tool:

| Tool | Files | Unique dictionaries | Parse column(s) | Pattern(s) |
|---|---:|---:|---|---|
| `adanovo` | 1 | 1 | `sequence` | `(?:^([+-]\d+\.\d+))|(\([+-]\d+\.\d+\))` |
| `alphapept` | 2 | 1 | `Modified sequence` | `([a-z]+)` |
| `casanovo` | 1 | 1 | `sequence` | `([\d+-.]+)` |
| `deepnovo` | 1 | 1 | `sequence` | `\(([^)]+)\)` |
| `diann` | 8 | 1 | `Sequence` | `\(([^()]*)\)` |
| `fragpipe` | 2 | 1 | `Modified Sequence` | `(?<=\[).+?(?=\])` |
| `fragpipe_DIA` | 5 | 2 | `Modified Sequence` | `(?<=\[).+?(?=\])`, `\[([^]]+)\]` |
| `i2massChroQ` | 2 | 1 | `proforma` | `(?<=\[).+?(?=\])` |
| `maxdia` | 5 | 2 | `Modified sequence` | `\([^()]*\)|\([^()]*\([^()]*\)[^()]*\)` |
| `maxquant` | 2 | 1 | `Modified sequence` | `\([^()]*\)|\([^()]*\([^()]*\)[^()]*\)` |
| `metamorpheus` | 2 | 1 | `Modified sequence` | `\[(.*?)\]` |
| `msaid` | 5 | 1 | `Sequence` | `\[(.*?)\]` |
| `peaks` | 9 | 1 | `Sequence` | `(?<=\().+?(?=\))` |
| `pepnet` | 1 | 1 | `sequence` | `\(([^)]+)\)` |
| `pihelixnovo` | 1 | 1 | `sequence` | `([\d+-.]+)` |
| `piprimenovo` | 1 | 1 | `sequence` | `([\d+-.]+)` |
| `sage` | 3 | 1 | `Sequence` | `(?<=\[).+?(?=\])`, `\[([^]]+)\]` |
| `spectronaut` | 6 | 1 | `Sequence` | `\[(.*?)\]` |
| `wombat` | 3 | 1 | `Sequence` | `(?<=\[).+?(?=\])` |

Initial APB migration should cover only tools with packaged APB rules:

- DIA-NN: ProteoBench `parse_settings_diann.toml`
- FragPipe: ProteoBench `parse_settings_fragpipe.toml` and relevant
  `parse_settings_fragpipe_DIA.toml`
- MaxQuant / MaxDIA: ProteoBench `parse_settings_maxquant.toml` and
  `parse_settings_maxdia.toml`
- PEAKS: ProteoBench `parse_settings_peaks.toml`
- Spectronaut: ProteoBench `parse_settings_spectronaut.toml`
- WOMBAT: ProteoBench `parse_settings_wombat.toml`

De novo tools and APB-unsupported quant tools should be inventoried but deferred
unless explicitly added to scope.

## Proposed APB Schema Additions

Add an optional top-level `[modifications]` section to APB parsing-rule TOMLs.

Suggested schema:

```toml
[modifications]
source_column = "Modified.Sequence"
sequence_column = "Stripped.Sequence"
output_column = "proforma_sequence"
parser = "token_regex"
token_pattern = "\\(([^()]*)\\)"
token_position = "after_residue"
case_sensitive = false
unknown_policy = "preserve"

[[modifications.map]]
token = "(unimod:35)"
name = "Oxidation"
accession = "UNIMOD:35"
target = "M"

[[modifications.map]]
token = "(unimod:4)"
name = "Carbamidomethyl"
accession = "UNIMOD:4"
target = "C"
```

Schema fields:

- `source_column`: vendor column to parse.
- `sequence_column`: optional stripped sequence column used for validation.
- `output_column`: normalized sequence column to write into the DataFrame before
  AnnData assembly, default `proforma_sequence`.
- `parser`: parsing mode. Initial values:
  - `token_regex`: extract vendor tokens with regex and map tokens.
  - `already_proforma`: validate/canonicalize existing ProForma.
  - `separate_mod_column`: use base sequence plus separate modification column.
- `token_pattern`: regex used by `token_regex`.
- `token_position`: how tokens attach to residues:
  - `before_residue`
  - `after_residue`
  - `n_term`
  - `c_term`
  - `embedded`
  - `unknown`
- `case_sensitive`: whether token matching is case-sensitive.
- `unknown_policy`: `preserve`, `drop`, or `error`.
- `modifications.map`: token-to-normalized-modification mapping.

The mapping is the core of the migration. Do not migrate ProteoBench entries as
plain `token -> name` dictionaries. Migrate them as structured identity records:

- `token`: exact vendor token or mass token found in the result file
- `name`: canonical display name
- `accession`: controlled vocabulary accession, preferably Unimod when known
- `target`: residue or terminus target when known
- `position`: localization constraint such as `Anywhere`, `Protein N-term`,
  `Peptide N-term`, `C-term`, or `unknown`
- `mass_delta`: optional numeric delta for mass-token mappings

Use **Unimod** as the preferred accession vocabulary for search-engine
modifications because it is also the practical SDRF-Proteomics target for these
examples. Preserve PSI-MOD accessions when the source data already provides them
or when no Unimod equivalent is known. Do **not** map directly to UniProt as the
primary target; UniProt features may be useful downstream, but they are not the
right interchange target for search settings or ProForma peptidoforms.

Example FragPipe migration from ProteoBench:

```toml
[modifications]
source_column = "Modified Sequence"
parser = "token_regex"
token_pattern = "(?<=\\[).+?(?=\\])"
token_position = "after_residue"
case_sensitive = false
unknown_policy = "preserve"
output_column = "proforma_sequence"

[[modifications.map]]
token = "57.0215"
name = "Carbamidomethyl"
accession = "UNIMOD:4"
target = "C"
position = "Anywhere"
mass_delta = 57.0215

[[modifications.map]]
token = "57.0216"
name = "Carbamidomethyl"
accession = "UNIMOD:4"
target = "C"
position = "Anywhere"
mass_delta = 57.0216

[[modifications.map]]
token = "15.9949"
name = "Oxidation"
accession = "UNIMOD:35"
target = "M"
position = "Anywhere"
mass_delta = 15.9949

[[modifications.map]]
token = "-17.026548"
name = "Gln->pyro-Glu"
accession = "UNIMOD:28"
target = "Q"
position = "N-term"
mass_delta = -17.026548

[[modifications.map]]
token = "-18.010565"
name = "Glu->pyro-Glu"
accession = "UNIMOD:27"
target = "E"
position = "N-term"
mass_delta = -18.010565

[[modifications.map]]
token = "42.0106"
name = "Acetyl"
accession = "UNIMOD:1"
target = "N-term"
position = "N-term"
mass_delta = 42.0106
```

This structured mapping supports multiple exports:

- ProForma for localized peptidoform strings, e.g. `M[UNIMOD:35]PEPTIDE`.
- SDRF searched-modification metadata, e.g.
  `NT=Oxidation;AC=UNIMOD:35;MT=variable;TA=M;PP=Anywhere`.
- Internal AnnData provenance in `uns` without losing the original vendor token.

One caveat: quantification result columns usually identify observed/localized
modification occurrences, but do not reliably say whether a modification was
searched as fixed or variable. `MT=fixed|variable` should therefore come from
parameter-file parsing when available. Result-column normalization can keep
`mod_type = "unknown"` unless the tool rule has authoritative knowledge.

Keep old ProteoBench names (`before_aa`, `isalpha`, `isupper`) out of the APB
public schema unless they describe a real domain concept. They are parser
implementation details and should be replaced with clearer fields such as
`token_position` and `case_sensitive`.

## Internal Data Model

Use separate models for searched modifications and sequence occurrences.

### SearchedModification

For parameter files and SDRF export:

- `name`
- `accession`
- `mod_type`: `fixed`, `variable`, or `unknown`
- `target`
- `position`
- `mass_delta`
- `source`

This model exports to SDRF key=value strings.

### ModificationOccurrence

For modified peptide/peptidoform sequences:

- `name`
- `accession`
- `target_residue`
- `sequence_index`
- `position`
- `mass_delta`
- `source_token`

This model exports to ProForma when localization is known.

### ModifiedSequence

For quantification results:

- `stripped_sequence`
- `proforma_sequence`
- `occurrences`
- `source_sequence`
- `unknown_tokens`

## Implementation Plan

1. **Inventory ProteoBench TOMLs**
   - Extract all `[modifications_parser]` sections.
   - Deduplicate mapping dictionaries by tool.
   - Decide which mappings belong in APB packaged TOMLs now and which remain
     deferred.

2. **Extend APB TOML schema**
   - Add `ModificationRule`, `ModificationMapEntry`, and `Modifications`
     Pydantic models.
   - Add validation:
     - `token_regex` requires `source_column`, `token_pattern`, and mappings.
     - `already_proforma` requires `source_column` but no token map.
     - `separate_mod_column` requires source columns defined by the rule.
     - `unknown_policy` must be explicit.

3. **Add normalization code**
   - New module: `src/anndata_proteomics/modifications/`.
   - Suggested files:
     - `model.py`
     - `sdrf.py`
     - `proforma.py`
     - `vendor_tokens.py`
     - `apply_rules.py`
   - Reuse established libraries if available and practical; do not implement a
     full ProForma parser from scratch unless no dependency fits the needed
     scope.

4. **Apply normalization before conversion**
   - In the conversion pipeline, after `read_table` and before
     `recognize/convert`, apply modification normalization for the matched rule.
   - Add normalized columns, do not silently replace source vendor columns.
   - Update `axis.var_keys` for APB rules to use normalized ProForma columns
     where possible.

5. **Migrate APB packaged TOMLs**
   - DIA-NN: map `(unimod:35)`, `(unimod:1)`, `(unimod:4)`.
   - MaxQuant / MaxDIA: map `(ox)`, `(ac)`, `(oxidation (m))`,
     `(acetyl (protein n-term))`, and related MaxDIA variants.
   - FragPipe / FragPipe DIA: map numeric mass tokens such as `57.0215`,
     `15.9949`, `-17.026548`, `-18.010565`, `42.0106`.
   - PEAKS: map `+57.02`, `+15.99`, `-17.026548`, `-18.010565`, `+42.01`.
   - Spectronaut: map `[oxidation (m)]`, `[acetyl (protein n-term)]`,
     `[carbamidomethyl (c)]`.
   - WOMBAT: map numeric tokens such as `160`, `147`, `80`, `111`, `43`.

6. **Connect parameter parsing**
   - Reuse the same canonical modification map for parameter parser outputs.
   - Parameter parsers should emit `SearchedModification` objects.
   - Add SDRF export helpers that render `comment[modification parameters]`
     values from `SearchedModification`.

7. **Backfill tests**
   - Unit-test vendor token parsing for each APB-supported tool.
   - Unit-test ProForma output on known examples.
   - Unit-test SDRF export for fixed and variable searched modifications.
   - Add converter tests proving `var` identity uses normalized modified
     sequence columns where configured.
   - Add regression tests against ProteoBench fixtures where available.

8. **ProteoBench compatibility**
   - After APB behavior is stable, make ProteoBench consume APB modification
     normalization instead of its local TOML parser logic.
   - Keep a compatibility adapter during transition.
   - Delete downstream duplication only after APB and ProteoBench tests agree.

## Design Decisions To Resolve Before Coding

- Should APB use `ProForma` as the canonical `var` modified-sequence column for
  all tools where localization is available?
- How should ambiguous mass-only tokens be handled when a mass could map to
  multiple modifications?
- Should APB canonicalize to Unimod accessions when known and preserve PSI-MOD
  accessions only when Unimod is unavailable?
- Should unknown modification tokens make conversion fail by default, or should
  they be preserved with warnings?
- Should de novo modification parsing be part of this migration or a later
  migration?

## Verification

Targeted APB checks:

```bash
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider \
  tests/test_modifications_*.py \
  tests/test_rule_validate.py \
  tests/test_converters_e2e.py
```

Compatibility checks after ProteoBench delegation:

```bash
pytest ProteoBench/test/test_parse_params_*.py
pytest ProteoBench/test/test_*parsing*.py
```

## References

- SDRF-Proteomics specification:
  https://sdrf.quantms.org/specification.html
- SDRF terms reference:
  https://sdrf.quantms.org/sdrf-terms.html
- HUPO-PSI ProForma:
  https://www.psidev.info/proforma
