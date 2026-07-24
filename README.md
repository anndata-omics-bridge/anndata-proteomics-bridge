# anndata-proteomics (APB)

Convert proteomics quantification output into **AnnData / MuData** using declarative JSON rules.

**Documentation:** <https://anndata-omics-bridge.github.io/anndata-proteomics-bridge/>

- **Declarative, not bespoke.** Every vendor × quantification-level is a small JSON rule shipped inside the package — adding or fixing a converter means editing a `.json`, not writing tool-specific Python.
- **One file → a multi-level MuData.** A single vendor export is converted into a MuData whose modalities are the quantification levels it provides (`ion` / `fragment` / `peptidoform` / `protein`) on a shared run axis — or a single-level AnnData when you ask for one level.
- **Standardised content.** Peptide modifications are normalised to **ProForma**; a per-vendor parser reads the vendor **parameter file** (enzyme, FDR, tolerances, …) into one typed record under `uns['search_parameters']`.
- **Enrichable and validated.** Join sample metadata onto `obs` (`apb annotate`); `apb fasta` adds protein annotation and automatically checks every peptide-derived feature against the supplied FASTA with Aho--Corasick.
- **ProteoBench scoring.** `apb proteobench` computes HYE intermediate statistics and compatible score JSON directly from converted ion or peptidoform matrices.
- **Interoperable.** Writes plain `.h5ad` / `.h5mu`, readable from Python (`anndata` / `mudata` / `scanpy`) and R (`anndataR`).

> **New to AnnData?**  It's the standard container for an annotated data matrix — observations (`obs`, here MS runs) × variables (`var`, here peptides/proteins), with multiple measurement `layers`, dimensionality-reduction slots (`obsm`/`varm`), and free-form metadata (`uns`). **MuData** bundles several AnnData objects as *modalities*. See [anndata.readthedocs.io](https://anndata.readthedocs.io) and [mudata.readthedocs.io](https://mudata.readthedocs.io).

**Six vendors, four quantification levels:**

| Vendor | Level | Shape | Version |
|---|---|---|---|
| DIA-NN | ion | long | 1.x, 2.x |
| DIA-NN | fragment | long | 1.x |
| DIA-NN | protein | long | 1.x, 2.x |
| Spectronaut | ion | long | 19.x, 20.x |
| Spectronaut | fragment | long | 19.x, 20.x |
| Spectronaut | protein | long | 19.x, 20.x |
| MaxQuant | ion | long | 2.6.7.0 |
| FragPipe | ion | wide | 22.1-build02 |
| PEAKS | ion | wide | 13 |
| WOMBAT | peptidoform | wide | 0.9.11 |

*Shape* = how the vendor lays out the table: **long** (one row per run × feature) or **wide** (samples as columns, one row per feature). *Version* is matched against the software version parsed from the parameter file; DIA-NN ships version-specific rules (`v1/`, `v2/`).

## Inputs per format

Each conversion takes the **quant file** (the numbers) plus the **parameter file** (the search settings). A per-vendor parser reads the parameter file into one typed record — enzyme, FDR, tolerances, modifications, match-between-runs, … — stored under `uns['anndata_proteomics']['search_parameters']`. It is required because it does double duty: the **software version** it reports selects the rule variant (e.g. DIA-NN v1 vs v2), and the **enzyme** it reports drives `apb fasta`'s theoretical peptide counts.

| Vendor | Quant file | Parameter file (required) |
|---|---|---|
| DIA-NN | `report.tsv` / `report.parquet` | `report.log.txt` (DIA-NN run log) |
| Spectronaut | long report export (`.tsv`) | `…ExperimentSetupOverview….txt` (setup export) |
| MaxQuant | `evidence.txt` | `mqpar.xml` |
| FragPipe | `combined_modified_peptide.tsv` | `fragpipe.workflow` |
| PEAKS | exported peptide CSV (`.csv`) | `parameters.txt` (PEAKS parameter export) |
| WOMBAT | standardised output (`.csv`) | `params.yaml` (WOMBAT-P) |

> The parameter parser supports more search engines than there are conversion rules today (also Sage, AlphaPept, MetaMorpheus, MSAID) — those can be paired with a conversion rule as the rules land.

## Install

```bash
uv sync --frozen
```

For development, install all locked quality and documentation tools:

```bash
uv sync --frozen --extra dev --group docs
uv run pre-commit run --hook-stage pre-commit --all-files
uv run pre-commit run --hook-stage pre-push --all-files
```

See [docs/development.md](docs/development.md) for the security audit and
individual checks.

## Integration test data

`test_data_download/` is the single local ProteoBench cache used by APB's
integration tests. Build it with the dedicated CLI:

```bash
apb-testdata catalog
apb-testdata select
apb-testdata download
apb-testdata annotations  # module TOMLs plus golden-verified per-tool scoring TOMLs
apb-testdata fasta
```

The cache and its CSV manifests are generated and remain outside git.

## Command-line interface

The umbrella CLI is `apb` (the installed Python package is `anndata-proteomics`).
`convert` creates the common input; `annotate`, `fasta`, and `proteobench` are
independent enrichments that can be applied in any order. `validate` and `list`
maintain parsing rules.

### Convert

```bash
# Default: convert every level the file/version provides → a multi-level MuData (.h5mu).
apb convert report.tsv --params report.log.txt

# A single level → a single-level AnnData (.h5ad).
apb convert report.tsv ion     --params report.log.txt
apb convert report.tsv protein --params report.log.txt
```

The parameter file gives the **software version**, which selects the matching software-version document (e.g. DIA-NN `v1` vs `v2`); the data columns must then match one or more levels in that document. The vendor is auto-detected from the column headers — override with `--software <slug>` (the rule-folder slug, e.g. `diann`). Pass `--rule-config my_rules.json` to use an external document; add `LEVEL` to write that level as `.h5ad`, or omit it to write all matching levels as `.h5mu` (including a one-modality MuData for single-level vendors). Output defaults next to the input. `--output` accepts an extensionless basename; APB appends the suffix matching the object it produces.

### Annotate `obs` with sample metadata

```bash
apb annotate data.h5mu module_settings.toml     # writes data.annotated.h5mu
```

Joins a ProteoBench module's top-level `[[samples]]` records onto `obs` (the run axis,
shared across MuData modalities). APB also reads annotation-only TOML using
`[[obs.samples]]`, plus CSV and TSV sample tables. The `raw_file` field is matched against
`obs_names` by default; every other field becomes an `obs` column. Annotation tables are
loaded directly and are not Pydantic models.

```toml
[[samples]]
raw_file = "LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_01"
sample_name = "Condition_A_Sample_Alpha_01"
condition = "A"

[[samples]]
raw_file = "LFQ_Orbitrap_AIF_Condition_B_Sample_Alpha_01"
sample_name = "Condition_B_Sample_Alpha_01"
condition = "B"
```

### Annotate and validate against FASTA

```bash
apb fasta data.h5mu proteome.fasta              # writes data.fasta.h5mu
apb fasta data.h5mu human.fasta crap.fasta      # multiple FASTA files
apb fasta data.h5mu proteome.fasta --no-validate # protein annotation only
```

Protein layers receive a prolfquapp-style annotation in `varm['fasta']`
(`fasta.id`, `fasta.header`, `protein_length`, `nr_peptides`, `gene_name`, and
decoy/contaminant classification). The join uses the leading protein-group
accession; `nr_peptides` is the **theoretical** in-silico digest count using the
stored search enzyme (override with `--cleavage`, `--min-length`,
`--max-length`).

Validation is enabled by default for every ion, fragment, peptidoform, or
peptide modality. `varm['fasta_validation']` reports whether the peptide occurs,
the total number of occurrence sites, the number and IDs of distinct matching
proteins, whether the reported leading protein exists, and whether the peptide
occurs in that leading protein. Unmatched features and their quantifications are
never removed. In MuData, representable peptide-feature → protein-feature edges
are added to the MuLink-compatible `varp['feature_mapping']` sparse matrix.

Decoy and contaminant patterns are inferred from raw FASTA IDs and persisted in
`uns['anndata_proteomics']['fasta_config']`; `--decoy-pattern` and
`--contaminant-pattern` override inference. Classification never filters
quantified rows or FASTA records.

### Compute ProteoBench scores

```bash
apb proteobench data.h5mu module_settings.toml parse_settings_diann.toml
# writes data.proteobench.h5mu
```

The module TOML is the experiment-design contract: it supplies the scoring
level, run-to-condition samples, species mapping, and expected ratios. The
per-tool ProteoBench TOML resolves vendor source columns through the conversion
rule stored in the object. Neither `apb annotate` nor `apb fasta` is required;
running either before or after scoring preserves the other enrichments.

For a standalone AnnData, the module level must match the object. For MuData,
only that modality is scored. Feature-aligned means, standard deviations, CVs,
observation counts, species assignments, epsilon, and precision values are
stored in `varm['proteobench']`. Compatible JSON field names and threshold
layout are retained in `uns['proteobench']['scores']`; resolved source roles and
compact protein-mapper provenance live beside it in `column_roles` and
`protein_mapping`.

The managed test-data registry currently advertises golden-verified settings
for DIA-NN AIF, Astral, diaPASEF, ZenoTOF, and single-cell ion modules, plus the
WOMBAT DDA peptidoform module. The library validates any supplied ProteoBench
TOMLs, but combinations outside that list are not claimed as golden-compatible.

### Inspect / maintain rules

```bash
apb list                      # list packaged parsing rules
apb validate                  # validate all packaged rules (or: apb validate my_rule.json)
apb export-schema             # regenerate source-document and effective-rule JSON Schemas
```

## Adding a new conversion (JSON)

A parsing-rule document is `rules.json` under `src/anndata_proteomics/parsing_rules/<vendor>/` (version-specific documents go in a `v1/`, `v2/`, … subfolder). One document covers one existing software-version group and contains a shared `base` plus a `levels` object. The full schema is in [docs/json_schema.md](docs/json_schema.md); validate your draft with `apb validate path/to/rules.json`.

Every document opens with `schema_version`, `file_version`, `software_name`, and `software_version` (a regex matched against the version from the parameter file). The keys under `levels` define the quantification levels; level fragments do not repeat that field.

**Minimal long rule** (one row per run × feature):

```json
{
  "schema_version": "0.1",
  "file_version": "1",
  "software_name": "MyTool",
  "software_version": "^1\\..*",
  "base": {
    "input_shape": "long",
    "axis": {"obs_keys": ["Run"], "duplicates": {"mode": "error"}},
    "columns": {"obs": {"select": {"Run": "R.FileName"}}}
  },
  "levels": {
    "ion": {
      "axis": {"var_keys": ["Precursor_Id"], "x_layer": "Intensity"},
      "columns": {"var": {"select": {"Precursor_Id": "PEP.Id"}}},
      "layers": [{"name": "Intensity", "source": "PEP.Quantity"}]
    }
  }
}
```

**Minimal wide rule** (samples as columns, one row per feature):

```json
{
  "schema_version": "0.1",
  "file_version": "1",
  "software_name": "MyTool",
  "software_version": "^22\\..*",
  "base": {
    "input_shape": "wide",
    "axis": {"obs_keys": ["sample"], "duplicates": {"mode": "error"}},
    "columns": {"obs": {"select": {"sample": "<sample>"}}}
  },
  "levels": {
    "ion": {
      "axis": {"var_keys": ["Precursor_Id"], "x_layer": "Intensity"},
      "columns": {"var": {"select": {"Precursor_Id": "Peptide", "Charge": "Charge"}}},
      "layers": [{"name": "Intensity", "source": "^(?P<sample>.+) Intensity$"}]
    }
  }
}
```

Two further objects — `columns.var.compute` (ProForma derivation) and `modifications` (vendor mod-token mapping) — appear on most shipped rules and get their own section just below. Layers can also be factor-encoded (`encoding_mode = "factor"` with a `categories` map; e.g. FragPipe's `Match Type`).

`apb` discovers `rules.json` automatically — no registry edits. The source and every effective level are validated by Pydantic.

**Base and levels.** Shared objects live in the document's `base`; level-specific axis, columns, layers, and fragment behavior live under `levels.<level>`. APB deep-merges each level over the base before effective-rule validation. There are no inheritance paths or external base files.

## ProForma sequences & modifications

Most shipped rules standardise sequence identifiers and peptide modifications so that features are comparable across vendors (protein-level rules, which have no peptide sequence, do not). Two JSON fields do this.

**Computed columns** — `columns.var.compute` derives standard `var` columns from selected ones. The `how` recipes are: `proforma_sequence` (vendor modified sequence → [ProForma 2.0](https://github.com/HUPO-PSI/ProForma)), `stripped_sequence` (sequence with modifications removed), and `proforma_ion` (peptidoform + charge → a precursor-ion id). These become the `var_keys` / `x_layer` targets.

```json
"compute": [
  {"name": "ProForma_peptidoform", "from": ["Modified_Sequence"], "how": "proforma_sequence"},
  {"name": "ProForma_peptide", "from": ["Modified_Sequence"], "how": "stripped_sequence"},
  {"name": "ProForma_ion", "from": ["ProForma_peptidoform", "Precursor_Charge"], "how": "proforma_ion"}
]
```

**Modification mapping** — `modifications` turns a vendor's modified-sequence column into a normalised ProForma string by mapping each vendor mod token to a UNIMOD accession. `parser = "token_regex"` extracts tokens with `token_pattern`; each `map` item maps one token; `unknown_policy` decides what happens to unmapped tokens (`preserve` keeps them verbatim).

```json
"modifications": {
  "source_column": "Modified.Sequence",
  "parser": "token_regex",
  "token_pattern": "\\(([^()]*)\\)",
  "token_position": "after_residue",
  "unknown_policy": "preserve",
  "output_column": "proforma_sequence",
  "map": [{"token": "UniMod:35", "accession": "UNIMOD:35"}]
}
```

Vendors that encode modifications as **mass deltas** rather than UniMod names use the same mechanism with a different pattern and map keys — e.g. FragPipe writes `M[15.9949]`, so `token_pattern = "\\[([^\\]]+)\\]"` and `token = "15.9949"` → `accession = "UNIMOD:35"`.

## Scope

APB is a pure library plus the `apb` CLI. It ships no GUI.

## Limitations & next steps

- Conversion coverage is one software-version document per existing version group, with the levels listed above; other versions may parse but are untested.
- MuLink relationships can only target protein features already present in a
  MuData. Full FASTA match IDs remain in `varm['fasta_validation']` even when no
  corresponding protein feature exists.
- HYE ion and peptidoform scoring is implemented; PYE/plasma, de-novo, and
  entrapment metrics remain future extensions.

## Documentation

Browse the generated documentation site:

- Published site: [anndata-omics-bridge.github.io/anndata-proteomics-bridge](https://anndata-omics-bridge.github.io/anndata-proteomics-bridge/)
- Local build: run `make docs`, then open `public/index.html`

Source pages:

- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) — module map, public API, data flow
- [docs/json_schema.md](docs/json_schema.md) — JSON parsing-rule schema spec
- [docs/parameter_parsers.md](docs/parameter_parsers.md) — vendor parameter-file parsers
- [docs/parsing_architecture.md](docs/parsing_architecture.md) — subsystem UML / diagrams
