# anndata-proteomics (APB)

Convert proteomics quantification output into **AnnData / MuData** using declarative TOML parsing rules.

- **Declarative, not bespoke.** Every vendor × quantification-level is a small TOML rule shipped inside the package — adding or fixing a converter means editing a `.toml`, not writing tool-specific Python.
- **One file → a multi-level MuData.** A single vendor export is converted into a MuData whose modalities are the quantification levels it provides (`ion` / `fragment` / `peptidoform` / `protein`) on a shared run axis — or a single-level AnnData when you ask for one level.
- **Standardised content.** Peptide modifications are normalised to **ProForma**; a per-vendor parser reads the vendor **parameter file** (enzyme, FDR, tolerances, …) into one typed record under `uns['search_parameters']`.
- **Enrichable.** Join sample metadata onto `obs` (`apb annotate`) and attach FASTA-derived protein annotation — theoretical (enzyme-aware) peptide counts, gene names, protein length — onto the protein modality's `varm['fasta']` (`apb fasta`).
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
uv venv --python 3.13
source .venv/bin/activate
uv pip install -e '.[dev]'   # drop [dev] if you only need the CLI (it adds pytest + ruff)
```

## Command-line interface

The umbrella CLI is `apb` (the installed Python package is `anndata-proteomics`). Typical flow: **convert** a vendor file → **annotate** / **fasta** to enrich it → **validate** / **list** to manage rules.

### Convert

```bash
# Default: convert every level the file/version provides → a multi-level MuData (.h5mu).
apb convert report.tsv --params report.log.txt

# A single level → a single-level AnnData (.h5ad).
apb convert report.tsv ion     --params report.log.txt
apb convert report.tsv protein --params report.log.txt
```

The parameter file gives the **software version**, which selects the rule variant (e.g. DIA-NN `v1` vs `v2`); the data columns must then match that rule. The vendor is auto-detected from the column headers — override with `--software <slug>` (the rule-folder slug, e.g. `diann`). Pass `--rule-toml my_rule.toml` to override rule selection entirely (single level, version-agnostic; `--params` then optional). A vendor that exposes only one level writes a `.h5ad` even without a level argument. Output defaults to `<stem>.h5mu` (MuData) or `<stem>.h5ad` (single level) next to the input; override with `--output`.

### Annotate `obs` with sample metadata

```bash
apb annotate data.h5mu annotation.toml          # writes data.annotated.h5mu
```

Joins the records in the annotation TOML onto `obs` (the run axis, shared across MuData modalities). Each record's `key_field` is matched per `match_on` (`"index"` → `obs_names`, else an `obs` column); every other field in the record becomes an `obs` column. Example from the ProteoBench `quant_lfq_ion_DIA_AIF` module (its `module_settings.toml [[samples]]` table, translated into the APB schema — 6 runs, two conditions; trimmed here):

```toml
schema_version = "0.1"

[obs]
match_on  = "index"      # match raw_file against obs_names (the run/file identifier)
key_field = "raw_file"

[[obs.samples]]
raw_file    = "LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_01"
sample_name = "Condition_A_Sample_Alpha_01"
condition   = "A"

[[obs.samples]]
raw_file    = "LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_02"
sample_name = "Condition_A_Sample_Alpha_02"
condition   = "A"

[[obs.samples]]
raw_file    = "LFQ_Orbitrap_AIF_Condition_B_Sample_Alpha_01"
sample_name = "Condition_B_Sample_Alpha_01"
condition   = "B"

# … one record per run (B_02, B_03, …)
```

### Annotate the protein layer from FASTA

```bash
apb fasta data.h5mu proteome.fasta              # writes data.annotated.h5mu
apb fasta data.h5mu human.fasta crap.fasta      # multiple FASTA files
```

Builds a prolfquapp-style protein annotation (`fasta_id`, `fasta_header`, `protein_length`, `nr_peptides`, `gene_name`) and attaches it to the protein modality's `varm['fasta']`. The join is on the leading accession of each protein group; `nr_peptides` is the **theoretical** in-silico digest count using the enzyme from `uns['search_parameters']` (override with `--cleavage`, `--min-length`, `--max-length`). Other flags: `--match-on`, `--no-is-uniprot`, `--decoy-pattern`.

### Inspect / maintain rules

```bash
apb list                      # list packaged parsing rules
apb validate                  # validate all packaged rules (or: apb validate my_rule.toml)
apb export-schema             # regenerate parse_rule.schema.json from the pydantic models
```

## Adding a new conversion (TOML)

A parsing rule is a TOML file under `src/anndata_proteomics/parsing_rules/<vendor>/parse_<software>_<level>.toml` (version-specific rules go in a `v1/`, `v2/`, … subfolder). It declares the table shape, which columns become `obs` / `var`, and which columns become measurement `layers`. The full schema is in [docs/toml_schema.md](docs/toml_schema.md); validate your draft with `apb validate path/to/rule.toml`.

Every rule opens with the same header: `schema_version`, `software_name`, `software_version` (a regex matched against the version from the parameter file — the `v*/` folder is chosen from this), `input_shape`, and `quantification_level`. (`file_version` is the rule's own revision, independent of the software version.)

**Minimal long rule** (one row per run × feature):

```toml
schema_version = "0.1"
file_version   = "1"
software_name  = "MyTool"
software_version = "^1\\..*"        # regex matched against the parsed software version
input_shape    = "long"
quantification_level = "ion"        # ion | fragment | peptidoform | protein

[axis]
obs_keys = ["Run"]                  # column(s) identifying a run    → obs index
var_keys = ["Precursor_Id"]         # column(s) identifying a feature → var index
x_layer  = "Intensity"             # which layer becomes X

[axis.duplicates]
mode = "error"                      # what to do if (run, feature) repeats

[columns.obs.select]                # rename input columns → obs
Run = "R.FileName"

[columns.var.select]               # rename input columns → var
Precursor_Id      = "PEP.Id"
Stripped_Sequence = "PEP.StrippedSequence"

[[layers]]                          # one matrix per measurement column
name   = "Intensity"
source = "PEP.Quantity"
```

**Minimal wide rule** (samples as columns, one row per feature):

```toml
schema_version = "0.1"
file_version   = "1"
software_name  = "MyTool"
software_version = "^22\\..*"
input_shape    = "wide"
quantification_level = "ion"

[axis]
obs_keys = ["sample"]
var_keys = ["Precursor_Id"]
x_layer  = "Intensity"

[axis.duplicates]
mode = "error"

[columns.obs.select]
sample = "<sample>"                 # the <sample> token = names captured from the layer headers

[columns.var.select]
Precursor_Id = "Peptide"
Charge       = "Charge"

[[layers]]
name   = "Intensity"
source = "^(?P<sample>.+) Intensity$"   # regex over column headers; (?P<sample>) captures the run name
```

Two further blocks — `[[columns.var.compute]]` (ProForma derivation) and `[modifications]` (vendor mod-token mapping) — appear on most shipped rules and get their own section just below. Layers can also be factor-encoded (`encoding_mode = "factor"` with a `categories` map; e.g. FragPipe's `Match Type`).

`apb` discovers the file automatically — no registry edits. A test enforces that the filename's level token matches the in-TOML `quantification_level`.

**Vendor base files.** For a vendor with several level rules (DIA-NN, Spectronaut), the blocks its levels share — `[modifications]`, the run-axis keys, common scalars — live in a **base file** `<vendor>/<vendor>.toml` that every leaf inherits automatically at load time (convention, no `extends` key; see [docs/toml_schema.md](docs/toml_schema.md)). A new leaf for such a vendor declares only its level-specific content; single-format vendors stay one self-contained file.

## ProForma sequences & modifications

Most shipped rules standardise sequence identifiers and peptide modifications so that features are comparable across vendors (protein-level rules, which have no peptide sequence, do not). Two TOML blocks do this.

**Computed columns** — `[[columns.var.compute]]` derives standard `var` columns from selected ones. The `how` recipes are: `proforma_sequence` (vendor modified sequence → [ProForma 2.0](https://github.com/HUPO-PSI/ProForma)), `stripped_sequence` (sequence with modifications removed), and `proforma_ion` (peptidoform + charge → a precursor-ion id). These become the `var_keys` / `x_layer` targets.

```toml
[[columns.var.compute]]
name = "ProForma_peptidoform"
from = ["Modified_Sequence"]
how  = "proforma_sequence"

[[columns.var.compute]]
name = "ProForma_peptide"
from = ["Modified_Sequence"]
how  = "stripped_sequence"

[[columns.var.compute]]
name = "ProForma_ion"
from = ["ProForma_peptidoform", "Precursor_Charge"]
how  = "proforma_ion"
```

**Modification mapping** — `[modifications]` turns a vendor's modified-sequence column into a normalised ProForma string by mapping each vendor mod token to a UNIMOD accession. `parser = "token_regex"` extracts tokens with `token_pattern`; each `[[modifications.map]]` maps one token; `unknown_policy` decides what happens to unmapped tokens (`preserve` keeps them verbatim).

```toml
# DIA-NN: mods look like AAC(UniMod:4)DEM(UniMod:35)K
[modifications]
source_column  = "Modified.Sequence"
parser         = "token_regex"
token_pattern  = "\\(([^()]*)\\)"
token_position = "after_residue"
unknown_policy = "preserve"
output_column  = "proforma_sequence"

[[modifications.map]]
token     = "UniMod:35"     # vendor token
accession = "UNIMOD:35"     # → oxidation
```

Vendors that encode modifications as **mass deltas** rather than UniMod names use the same mechanism with a different pattern and map keys — e.g. FragPipe writes `M[15.9949]`, so `token_pattern = "\\[([^\\]]+)\\]"` and `token = "15.9949"` → `accession = "UNIMOD:35"`.

## Interactive tools

apb is a pure library + `apb` CLI — it ships **no** GUI and does not depend on marimo. All
interactive [marimo](https://marimo.io) tooling lives in the sibling **`apb_studio`** package (in
its own repo), which drives apb entirely through the `apb` CLI:

- the **test-data browser** (`apb_studio … ui/test_tool.py`, `make test-tool`) — browse the
  ProteoBench corpus, convert a dataset (shelling out to `apb convert`), and inspect the result;
- the **corpus dashboard** + Snakemake pipeline (`make ui`) for whole-corpus coverage.

apb_studio consumes apb as a sibling install and imports only apb's pure read-only helpers for
catalog/metadata — conversion itself always runs via the CLI.

## Limitations & next steps

- Conversion coverage is one rule per vendor/level at the versions listed above; other versions may parse but are untested.
- Protein-level annotation (`apb fasta`) targets a protein AnnData or the `protein` modality of a MuData; other inputs are rejected.
- Per-tool `uns['<app_name>']['column_roles']` writeback (the tool-specific view ADR) is not yet populated — only `uns['anndata_proteomics']` is written.

## Documentation

Browse the generated documentation site:

- Published site: [anndata-omics-bridge.github.io/anndata-proteomics-bridge](https://anndata-omics-bridge.github.io/anndata-proteomics-bridge/)
- Local build: run `docs/render_docs.sh`, then open [public/index.html](public/index.html)

Source pages:

- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) — module map, public API, data flow
- [docs/toml_schema.md](docs/toml_schema.md) — TOML parsing-rule schema spec
- [docs/parameter_parsers.md](docs/parameter_parsers.md) — vendor parameter-file parsers
- [docs/parsing_architecture.md](docs/parsing_architecture.md) — subsystem UML / diagrams
