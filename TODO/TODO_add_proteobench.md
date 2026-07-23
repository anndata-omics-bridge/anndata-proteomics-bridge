# TODO: add ProteoBench scoring to APB and APB Studio

> Approved implementation tracker. The HYE-first implementation is complete;
> PYE/plasma and additional vendor compatibility remain follow-up work.

## Goal

Move the reusable quantitative ProteoBench calculations into APB under
`src/anndata_proteomics/proteobench/`, compute them directly from converted
AnnData/MuData objects, and expose the operation as an APB Studio pipeline stage.

The first implementation covers quantitative LFQ HYE-style modules at the `ion`
and `peptidoform` levels. It must:

- compute the ProteoBench intermediate feature statistics (condition means and
  standard deviations, CVs, observation counts, fold changes, species, epsilon,
  and precision epsilon);
- store the feature-aligned derived columns in `varm`;
- compute the thresholded ProteoBench score dictionary;
- retain ProteoBench's JSON field names and nested `results` format;
- reproduce the checked-in/downloaded ProteoBench result JSONs and
  `result_performance.csv` files for the same inputs;
- use a matrix-native implementation so APB does not rebuild ProteoBench's large
  long-form table merely to group and pivot it again.

PYE/plasma metrics are the next quantitative variant because they reuse the same
intermediate table. De-novo and entrapment scoring are separate future work and
must not be pulled into this migration.

## Implementation result (2026-07-22)

- `apb proteobench` scores a converted or independently enriched H5AD/H5MU from
  the required module and per-tool TOMLs.
- Feature intermediates are stored at `varm["proteobench"]`; the compatible
  score document is stored at `uns["proteobench"]["scores"]`.
- `uns["proteobench"]["protein_mapping"]` records the exact module
  `species_mapper`, accession-mapper asset hash/size, and matched/unmatched token
  counts without embedding the 38k-row mapper in every result.
- The DIA-NN Astral golden test compares the full legacy intermediate and every
  score threshold. A wider audit reproduced exact threshold feature counts for
  DIA-NN AIF, Astral, diaPASEF, ZenoTOF, and single-cell, plus WOMBAT
  peptidoform. Only these combinations are advertised by the managed settings
  registry.
- FragPipe Q Exactive/Astral and MaxQuant Astral are now advertised as
  golden-compatible. APB conversion produces the complete protein value that
  ProteoBench passes to its shared accession mapper:
  - FragPipe conversion previously dropped `Mapped Proteins` instead of joining it to
    `Protein`;
  - MaxQuant conversion previously left missing `Proteins` values instead of filling
    them from `Leading proteins`.
  The declarative conversion fix plus MaxQuant duplicate aggregation reproduces
  every golden intermediate value and score metric.
- Spectronaut is not part of this protein-conversion gap. Its short accessions
  map with zero unmatched tokens, and the compatible Astral fixtures checked
  during the diagnosis reproduce exact threshold counts. Any other
  Spectronaut discrepancy must be isolated independently before changing its
  conversion rule.
- On DIA-NN Astral `269bb831`, plain H5AD load + APB scoring took 8.94 s versus
  11.32 s for the pinned ProteoBench raw parse + intermediate + datapoint path.
  FASTA-enriched input produced identical scores and no scoring speedup
  (8.71 s scoring versus 8.72 s plain; load 0.28 s versus 0.22 s), so FASTA
  remains an independent optional enrichment.

## Ownership and dependency direction

- APB owns the reusable scoring implementation, storage contract, validation,
  and CLI.
- APB Studio orchestrates the APB CLI and renders stored results. It must not
  calculate proteomics metrics.
- APB must not depend on the `proteobench` Python package. Migrate the relevant
  formulas and test their compatibility; do not copy the ProteoBench UI,
  submission, GitHub, plotting, or parser object hierarchy.
- Keep `src/anndata_proteomics/proteobench/__init__.py` empty, matching the APB
  package convention.

Relevant ProteoBench source areas are:

- `proteobench/score/quantscoresHYE.py` — intermediate statistics and epsilon;
- `proteobench/datapoint/quant_datapoint.py` — accuracy, precision, CV, variance,
  feature-count, and ROC-AUC metrics;
- `proteobench/score/quantscoresPYE.py` and `QuantDatapointPYE` — later plasma
  extension;
- `proteobench/modules/quant/benchmarking.py` — orchestration semantics only.

Do not migrate `ScoreBase`, `DatapointBase`, timestamp-generated IDs, repository
append/deduplication, plotting filters, or web-interface code. If source is copied
rather than independently adapted to APB's matrix representation, retain the
required Apache-2.0 attribution and record the source revision.

### Source baseline

Implementation pins ProteoBench v0.17.0 commit
`fc95e712ca0466485814d3895087a048cfc0d2b0`; provenance is stored with each
scored object and the redistributed mapper includes its notice and Apache-2.0
license.

The ProteoBench checkout was fetched while writing this plan. The current
checkout (`90773d2d`) is ahead 3 / behind 67 relative to
`origin/intermediate_format_interface`; current `origin/main` is `e1b132a1`
(v0.17.2). Before implementation, deliberately select and record one source
revision. Golden JSON compatibility remains authoritative, including each
fixture's recorded `proteobench_version`.

## Pipeline independence

`convert` is the only prerequisite. The three enrichment operations are peers:

```text
convert
  |-- annotate      -> owns sample annotations in obs
  |-- fasta         -> owns FASTA-derived varm/varp data
  `-- proteobench   -> owns ProteoBench varm/uns data
```

Each operation must preserve fields written by the others and must be runnable
before or after either of them. A command applied to an accumulating artifact
therefore produces the same logical result for every ordering of `annotate`,
`fasta`, and `proteobench`.

Independent Studio DAG children of `convert` produce separate artifacts; they
do not implicitly merge. If Studio needs one cumulative artifact, model that as
an explicit composition target which chains the selected enrichments in any
order. Do not create a hidden merge or make one enrichment a prerequisite of
another.

## Existing APB inputs

The current converted `.h5ad` / `.h5mu` outputs plus the ProteoBench module TOML
and per-tool parsing TOML provide the complete scoring input:

- the designated quantitative matrix is `X`;
- runs are rows and quantified features are columns;
- converted `obs` retains the source run identifier that can be matched to each
  `[[samples]].raw_file` entry in the module TOML;
- the quantification level and software live under
  `uns["anndata_proteomics"]`;
- reported protein identifiers are retained in `var`, although their current
  column names differ by vendor;
- the downloaded ProteoBench `module_settings.toml` contains
  the `[[samples]]` run-to-condition design, `species_mapper`,
  `species_expected_ratio`, `min_count_multispec`, and quantitative level.
- the ProteoBench per-tool TOML contains the raw-vendor-to-ProteoBench mapper and
  contaminant/decoy interpretation.

**Neither annotation nor FASTA is a scoring input.** ProteoBench gets the sample
design from the required module TOML and scores species from reported protein
identifiers; it never loads a FASTA. `apb proteobench` must therefore work
directly after `convert`, must not require `obs["condition"]`, and must not inspect
`varm["fasta"]` or `varm["fasta_validation"]`. It must produce the same scores
when annotation and/or FASTA enrichment already exists or is added later.

The canonical local test-data cache also has both compatibility oracles:

- `<fixture>/result_performance.csv` for the full intermediate table;
- `<repo>-main/<intermediate_hash>.json` for the score payload.

## Required input contract

Do not add vendor-name branches or ProteoBench-specific roles to APB conversion
rules. Conversion remains consumer-neutral.

### 1. Load and resolve the per-tool ProteoBench settings

Require the ProteoBench per-tool TOML as a separate CLI/library input. Its
`[mapper]` identifies the exact raw vendor column standardized as `Proteins` and
its `[general]` block defines contaminant/decoy handling.

The converted object already stores its complete effective APB rule as JSON at
`uns["anndata_proteomics"]["rule_json"]`. Resolve each raw source name from the
per-tool mapper through `columns.var.select`, `columns.obs.select`, and the
declared layers to its exact APB location. Do not use sanitization heuristics or
software-name switches. Write the resolved mapping to:

```python
target.uns["proteobench"]["column_roles"]
```

For DIA-NN this resolves `Protein.Ids -> Protein_Ids`; it must not substitute
`Protein_Group`. Species flags are computed by applying the module TOML's
`species_mapper` substrings to that reported-protein-ID field. Do not replace
reported assignments with theoretical FASTA matches.

Audit every supported APB vendor/level against the actual per-tool mapper. If a
referenced raw column is not retained by conversion, retain that exact source
column in the relevant APB rule before declaring the combination scoreable.
Do not add speculative columns: the currently fetched MaxQuant per-tool TOML,
for example, does not map `Reverse` or `Potential contaminant`.

## Follow-up: produce ProteoBench's complete protein value during conversion

### Requirements

The converted AnnData must contain the complete tool-reported protein value
that ProteoBench standardizes as `Proteins` before applying `mapper.csv`.

This follow-up must:

- fix the problem at conversion time so every downstream consumer sees the
  completed protein assignment;
- preserve the existing tool-specific APB output column names (`Protein` for
  FragPipe and `Proteins` for MaxQuant);
- preserve the ProteoBench per-tool TOMLs unchanged:
  - FragPipe continues to declare `Protein = "Proteins"`;
  - MaxQuant continues to declare `Proteins = "Proteins"`;
- preserve the generic ProteoBench role resolver and shared scorer unchanged;
- preserve raw leading/mapped-protein source columns when selected;
- keep conversion consumer-neutral: the new operations are generic parsing-rule
  computations, not software-name branches or ProteoBench-specific roles;
- remain independent of annotation and FASTA.

The shared accession mapper continues to use the second `mapper.csv` column
(`gene_name`) as its lookup key. It maps each short reported token to the
corresponding description, leaves already-complete or unknown tokens unchanged,
and then applies the module TOML's species substrings to the combined protein
string.

### Verified root cause

| Fixture | Current threshold-1 difference | Conversion information lost | In-memory restoration |
| --- | ---: | --- | ---: |
| FragPipe DDA Q Exactive `b9f217a2` | +557 | `Mapped Proteins` omitted | 0 |
| FragPipe DDA Astral `54db0c0b` | +476 | `Mapped Proteins` omitted | 0 |
| MaxQuant DDA Astral `e5a709af` | -56 | missing `Proteins` not filled | 0 |

The two FragPipe raw files contain 5,610 and 4,078 rows, respectively, with
non-empty `Mapped Proteins`. Some assignments cross species; for example a
YEAST primary protein can have a HUMAN mapped protein. ProteoBench joins both,
marks the feature multi-species, and removes it. APB currently retains only the
primary protein and therefore keeps the feature.

The MaxQuant raw file contains 188 rows with missing `Proteins`, collapsing to
56 converted ion features. ProteoBench fills those values from
`Leading proteins` before mapping. Applying the same fallback to the converted
feature metadata removes the complete `-56` threshold-1 discrepancy and makes
all thresholds exact.

Full score comparison during implementation exposed a second MaxQuant
conversion mismatch that feature counts do not reveal: MaxQuant evidence can
contain repeated rows for the same run and ion. ProteoBench sums those
intensities before computing condition statistics; the APB MaxQuant rule's
previous duplicate policy retained one value. For example, one canonical ion
has raw intensities `3,321,800` and `328,300` in one run, and the golden
intermediate contains their sum, `3,650,100`. The MaxQuant rule must therefore
use APB's existing `duplicates.mode = "aggregate"` policy.

Spectronaut is the counter-check: its converted `PG_ProteinGroups` contains
short accessions, but the shared mapper resolves every observed token. The
protein mapper is therefore not the cause of a remaining Spectronaut score
difference.

### Design

Extend the existing parsing-rule `columns.var.compute` mechanism with two
generic string operations:

1. `coalesce`
   - takes two or more selected columns in priority order;
   - returns the first non-null value;
   - does not treat an empty string as null, matching ProteoBench's Pandas
     `fillna` behavior.
2. `join_nonempty`
   - takes two or more selected string columns;
   - skips null/empty inputs;
   - joins the remaining values using a required rule-declared separator;
   - retains the order declared in `from`.

Allow a computed column to replace a selected column with the same name. This
lets the conversion rule keep the raw-source mapping used by the per-tool TOML
while materializing the completed value in place:

```json
{
  "select": {
    "Proteins": "Proteins",
    "Leading_Proteins": "Leading proteins",
    "Leading_Razor_Protein": "Leading razor protein"
  },
  "compute": [
    {
      "name": "Proteins",
      "from": ["Proteins", "Leading_Proteins"],
      "how": "coalesce"
    }
  ]
}
```

FragPipe uses the same pattern:

```json
{
  "select": {
    "Protein": "Protein",
    "Mapped_Proteins": "Mapped Proteins"
  },
  "compute": [
    {
      "name": "Protein",
      "from": ["Protein", "Mapped_Proteins"],
      "how": "join_nonempty",
      "separator": ","
    }
  ]
}
```

The materialization order remains select first, compute second, so the compute
overwrites the selected raw value before axis construction. Deduplicate
`ColumnGroup.names` while preserving declaration order so an in-place computed
column appears exactly once in `var`.

Do not introduce a new consensus protein column. Do not change the module TOML,
per-tool TOML, `resolve_roles()`, mapper, species matching, or scoring formulas.

### Implementation plan

- [x] Extend `ColumnComputeMode` and `ColumnCompute` in
      `src/anndata_proteomics/rules/schema.py` with `coalesce`,
      `join_nonempty`, and the separator contract.
- [x] Update computed-column validation:
      - permit these generic output names;
      - permit intentional replacement of an already selected output;
      - validate source counts and require `separator` only for
        `join_nonempty`;
      - keep the existing fixed ProForma names and level-specific validation.
- [x] Implement vectorized `coalesce` and `join_nonempty` in
      `src/anndata_proteomics/converters/assemble.py`.
- [x] Deduplicate materialized column names so in-place computed outputs produce
      one `var` column.
- [x] Update the FragPipe ion rule to retain `Mapped Proteins` and compute the
      completed `Protein` value in place.
- [x] Update the MaxQuant ion rule to retain exact raw `Leading proteins` and
      compute the completed `Proteins` value in place. Keep
      `Leading razor protein` separately where already retained, and aggregate
      repeated evidence rows for the same run/ion.
- [x] Regenerate the committed parsing-rule JSON Schemas.
- [x] Add schema and materialization unit tests covering:
      - first value present;
      - first value null;
      - all values null;
      - empty strings versus null;
      - deterministic join order and separator;
      - in-place replacement;
      - H5AD categorical/string round-trip.
- [x] Add converter regression tests asserting the completed protein values for
      representative FragPipe and MaxQuant rows.
- [x] Extend the golden harness with FragPipe Q Exactive, FragPipe Astral, and
      MaxQuant Astral:
      - exact included feature identities and species flags;
      - exact `nr_feature` for thresholds 1 through 6;
      - all score metrics within the established floating tolerance;
      - unchanged resolved ProteoBench `column_roles`.
- [x] Advertise those tool/module pairs only after the full golden comparisons
      pass.
- [x] Audit any remaining Spectronaut discrepancy separately; do not modify its
      protein columns unless a failing protein-assignment comparison proves a
      conversion defect.

### Acceptance

- Converted FragPipe `var["Protein"]` contains primary and mapped proteins in
  the same logical order used by ProteoBench.
- Converted MaxQuant `var["Proteins"]` contains the ProteoBench fallback value
  whenever raw `Proteins` is null.
- The per-tool TOMLs, `resolve_roles()`, mapper, and scorer require no changes.
- The three verified fixtures reproduce every golden threshold count and score
  metric.
- DIA-NN, WOMBAT, and currently compatible Spectronaut results remain
  unchanged.

### 2. Load the scoring subset of module settings as a typed APB model

Add a small, strict model for:

- the `[[samples]]` run identifier and condition mapping;
- species flag to species-name mapping;
- expected A/B ratio and optional color by species;
- `min_count_multispec`;
- quantitative level;
- default cutoff and maximum observation threshold when they differ from the
  HYE defaults.

Read these fields from the same ProteoBench module TOML already used for sample
annotation. The scorer must align converted observation identifiers with
`[[samples]].raw_file` itself and derive its condition masks locally. Do not
duplicate module values in APB Studio, introduce a second configuration file,
or require `apb annotate` to copy those values into `obs` first. If compatible
sample annotations are already present, validate them against the module TOML;
do not use them as an alternative source of truth. Keep ordinary
`annotation.loader` focused on writing the sample table; the ProteoBench module
owns the scoring-specific reader.

## APB storage contract

Resolve exactly one target AnnData: the input itself for a standalone object, or
the modality named by the module level for MuData.

### Intermediate values

Store the derived, feature-aligned table at:

```python
target.varm["proteobench"]
```

Its index must equal `target.var_names` in the same order. Use ProteoBench column
names where they already exist:

- `log_Intensity_mean_<condition>` and `log_Intensity_std_<condition>`;
- `Intensity_mean_<condition>` and `Intensity_std_<condition>`;
- `CV_<condition>`;
- `log2_A_vs_B`;
- `nr_observed`;
- one boolean column per configured species plus `unique` and `species`;
- `log2_expectedRatio` and `epsilon`;
- `log2_empirical_median`, `log2_empirical_mean`,
  `epsilon_precision_median`, and `epsilon_precision_mean`.

Do not duplicate the per-run intensity matrix in `varm`; it already lives in
`X`. Provide one internal compatibility assembler that reconstructs the legacy
`result_performance.csv` column order from `X`, `obs`, `var`, and the stored
derived table. Use that assembler for regression tests and the legacy
`intermediate_hash`.

### Score JSON

Store the JSON-compatible score mapping at:

```python
target.uns["proteobench"]["scores"]
```

The decoded object must retain the ProteoBench score layout and names, notably:

```json
{
  "intermediate_hash": "...",
  "results": {
    "1": {"median_abs_epsilon_global": 0.0, "CV_median": 0.0},
    "2": {}
  },
  "median_abs_epsilon_global": 0.0,
  "mean_abs_epsilon_global": 0.0,
  "median_abs_epsilon_eq_species": 0.0,
  "mean_abs_epsilon_eq_species": 0.0,
  "median_abs_epsilon_precision_global": 0.0,
  "mean_abs_epsilon_precision_global": 0.0,
  "median_abs_epsilon_precision_eq_species": 0.0,
  "mean_abs_epsilon_precision_eq_species": 0.0,
  "nr_feature": 0
}
```

Requirements:

- threshold keys are strings, as in the repository JSONs;
- retain exact metric spelling and case (`CV_median`, `roc_auc`, etc.);
- top-level metrics are projections from the configured default cutoff;
- convert NumPy/Pandas scalars to ordinary JSON-compatible values and non-finite
  values to `None`/JSON `null`; never emit non-standard `NaN` tokens;
- use the legacy SHA-1 over the reconstructed intermediate's `to_string()` for
  `intermediate_hash` while compatibility is required;
- keep submission-only metadata (`submission_comments`, dataset URL, PR state,
  etc.) outside the core scorer. Existing APB search parameters may be merged
  into a full submission document later without changing the score block.

Record APB-side provenance and the storage schema version separately under
`uns["anndata_proteomics"]`; do not add APB-only keys inside the compatible
ProteoBench score document.

On rerun, fail clearly if either `varm["proteobench"]` or
`uns["proteobench"]["scores"]` is already present instead of silently
overwriting a previous scoring run. Preserve any other entries already present
in the `uns["proteobench"]` namespace.

## Computation design

### Matrix-native intermediate

Port the formulas, not the long-table implementation:

1. Validate the target level, unique feature names, resolved per-tool columns,
   complete one-to-one
   alignment between converted run identifiers and `[[samples]].raw_file`,
   and required A/B conditions.
2. Read `X` as the ProteoBench intensity matrix. A valid observation is finite
   and strictly greater than zero.
3. For each condition, slice rows once and compute per feature with NumPy/SciPy:
   positive count, arithmetic mean/std, mean/std of log2 intensity, and CV.
   Match Pandas' sample standard deviation (`ddof=1`) and missing-value behavior.
4. Compute `nr_observed` across all runs and `log2_A_vs_B` from condition log
   means.
5. Assign species from the reported-protein-ID strings using `species_mapper`;
   exclude decoy/contaminant and over-`min_count_multispec` features, then retain
   exactly one-species features for epsilon.
6. Compute expected-ratio epsilon and per-species empirical centers/precision.
7. Write the var-aligned derived frame once.

Support dense arrays first but keep reductions chunked so peak memory is bounded.
Use SciPy sparse reductions when `X` is sparse; do not densify the complete
matrix. Accumulate in float64 to match ProteoBench/Pandas numerical behavior.

### Aggregate scores

For every threshold `1..max_nr_observed`, reproduce:

- median/mean absolute epsilon globally, equally weighted by species, and per
  species;
- median/mean empirical log2 centers and absolute precision epsilon globally,
  equally weighted by species, and per species;
- `CV_median`, `CV_q75`, `CV_q90`, and `CV_q95` from the two conditions;
- `variance_epsilon_global` using Pandas-compatible sample variance;
- `nr_feature`;
- absolute-fold-change `roc_auc` with the unchanged species inferred from the
  expected ratio nearest 1:1.

Avoid adding scikit-learn solely for ROC-AUC. APB already depends on SciPy; use
`scipy.stats.rankdata(method="average")` and the equivalent Mann-Whitney/AUC
formula, verified against the golden JSONs including tied values.

PYE later adds its current spike-in error, spike-in count, human-plasma dynamic
range, and plasma epsilon fields without changing the HYE intermediate contract.

## Minimal APB implementation surface

Proposed files:

```text
src/anndata_proteomics/proteobench/
  __init__.py        # empty
  config.py          # typed module-settings subset
  intermediate.py    # matrix reductions + legacy-frame assembler/hash
  metrics.py         # HYE metric dictionary; later PYE extension
  pipeline.py        # resolve target, validate, store results
```

Expose one library operation from `pipeline.py`, tentatively:

```python
score_quantification(obj, module_settings, tool_settings) -> obj
```

Do not add layer/vendor override parameters in the first version. `X` and the
standardized converted APB representation are the data contract; the module
TOML is the experiment-design and scoring contract.

Add the CLI command:

```text
apb proteobench <converted-or-enriched.h5ad|h5mu> <module-settings.toml> \
  <per-tool.toml> \
  --output <scored.h5ad|h5mu>
```

The command loads the container, calls the library operation, and writes the same
container type non-destructively. Extend APB's `describe()` output to return the
stored score payload when present so consumers never need to know the HDF5
encoding detail.

## APB Studio integration

Enable this only after core golden parity is green.

1. Add a `proteobench` stage to `config/registry.yaml`:
   - depends only on `convert`;
   - uses the module TOML through a clearly named `module_settings` resource;
     reuse/alias the existing resource path rather than copying the file;
   - resolves the matching per-tool TOML as a separate `tool_settings` resource;
   - invokes
     `apb proteobench {input} {module_settings} {tool_settings} --output {output}`;
   - writes a distinct `*.proteobench.h5ad/.h5mu` artifact.
2. Add a registry field for supported quantitative branches/levels and teach
   target expansion to omit the stage for protein/fragment branches. Do not run
   an unsupported command and convert the expected skip into `FAILED`.
3. Add the matching Snakemake rule, artifact regex, failure marker, provenance
   sidecar, and DAG/status tests. The dashboard's stage and basket columns should
   continue to come from the registry.
4. Let the AnnData browser show the stored ProteoBench JSON through APB's
   `describe()` result. Studio renders only; it does not recalculate or reshape
   metrics.
5. Make annotation, FASTA, and ProteoBench independent children of `convert`.
   Each stage is blocked only by its own missing resource; the other enrichment
   targets remain runnable. Keep the current `BLOCKED` contract.
6. If a cumulative Studio artifact is required, add an explicit composition
   target after the independent stages are working. It must run the ordinary APB
   commands in a declared order and rely on their order-independent storage
   contracts, not duplicate enrichment logic or merge HDF5 internals.

## Test plan

### Unit tests in APB

- exact condition mean/std/CV and `nr_observed` on a small hand-computed matrix;
- zeros, negative values, infinities, NaNs, one replicate (`ddof=1`), and a
  missing condition;
- species assignment, no-species, multi-species, contaminant, and decoy cases;
- exact per-tool raw-column resolution through stored `rule_json` and clear
  failures for genuinely unretained mapped columns;
- generic conversion computations for null coalescing and ordered joining,
  including in-place replacement of the selected protein output;
- FragPipe primary-plus-mapped protein conversion and MaxQuant
  missing-`Proteins` fallback;
- epsilon and precision-center formulas for two and three species;
- every thresholded score key and top-level default-cutoff projection;
- ROC-AUC missing-class and tied-score cases;
- dense/sparse equivalence;
- `varm` index alignment when feature order is non-sorted;
- H5AD/H5MU round-trip of both storage keys, including the nested score mapping;
- collision guard on rerun;
- direct scoring of a converted object with no prior annotation or FASTA;
- equality after `annotate` and/or `fasta`, including all operation orders on a
  small fixture;
- clear failure for missing/duplicate/unmatched module-TOML sample rows and for
  disagreement with already-present sample annotations;
- CLI smoke tests and `describe()` exposure.

### Golden compatibility harness

Use the existing canonical cache; do not create another fixture downloader or
copy the results repositories into APB.

For a representative fixture, and then parametrically across supported
vendors/modules:

1. load the raw APB converted object (no annotation or FASTA stage required);
2. score it with that module's downloaded TOML and matching per-tool TOML;
3. compare the reconstructed legacy intermediate to
   `<fixture>/result_performance.csv`:
   - identical feature identifiers, row order, column order, and shape;
   - booleans/strings exactly equal;
   - floats with an explicitly small tolerance;
4. compare the decoded stored score payload to the corresponding
   `<intermediate_hash>.json`:
   - all `results` thresholds and metric keys;
   - top-level default-cutoff score fields;
   - exact `nr_feature`;
   - exact `intermediate_hash` where the reference Pandas rendering is
     reproducible;
5. report missing APB standardized inputs as a concrete vendor-rule gap, then fix
   that upstream before declaring the vendor supported.

Mark full-cache cases as integration tests and skip with a precise message when
the gitignored cache is absent. Keep a small committed synthetic unit fixture so
the core calculations always run in CI.

### Performance verification

On at least one large DIA-NN Astral fixture, compare the migrated implementation
with the pinned ProteoBench implementation using the same input:

- wall time;
- peak resident memory;
- intermediate and JSON parity.

Record the command and measurements in the implementation PR/`CHANGES.md`.
Performance work is accepted only after parity. The expected gain comes from
condition-wise matrix reductions and avoiding the long DataFrame, repeated
groupbys, pivot, merge, and duplicated raw-intensity columns.

### APB Studio tests

- target expansion adds scoring only to eligible MuData/ion/peptidoform
  branches;
- command rendering reuses the module and per-tool TOMLs and consumes the
  `convert`-stage output;
- missing module or per-tool settings yields `BLOCKED`, not `FAILED`, without
  blocking annotation or FASTA targets;
- annotation, FASTA, and scoring targets are independent, and scoring results
  are identical for every enrichment order;
- scoring output, log, failure marker, provenance, baskets, clean behavior, and
  run snapshot round-trip;
- stored scores appear in the selected-container detail view.

## Implementation order

- [x] Approve the storage keys, HYE-first scope, and score-only JSON boundary.
- [x] Pin the ProteoBench source revision and record attribution obligations.
- [x] Add golden tests against one DIA-NN fixture.
- [x] Implement typed module and per-tool settings plus exact `rule_json`
      resolution; retain only genuinely missing source columns in APB rules.
- [x] Implement the vectorized HYE intermediate and legacy assembler/hash.
- [x] Implement aggregate HYE scores and JSON serialization.
- [x] Add AnnData/MuData storage, CLI, provenance, summary, and round-trip tests.
- [x] Audit golden parity and advertise only verified module/vendor pairs.
- [x] Benchmark and optimize only while retaining parity.
- [x] Add the registry-driven APB Studio stage and UI rendering.
- [x] Update `docs/ARCHITECTURE.md`, README/CLI docs, and `CHANGES.md` in both
      touched repositories.
- [x] Add generic conversion-rule `coalesce` and `join_nonempty` computed
      columns with in-place output replacement.
- [x] Complete FragPipe and MaxQuant protein values during conversion and add
      their golden parity fixtures.
- [x] Advertise FragPipe and MaxQuant only after all threshold counts and score
      metrics match.
- [ ] Add PYE/plasma as a follow-up using the same intermediate contract.

## Done when

- one `apb proteobench` command turns a converted APB object plus its required
  module and per-tool TOMLs into a scored object with aligned
  `varm["proteobench"]`, resolved
  `uns["proteobench"]["column_roles"]`, compact
  `uns["proteobench"]["protein_mapping"]`, and a
  compatible JSON mapping at `uns["proteobench"]["scores"]`;
- the command has no annotation or FASTA precondition and returns the same
  scores for every ordering of optional annotation and FASTA enrichment;
- the representative golden intermediate and score JSON match, followed by all
  vendors declared supported;
- FragPipe and MaxQuant protein completion is owned by their declarative
  conversion rules, with no vendor logic added to the ProteoBench scorer;
- APB Studio runs the stage only on eligible branches and displays its stored
  scores;
- APB and APB Studio tests and Ruff checks pass;
- no metric implementation, vendor switch, second module map, or duplicate
  fixture downloader has been added to APB Studio.
