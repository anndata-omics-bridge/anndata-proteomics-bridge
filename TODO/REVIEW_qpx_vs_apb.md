# Review: QPX lessons for APB

Date: 2026-06-26

## Scope

Compared the local `qpx/` checkout against current APB architecture and planning notes in `apb/`.
This is not a proposal to merge designs wholesale. QPX and APB solve adjacent but different
problems:

- QPX is a Parquet-first proteomics dataset format with SQL/query/view tooling and optional
  AnnData/MuData export.
- APB is a TOML-driven vendor-output-to-AnnData/MuData converter, with one parsing rule per
  vendor quantification level.

## Findings

### 1. Do not copy QPX's converter architecture into APB

QPX converters are code adapters. `BaseConverter` expects subclasses to load tool output into
DuckDB, transform with SQL into QPX schemas, and stream rows into writers. It also has a central
YAML mapping registry where QPX fields map to ordered candidate vendor columns.

That is useful for QPX because QPX owns a canonical Parquet schema. It is the wrong center of
gravity for APB. APB's explicit plan is the opposite: vendor-specific conversion should mostly live
in declarative TOML rules, with the Python converter remaining generic.

APB should keep the TOML rule model as the upstream source of truth. Copying QPX-style adapter
classes would reintroduce per-vendor Python logic that APB has deliberately removed.

Useful lesson: QPX's ordered candidate-column mapping is a good warning sign. APB should keep exact
vendor column references in TOML by default, but if real vendor variants force aliases, add them to
the TOML schema as explicit version-scoped alternatives, not as a hidden central YAML fallback.

### 2. QPX confirms APB's "real report-backed levels only" decision

QPX distinguishes stored data views from API views. Its docs explicitly say API views are computed
on demand by joining and aggregating primary views, while stored views are concrete files. This
matches the APB rule that a standalone AnnData level is valid only when the vendor output has real
quantitative layer columns for that level.

This supports APB's current decision not to manufacture peptide or peptidoform matrices from DIA-NN
precursor columns during parsing. Derived rollups should be a separate derivation pipeline or
computed view, not a parsing TOML.

### 3. QPX's ontology/provenance tables are the strongest transferable idea

QPX carries `ontology.parquet` as a machine-readable field-to-ontology mapping and `provenance`
metadata for processing decisions. APB already stores the rule JSON and parsed search parameters in
`uns['anndata_proteomics']`, but APB does not yet have an equivalent compact, queryable "column
roles / source columns / semantic terms" surface.

This maps directly to APB's not-yet-implemented item: per-tool
`uns['<app_name>']['column_roles']` writeback. That should be prioritized before inventing new
viewer-specific metadata, because it would make APB outputs more self-describing and would surface
automatically in the existing report/viewer path.

Concrete APB shape:

- For each AnnData modality, write a small `uns['anndata_proteomics']['column_roles']` record.
- Include layer name, original vendor column, axis placement (`obs`, `var`, `layer`), level, and
  whether the field is selected or computed.
- Later add ontology accessions only where APB already has them from `unimod_registry` or parameter
  parsing. Do not block the role table on full ontology coverage.

### 4. QPX's `pepmap` is useful prior art, but APB should not add it yet

QPX has a dedicated `pepmap.parquet` for deduplicated peptide-to-protein mappings, with uniqueness
flags and optional protein positions. This directly addresses the hard APB case already identified
in `TODO_to_mu_data.md`: peptide-to-protein is many-to-many and does not fit the simple child-column
to parent-index foreign-key pattern.

For APB, the lesson is not "add `pepmap` now". The existing APB plan is still right:

- Keep `Protein_Group`, `Protein_Ids`, `Protein_Names`, and `Genes` as TOML-defined `.var`
  metadata first.
- Do not put link tables in global `mdata.var`.
- Add a derived relation table only after real use cases need peptide/protein navigation or
  uniqueness analysis.

If that need appears, QPX's `pepmap` schema is the best starting point for the derived view.

### 5. QPX's MuData export is less precise than APB's level model

QPX builds MuData modalities named `precursors`, `proteins`, `expression`, and `differential`, and
stores a sparse global `mdata.varp["feature_mapping"]` between precursors and proteins. This is
reasonable for exporting from QPX's Parquet model, but it is weaker than APB's TOML-defined
level-specific `.var` links.

APB should not replace its current level model with QPX's generic feature mapping. APB already knows
the quantification level and computed identifiers (`ProForma_ion`, `ProForma_peptidoform`,
`ProForma_peptide`, `Protein_Group`) at parse time. Those should remain in each modality's `.var`
table as authoritative metadata.

Potentially useful later: a derived sparse adjacency view could be generated from APB `.var` links
for graph-style navigation, but it should be derived from TOML-defined columns, not primary storage.

### 6. QPX's large-data choices expose APB's scaling risk

QPX is built around DuckDB, Arrow, Parquet, batched writers, compression, S3 registration, and
partitioned datasets. APB currently pivots vendor tables into AnnData matrices directly. That is
fine for APB's current converter goal, but the fragment-level notes already show the pressure point:
full DIA-NN fragment conversion can explode rows and peak at several GB.

APB should not become QPX. But APB should borrow the engineering pattern:

- Add explicit big-file guardrails for fragment and MuData conversions.
- Add row-capped/per-run conversion modes in the UI/test tool.
- If full-scale fragment remains important, implement chunked explode/scatter rather than
  materializing the full exploded frame.
- Consider sparse output only when the matrix is actually sparse; the APB note says the DIA-NN
  fragment matrix was about 90 percent dense, so sparse is not automatically a win.

### 7. QPX's tests validate converted artifacts, not just converter execution

QPX integration tests run converters once and then validate produced Parquet schemas and plausible
values. APB already has similar E2E coverage for packaged TOMLs and MuData proof tests. The lesson
is to keep artifact-level assertions as APB grows:

- Validate stored `uns` provenance/search parameters.
- Validate `column_roles` once implemented.
- Validate MuData round-trips and prefixed `var_names`.
- Validate that no standalone level appears without a real `x_layer`.

## What APB should do next

1. Keep TOML-first conversion. Do not add QPX-style per-vendor converter adapter classes.
2. Implement the planned `column_roles` / source-column writeback in `uns['anndata_proteomics']`.
3. Keep report-backed-level discipline: no derived peptide/peptidoform matrices in parse rules.
4. Treat protein/peptide relation tables as derived views. Use QPX `pepmap` as prior art only when
   the use case becomes real.
5. Add fragment/MuData guardrails before more large real-data demos.
6. Continue artifact-level tests around `.h5ad`/`.h5mu`, especially provenance and level links.

## Local evidence checked

- `qpx/README.md`
- `qpx/docs/spec/index.md`
- `qpx/docs/spec/anndata.md`
- `qpx/docs/spec/ontology.md`
- `qpx/docs/spec/pepmap.md`
- `qpx/qpx/converters/base.py`
- `qpx/qpx/converters/column_mappings.yaml`
- `qpx/qpx/mudata.py`
- `qpx/qpx/dataset.py`
- `qpx/qpx/writers/base.py`
- `apb/README.md`
- `apb/docs/ARCHITECTURE.md`
- `apb/TODO/Archive/TODO_to_mu_data.md`
- `apb/TODO/TODO_ui_test_tool.md`
- `apb/TODO/TODO_modification_homogenization_design.md`
