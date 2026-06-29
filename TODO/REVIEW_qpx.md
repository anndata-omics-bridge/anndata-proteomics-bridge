# QPX as Prior Art for APB

**Date:** 2026-06-26
**Verified against:** `qpx` at commit `63d6382` (2026-05-06).
**Method:** read the local `qpx/` checkout (code, schemas, docs) and cross-checked the
project's external lineage. Claims below are grounded in code, not upstream marketing.

---

## TL;DR

**Do not replace APB with QPX.** They solve adjacent but different problems, and APB's core
job — *vendor output → AnnData/MuData driven by a versioned declarative rule* — is not what QPX
does.

- **QPX** is a Parquet-first, multi-vendor proteomics *dataset and query format*. You convert a
  tool's output into a directory of normalised Parquet tables, then query them with DuckDB. It
  ships ontology/provenance tables and an optional MuData export.
- **APB** is a TOML-first *converter*: the parsing contract is a versioned vendor TOML rule plus a
  parsed parameter file, and the output is AnnData/MuData with explicit quantification-level
  objects and exact vendor column roles preserved in the rule.

**What to take from QPX:** its provenance / ontology / field-role metadata design, and its
large-data engineering. Treat QPX as a possible *downstream / interop target* (APB → QPX export),
not as a drop-in replacement.

---

## 1. What QPX Is

QPX stands for **Q**uantitative **P**roteomics e**X**change (the repo title also gives "aka
quantms.io"). It is developed by the **bigbio** group — the PRIDE Team at EMBL-EBI
(Yasset Perez-Riverol et al.) — and is the format/storage layer of the broader **quantms**
ecosystem. Concretely it is:

- a **Parquet-first dataset format**: one experiment becomes a directory of normalised Parquet
  tables (`feature`, `pg`, `psm`, `sample`, `run`, plus metadata tables);
- a **DuckDB query layer** over those tables (`qpxc query sql ...`);
- a **CLI** (`qpxc`) with `convert` / `transform` / `query` / `info` / `validate` / `ontology`
  command groups;
- a set of **per-vendor Python converters** (one adapter family per tool);
- an **optional MuData/scverse export** for quantification results.

It is built squarely on community standards: **SDRF-Proteomics** for experimental design,
**ProForma 2.0** for peptidoform notation, and **PSI-MS / PRIDE controlled-vocabulary** terms for
the ontology tables. Its sibling tools are `mokume` (feature → peptide/protein quantification) and
`ibaqpy` (SDRF-driven absolute quantification). There is **no dedicated QPX/quantms.io paper** yet;
the format is documented through the surrounding ecosystem (see *References*).

**One-line contrast:** QPX normalises *many tools* into *one queryable Parquet schema*; APB lifts
*one tool's output* into *AnnData/MuData* while preserving that tool's exact semantics in a
versioned rule.

---

## 2. Converter Coverage — QPX Is Not quantms-Only

QPX ships converters for seven sources. quantms/mzTab is the reference ecosystem, but the design is
explicitly multi-vendor (one Python adapter family per tool). Coverage differs per converter:

| Converter   | PSM | Feature | PG  | Pepmap | Sample/Run |
|-------------|:---:|:-------:|:---:|:------:|:----------:|
| MaxQuant    | yes | yes     | yes | no     | if SDRF    |
| FragPipe    | yes | yes     | yes | no     | if SDRF    |
| DIA-NN      | no  | yes     | yes | no     | if SDRF    |
| Spectronaut | no  | yes     | yes | no     | if SDRF    |
| quantms     | yes | yes     | yes | no     | if SDRF    |
| mzIdentML   | yes | no      | no  | yes    | if SDRF    |
| SDRF        | no  | no      | no  | no     | yes        |

*Spectronaut was added recently (commit `8fb32b5`) and its row is current. DIA-NN and Spectronaut
produce no PSM table. mzIdentML is identification-only (PSM + peptide-map, no quantification).*

---

## 3. The QPX Data Model — What a Dataset Looks Like on Disk

A conversion writes **many Parquet files into one directory**, one file per *view*, named
`{prefix}.{view}.parquet`. It is **not** one combined Parquet file, and it is **not** one file per
AnnData layer.

A full MaxQuant run (`msms.txt` + `evidence.txt` + `proteinGroups.txt` + SDRF) writes roughly:

```text
maxquant.psm.parquet         # always (if msms.txt given)
maxquant.feature.parquet     # always (if evidence.txt given)
maxquant.pg.parquet          # always (if proteinGroups.txt given)
maxquant.sample.parquet      # only if a valid SDRF is supplied
maxquant.run.parquet         # only if a valid SDRF is supplied
maxquant.ontology.parquet    # always
maxquant.provenance.parquet  # always
maxquant.dataset.parquet     # always
```

The quantitative structures auto-detect from the inputs unless `--structures` overrides them:
`msms.txt → psm`, `evidence.txt → feature`, `proteinGroups.txt → pg`. The `ontology`,
`provenance`, and `dataset` tables are written unconditionally; `sample`/`run` are written **only
when a valid SDRF is provided** (see §5).

### The three quantitative views

- **`feature`** — peptide/precursor-level records per run. Primary key
  `[sequence, charge, run_file_name, anchor_protein]`. Carries `peptidoform`, `modifications`,
  charge, m/z, mass error, scores, RT, ion mobility, missed cleavages, the two intensity columns,
  and protein-grouping fields (`pg_accessions` as structs with `start/end/pre/post`,
  `anchor_protein`, `unique`, `pg_global_qvalue`, gene-group `gg_*`).
- **`pg`** — protein-group records per run. Primary key `[anchor_protein, run_file_name]`. Carries
  `pg_accessions` (plain string list here), `pg_names`, gene-group fields, q-values, the two
  intensity columns, peptide/feature counts, sequence coverage, molecular weight.
- **`psm`** — PSM-level records where the converter supports them.

### Intensities are nested structs, not matrix layers

Each quantitative view stores **at most two intensity-bearing columns**, each a *list of structs*
with a variable number of named values per row:

- **`intensities`** — list of `{label: string, intensity: float32}` (the primary value).
- **`additional_intensities`** — optional list of
  `{label: string, intensities: [{intensity_name: string, intensity_value: float32}, …]}`
  (secondary values).

These are per-row nested values to be unnested/pivoted later in DuckDB — **not** AnnData-style
dense/sparse matrices. (Note: the struct field is now `label`; older datasets used `channel`, and
the MuData builder auto-detects which is present.)

What MaxQuant actually stores:

| Level / type            | `intensities`                                       | `additional_intensities`                              |
|-------------------------|-----------------------------------------------------|-------------------------------------------------------|
| LFQ feature             | one value labeled `LFQ` (prefers `LFQ intensity`, falls back to raw `Intensity`) | *none* (left `None`)            |
| LFQ protein group       | raw `Intensity <run>` labeled `LFQ`                 | `lfq` (← `LFQ intensity <run>`), `ibaq` (← `iBAQ <run>`) |
| TMT / iTRAQ (both)      | one value per channel                               | `corrected_reporter_intensity` per channel            |

### What QPX does *not* persist

QPX stores its normalised views, not the raw vendor tables. Dropped on conversion: unmapped vendor
columns, the raw table layout, search-engine parameters (e.g. MaxQuant's `mqpar.xml`, which is not
even an input — see §7), and any quantitative column not folded into the two intensity structs.
Exact vendor column names survive only inside the ontology field-provenance entries. **There is no
raw-table sidecar**, so a QPX dataset is not a lossless round-trip back to the original tool output.

---

## 4. How Conversion Works (engine internals)

- **DuckDB is the engine, end to end.** Vendor files are loaded with DuckDB's
  `read_csv_auto(..., delim='\t')` into named in-memory tables (e.g. MaxQuant → `evidence`, `msms`,
  `protein_groups`). Polars is **not** in the loading path — it appears only as an optional
  *output* materialiser (`to_polars()` on query results).
- **Conversion is mostly Python row transformation.** Adapters stream DuckDB Arrow batches
  (`fetch_record_batch`) into pandas and do the semantic work (ProForma parsing, decoy flags, m/z,
  intensity assembly) row by row in Python. The SQL at this stage is a trivial `SELECT *`.
- **Large mzTab files are split.** quantms/mzTab files ≥ 500 MB (and not gzipped) are split into
  temporary per-section files and loaded with DuckDB's native CSV reader into tables `metadata`,
  `proteins`, `peptides`, `psms`, and optional `msstats`. Smaller or gzipped files use an in-memory
  loader.
- **SQL is inline, but now built through a hardening helper.** There are no separate `.sql` files,
  but SQL is no longer raw f-strings: a `qpx/core/sql.py` module provides `sql_build`,
  `validate_identifier`, `validate_table`, and `escape_path` to allowlist identifiers and prevent
  SQL injection (DuckDB cannot parameterise identifiers; this satisfies Bandit B608 without
  `# nosec`). Values are bound with `?`/`$1` placeholders where possible.
- **Where SQL is genuinely central:** the quantms **LFQ** feature path is SQL-first — it builds
  DuckDB lookup tables (`_psm_lookup`, `_protein_qvalues`, `_protein_genes`, `_proforma_lookup`)
  and `LEFT JOIN`s them to the `msstats` table in one query, then assembles structs in Python. The
  MuData export is also SQL-heavy: it unnests the nested intensity structs back into matrices.
  (Isobaric quantms data instead uses a pandas `groupby` aggregation.)

---

## 5. The SDRF Dependency — "No SDRF, No Read?"

**No.** A missing SDRF does **not** stop QPX from reading or converting a vendor file. The
quantitative views (`psm` / `feature` / `pg`) and the `ontology` / `provenance` / `dataset` tables
are produced regardless. The **only** thing you lose is the `sample.parquet` and `run.parquet`
tables, which are produced *exclusively* from the SDRF. No converter aborts when the SDRF is absent
— they all degrade gracefully (and a malformed SDRF is caught, logged, and its partial files
deleted).

Two important structural facts make this clear:

1. **Run names come from the vendor file, never from the SDRF.** `run_file_name` is derived from
   the tool's own output (MaxQuant `Raw file`, DIA-NN `Run`, Spectronaut `R.FileName`,
   FragPipe experiment/ion column, quantms MSstats `Reference` / mzTab `ms_run`).
2. **`feature` and `pg` carry no sample/condition/replicate/fraction/instrument columns at all.**
   Those fields never live in the quantitative views — they live only in `sample.parquet` /
   `run.parquet`, joined back by `run_file_name`. So "without SDRF, which feature/pg fields go
   empty?" → effectively none; you simply don't get the two design tables.

Per-converter behaviour:

| Converter   | `--sdrf-file` | Without SDRF |
|-------------|---------------|--------------|
| MaxQuant    | optional      | quant + metadata tables produced; `sample`/`run` skipped |
| FragPipe    | optional (and unused by the feature adapter) | quant + metadata tables produced; `sample`/`run` skipped |
| Spectronaut | optional      | quant + metadata tables produced; `sample`/`run` skipped |
| DIA-NN      | **CLI-required** (`required=True`) but algorithmically soft — only used to read the enzyme name for missed-cleavage counts | if bypassed: feature/pg still produced; `sample`/`run` skipped; missed cleavages → `None` |
| quantms     | **CLI-required**; used for `sample`/`run` + enzyme | the real hard gate for feature/pg is **`--msstats-file`**, *not* SDRF — see below |

**quantms is the special case.** Its hard requirement for the quantitative views is the **MSstats
input file**, not the SDRF:

- Feature is produced only `if FEATURE in structures and self.msstats_file` — no MSstats, no
  `feature.parquet`.
- PG is generated *from* `feature.parquet` and **raises** if it is missing. So no MSstats → no
  feature → PG aborts.
- With an SDRF but no MSstats file, quantms yields **PSM only**.

**Implication for APB.** This is a real model difference. In APB, the AnnData `.obs` *is* the
sample/design annotation, attached per modality. In QPX the design lives in detached sidecar tables
keyed by run name; even with an SDRF, the sample linkage is *not* baked into the quantitative rows.
The MuData export (§6) reconstructs `obs` by joining `run.parquet` back in — so the richness of a
QPX-derived MuData's `obs` is exactly the richness of the SDRF you fed in. No SDRF → `obs` is little
more than the run name parsed from the vendor file.

---

## 6. MuData / scverse Export

QPX can build a `MuData` object from a dataset directory — this is its bridge to the scverse
ecosystem, and it is **quantification-only** (it does not round-trip the whole experiment).

**API gotcha:** the README advertises `ds.to_mudata()`, **but that method does not exist** in the
code. The working entry point is `qpx.mudata.build_mudata(dataset, intensity_label=None,
modalities=None)`. The README snippet would raise `AttributeError`.

```python
from qpx import Dataset
from qpx.mudata import build_mudata

ds = Dataset("qpx_output")
mdata = build_mudata(ds)          # NOT ds.to_mudata()
mdata.write("maxquant.h5mu")
```

What the builder produces:

- **`precursors`** AnnData ← `feature.parquet`; **`proteins`** AnnData ← `pg.parquet`.
- Optional **`expression`** ← a `*.pe.h5ad` / `*.pe.zarr` file if present; optional
  **`differential`** ← a `*.de.h5ad` / `*.de.zarr` file if present.
- **`X`** holds a *single* chosen intensity label (auto-detected as the first label, or passed
  explicitly). The QPX `additional_intensities` are **not** materialised as AnnData layers. The
  *only* modality that uses `.layers` is `differential` (it puts `log2FC` in `X` and adds
  `pvals_adj` / `scores` / `pvals` / `se` as layers).
- `obs` on both quant modalities is enriched from `run.parquet` (run/sample accession, replicate,
  fraction, instrument); `uns` is populated from `dataset.parquet`; protein `var` gets a
  `gene_name`; object columns are NaN-sanitised so h5py can write them.
- When both quant modalities are present, the builder adds
  `mdata.varp["feature_mapping"]` — a **symmetric boolean (N×N) adjacency** over the combined
  variable index linking precursors to their proteins (not a rectangular precursor×protein block).

---

## 7. MSstats Integration

MSstats enters QPX **only through the quantms converter** — DIA-NN, MaxQuant, FragPipe, and
Spectronaut never touch it. And QPX **does not run MSstats**; it *consumes* an MSstats input table
produced upstream (typically by the quantms / OpenMS workflow).

In the quantms path:

1. mzTab → DuckDB tables `metadata`, `proteins`, `peptides`, `psms`.
2. MSstats file → DuckDB table `msstats`.
3. SDRF → `sample` / `run`.
4. **PSM** output comes from the mzTab PSM section.
5. **Feature** output comes primarily from MSstats rows, enriched by mzTab PSM/protein metadata.
6. **PG** output is generated from the already-written `feature.parquet` plus mzTab protein
   metadata.

Columns the quantms feature path consumes from MSstats (exact names, resolved via the column
mapping): `PeptideSequence`/`peptidoform`, `ProteinName`, `Reference`/run, `Charge`, `Intensity`,
`Channel` (TMT/iTRAQ), `RetentionTime`. LFQ uses the SQL-first join path (§4); isobaric data is
grouped by feature and aggregated per channel.

**DIA-NN is a fully separate path:** it reads `report.tsv` *or* `report.parquet` (plus an optional
`pg_matrix.tsv` to enable PG) and does **not** use MSstats. So if an upstream quantms workflow
produced an MSstats table *from* DIA-NN output, QPX only ever sees the post-quantms MSstats table —
it does not model or preserve the original DIA-NN fragment-level path the way an APB vendor rule
would.

---

## 8. Fidelity Caveats

Things to know before treating a QPX export as faithful to the original output:

- **No lossless round-trip.** Unmapped columns, raw layout, and search parameters are dropped; there
  is no raw-table sidecar (§3).
- **Search parameters are not captured for MaxQuant.** `mqpar.xml` is not an input, so enzyme,
  tolerances, full modification config, and version are absent — `software_version` is written as
  `None`.
- **Some MaxQuant fields are duplicated, not independently sourced.** In the PG path `gg_names` is
  set equal to `gg_accessions`, and `pg_qvalue` is set equal to `global_qvalue`.
- **LFQ feature rows have no `additional_intensities`** (only PG-LFQ and TMT rows populate it).
- **`--standardized-intensities` is accepted but inert** for the MaxQuant PG path (a latent dead
  parameter).
- **README drift.** `ds.to_mudata()` (§6) is documented but unimplemented.

---

## 9. What This Means for APB

| QPX is stronger as…                          | APB is stronger at…                                              |
|----------------------------------------------|------------------------------------------------------------------|
| a Parquet query / interchange format         | direct vendor output → AnnData/MuData                            |
| a multi-vendor normalised dataset format     | versioned declarative TOML parsing rules                         |
| a provenance / ontology / field-role design  | preserving exact vendor column roles in the rule                 |
| large-data engineering prior art             | parsed parameter-file semantics                                  |
|                                              | explicit quantification-level AnnData objects                    |
|                                              | APB-defined `.var` links (`ProForma_ion`, `ProForma_peptidoform`, `ProForma_peptide`, `Protein_Group`, `Protein_Ids`, `Protein_Names`, `Genes`) |

**Recommended direction:**

1. **Keep APB.** Its converter model and TOML rules are not what QPX provides.
2. **Borrow QPX's provenance / field-role / ontology metadata design** — this is the most directly
   reusable idea.
3. **Consider an APB → QPX export** only if QPX-ecosystem compatibility or Parquet-scale query
   becomes a real requirement.
4. **Do not replace APB's TOML converter** with QPX's per-vendor Python adapter model.

---

## References (external lineage)

- **quantms pipeline** — Dai, Füllgrabe, Pfeuffer, … Perez-Riverol et al., *quantms: a cloud-based
  pipeline for quantitative proteomics…*, **Nature Methods** (2024).
  doi:10.1038/s41592-024-02343-1
- **SDRF-Proteomics** — Dai et al., *A proteomics sample metadata representation for multiomics
  integration and big data analysis*, **Nat. Commun.** 12:5854 (2021).
  doi:10.1038/s41467-021-26111-3 — HUPO-PSI specification since 2023.
- **ProForma 2.0** — LeDuc, Deutsch, Binz, … Vizcaíno, *Proteomics Standards Initiative's ProForma
  2.0…*, **J. Proteome Res.** 21(4):1189–1195 (2022). doi:10.1021/acs.jproteome.1c00771
- **QPX repo** — github.com/bigbio/qpx (formerly `quantms.io`). Companion libraries: `mokume`
  (feature→protein quantification), `ibaqpy` (SDRF-driven absolute quantification).
- *No dedicated QPX/quantms.io publication was found as of 2026-06; the format is documented via
  the ecosystem above.*

---

## Appendix: Local Evidence Checked

Verified at `qpx@63d6382`:

- `qpx/README.md`, `qpx/docs/spec/converter-coverage.md`, `qpx/docs/spec/file-naming.md`
- `qpx/qpx/cli/convert.py`
- `qpx/qpx/converters/base.py`, `orchestrator.py`
- `qpx/qpx/converters/maxquant/{converter,feature_adapter,psm_adapter,pg_adapter}.py`
- `qpx/qpx/converters/diann/{converter,base_adapter,feature_adapter,pg_adapter}.py`
- `qpx/qpx/converters/spectronaut/converter.py`, `qpx/qpx/converters/fragpipe/{converter,feature_adapter}.py`
- `qpx/qpx/converters/mztab.py`, `qpx/qpx/converters/quantms/{converter,feature_adapter,pg_adapter}.py`
- `qpx/qpx/mudata.py`, `qpx/qpx/dataset.py`
- `qpx/qpx/core/{engine,sql}.py`
- `qpx/qpx/core/data/schemas/{feature,pg,psm,types}.yaml`, `qpx/qpx/config/column_mappings.yaml`
