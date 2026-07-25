# Canonical AnnData ProteoBench scoring

Status: completed on 2026-07-24.

## Problem

`apb proteobench` currently requires a ProteoBench per-tool parsing TOML after
APB has already converted the vendor output. That repeats raw-format knowledge
downstream, makes scoring depend on a second vendor mapping, and can reject a
valid canonical AnnData object because the TOML names a different raw source
column than the APB rule used for `X`.

The concrete DIA-NN DDA Astral case exposes the problem:

- DIA-NN reports both `Ms1.Normalised` and `Precursor.Normalised`;
- the ProteoBench DDA Astral DIA-NN adapter uses `Ms1.Normalised`;
- APB DIA-NN v2 currently puts `Precursor.Normalised` in `X`;
- the two columns contain different DIA-NN quantities, so selecting the correct
  source belongs in conversion rather than in scoring.

## Approved design

Canonical AnnData is the scoring contract:

```text
vendor output -> APB parsing rule -> canonical AnnData/MuData -> scoring
```

- `X` is the quantitative matrix consumed by scoring.
- canonical ProForma feature columns are the scoring feature identifiers.
- the stored APB rule supplies the minimal semantic metadata needed to locate
  protein accessions; scoring must not rediscover vendor columns.
- the ProteoBench module TOML remains the external experiment contract for
  samples, conditions, species mappings, expected ratios, and thresholds.
- no ProteoBench per-tool parsing TOML or per-tool settings cache participates
  after conversion.

## Implementation

- [x] Change the DIA-NN v2 ion rule so `Ms1.Normalised` is the `x_layer`, while
      retaining `Precursor.Normalised` as an auxiliary reported layer.
- [x] Add the smallest rule-schema semantic role needed to identify the
      canonical protein-accession column and validate that it names a declared
      `var` column.
- [x] Populate that role in every packaged rule currently covered by
      ProteoBench scoring regressions.
- [x] Refactor role resolution, run alignment, intermediate construction, and
      `score_quantification()` to consume canonical AnnData plus module
      settings only.
- [x] Change `apb proteobench` to accept `DATA MODULE_SETTINGS` without a
      per-tool settings argument.
- [x] Remove per-tool settings downloads, lookup APIs, cache paths, and their
      documentation/tests from APB.
- [x] Remove `tool_settings` from APB Studio run snapshots, registry resources,
      target commands, prerequisite diagnostics, and tests.
- [x] Preserve ProteoBench-compatible scoring for existing audited DIA-NN,
      FragPipe, MaxQuant, and WOMBAT paths, with DIA-NN results intentionally
      recomputed from `Ms1.Normalised`.

## Verification

- [x] All packaged parsing-rule JSON and generated schemas validate.
- [x] Every cached DIA-NN v2 fixture still converts; `X` is sourced from
      `Ms1.Normalised`.
- [x] The cached DDA Astral DIA-NN submission scores without a tool TOML.
- [x] Canonical and golden ProteoBench scoring regressions pass.
- [x] APB and APB Studio fast and full quality gates pass.
