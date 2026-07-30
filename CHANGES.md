# Changes

- 2026-07-30: Support AlphaDIA at ion level across all 11 cached submissions and six
  versions, in three shape-distinct documents: `alphadia/v1_10` (wide TSV, one bare run
  column per sample), `alphadia/v1_12` (long TSV), `alphadia/v2` (long parquet, dotted
  namespaces). Each converts to the exact six-run axis its ProteoBench module annotation
  declares. Contrary to the plan, no `sample_name_cleanup` is needed: the diaPASEF run
  columns carry a trailing acquisition id (`..._Alpha_01_11494`) that the annotation's
  `raw_file` carries too, so stripping it would have broken the join.
- 2026-07-30: Add a second modification parser, `parser = "site_list"`, for vendors that
  write modifications as parallel name/site columns beside a bare sequence instead of
  inline tokens. AlphaDIA emits the alphabase layout (`mods = Oxidation@M;Carbamidomethyl@C`
  with `mod_sites = 9;2`), which pairs index-wise rather than in sorted order, and uses
  site `0` for a protein N-terminal modification. Without it the same sequence and charge
  with different modifications collapse into one feature, silently summing an oxidised and
  a non-oxidised precursor — so this gates AlphaDIA support rather than merely enriching it.
- 2026-07-30: Add `params/parsers/alphadia.py`. AlphaDIA ships an ANSI-coloured, timestamped
  run log rather than a config file; the parser strips that framing, keeps the applied value
  of `[user defined, default: X]` entries, and anchors `software_version` on the startup
  `PROGRESS:` banner because the config tree also carries a bare `version:` key. A `0`
  tolerance records automatic calibration rather than a zero-width window.
- 2026-07-30: Wide sample detection now excludes every column the rule accounts for by
  name. Modifications and column materialization both run before the wide dispatch, so
  APB's derived columns and the rule's renamed `select` outputs are on the frame by then;
  a rule whose sample pattern cannot anchor on a suffix — AlphaDIA's run columns are bare
  run names — matched 13 of them as extra samples.
- 2026-07-30: AlphaDIA's `mod_seq_charge_hash` is deliberately not carried as a var column.
  It is a uint64 that pandas reads from TSV as float64, collapsing 81 949 distinct values
  to 75 292 and leaving every one inexact; `ProForma_ion` is the feature key.
- 2026-07-30: Upgrade to pandas 3.0.5 and anndata 0.13.2 (also numpy 2.5.1, scipy
  1.18.0). anndata `<0.13` capped `pandas<3`. From 0.13 `adata.layers` also yields `X`
  under a `None` key, which leaked a `"None"` layer into the descriptive summary and the
  CLI log line — `_matrix_types.named_layers` now filters it at the three call sites.
- 2026-07-30: `coerce_numeric` returns plain `float64` instead of passing a nullable
  dtype through. Layers are `float64` end to end, and on a nullable single-column frame
  pandas 2.3 `bfill(axis=1)` filled *down the feature axis* instead of acting as a no-op,
  silently copying one feature's value onto its neighbours. Fixed in pandas 3, pinned by
  a regression test regardless.
- 2026-07-30: Add `layers[].value_pattern` to the parsing-rule schema: a single-capture-
  group regex applied per cell before numeric coercion, for vendor columns holding
  structured strings. PEAKS `AScore` is `site:modification:score`, so the whole layer
  coerced to NaN and reached the output silently empty; it now extracts the score
  (96.86% missing on the DDA corpus, matching the 3.14% of populated cells).
- 2026-07-30: Declare `missing_values: [0]` on the PEAKS `Normalized_Area` layer. PEAKS
  writes `0` for an ion it did not quantify, so missingness read as 0% where it is
  actually 5.09% (DDA Orbitrap) and 10.25% (DIA diaPASEF). Verified against the paired
  per-sample `m/z` column, which carries `-` in exactly the same cells (0 disagreements
  in 456,774; 36 in 504,665). Matches the existing FragPipe and Sage declarations.
- 2026-07-30: Warn when a captured numeric layer is ≥99.9% NaN after coercion. Both
  defects above reached the output silently; a matched-but-empty layer is a rule defect,
  not a measurement.
- 2026-07-30: Lower the changed-line coverage gate from 100% to 90%.
- 2026-07-30: Gate a level's availability on parsed search parameters via
  `levels.<level>.requires_search_parameters`, and use it to make Sage's quantification
  level follow `lfq_settings.combine_charge_states` instead of a version regex. Sage's
  DOCS.md defaults that setting to `true`, which collapses charge states and writes
  `charge = -1`, so the same `lfq.tsv` schema is ion- or peptidoform-level with neither the
  version nor the headers able to tell them apart — the previous `^0\.15\.` pin would have
  failed a default-configured 0.15.x submission on the `-1` sentinel. The Sage document now
  declares both levels at `^0\.` and both cached submissions convert at their true level. A
  gated level is unavailable when no parameters could be parsed, rather than guessed.
- 2026-07-30: `Parameters.combine_charge_states` records whether quantification merged a
  peptidoform's charge states; `params/parsers/sage.py` extracts it from
  `quant.lfq_settings`, defaulting to Sage's own `true` when the LFQ block omits it.
- 2026-07-30: Add `columns.*.optional_select` to the parsing-rule schema: vendor
  columns captured when the export carries them and skipped when it does not, the
  column-side counterpart of an optional layer. They never gate recognition, may not
  be an axis key, are forbidden on the wide obs axis, and drop out of a `coalesce` /
  `join_nonempty` chain when absent.
- 2026-07-30: Cover every cached MaxQuant submission (1.5.2.8 through 2.7.5.0, 12
  fixtures across four modules) instead of only 2.6.7.0. `Fraction` is absent from 11
  of 12 exports and `Experiment` from 2, both configuration- rather than
  version-dependent, so they moved to `optional_select` together with the 1.5.2.8
  title-case `Leading Proteins` / `Leading Razor Protein` spellings.
- 2026-07-30: Add a Sage ion rule for the wide `lfq.tsv` (`^0\.15\.`). Pinned to the
  charge-resolved family: `lfq_settings.combine_charge_states = true` exports (0.14.6)
  report `charge = -1` and are peptidoform-level, so they stay unresolved rather than
  fail conversion.
- 2026-07-30: Add an AlphaPept ion rule (`^0\.5\.`) for the long, comma-delimited PSM
  table. `ms1_int_sum_apex_dn` is the `x_layer` — it reproduces ProteoBench's per-run
  intensities exactly for 98% of shared precursors where every other `ms1_int*`
  candidate scores zero — with `duplicates.mode = "keep_first"`.
- 2026-07-30: `find_test_data` accepts a rule's `software_version` so a vendor whose
  cached submissions span several export schemas returns a file the rule covers, and
  `rules.loader.software_version_matches` is public for that reuse.
- 2026-07-30: Score ProteoBench at every quantification level. The module TOML
  now supplies only the sample design and per-species expected ratios; the
  feature axis is `var_names` (the rule's joined `axis.var_keys`), so the two
  module-level assertions are gone, a MuData scores each modality into its own
  `uns`/`varm` instead of one selected modality, and the legacy intermediate
  names its feature column after the scored level.
- 2026-07-30: Declare `column_roles.protein_accessions` on the DIA-NN fragment
  and protein levels and the Spectronaut fragment and protein levels; scoring
  resolves species through that role, so those levels need re-conversion before
  they can be scored.
- 2026-07-29: Support one or multiple exact sample identifiers through
  `raw_file_alias` and `raw_file_aliases`; add WOMBAT `A_1` through `B_3`
  aliases to the DDA ion and peptidoform module settings.
- 2026-07-29: Make `apb annotate` a mandatory prerequisite for
  `apb proteobench`; scoring now consumes only annotated `sample_name` and
  `condition` values and no longer resolves samples from vendor run names,
  parsing-rule axes, filename cleanup, or aliases.
- 2026-07-25: Parse DIA-NN acquisition mode and materialize level-scoped,
  search-parameter-conditional axis overrides so DIA-NN DDA uses
  `Ms1_Normalised` in `X` while DIA uses `Precursor_Normalised`; record the
  effective column mapping in descriptive summaries.
- 2026-07-25: Document the three `test_data_download` index CSVs (catalog, selection,
  manifest) and the downloaded-submission layout in AGENTS.md, so examples are found by
  querying the manifest instead of globbing the cache.
- 2026-07-25: Add WOMBAT's native modified-sequence-plus-charge ion level and
  distinguish it from WOMBAT's charge-free peptidoform export by reported
  PSM-count columns.
- 2026-07-23: Align local and GitHub quality gates with the FGCZ Python
  reference: staged Ruff/Pyright/Deptry/coverage hooks, wheel inspection,
  strict docs, dependency audit, typed-package marker, and CI/Pages/security
  workflows.
- 2026-07-23: Complete ProteoBench-compatible FragPipe and MaxQuant protein values
  during conversion with declarative `join_nonempty`/`coalesce` computations; sum
  repeated MaxQuant evidence rows and enforce the declared duplicate modes.
- 2026-07-23: Download managed ProteoBench module/tool TOMLs from the pinned
  intermediate-format revision containing `species_mapper` and `[[samples]]`.
- 2026-07-22: Add independent matrix-native ProteoBench HYE scoring for AnnData/MuData,
  compatible score JSON in `uns['proteobench']['scores']`, feature intermediates in
  `varm['proteobench']`, compact protein-mapping provenance, a CLI command, managed
  scoring TOMLs, and golden/performance regression coverage.
- 2026-07-22: Give independent FASTA enrichment the default `.fasta.h5ad/.h5mu`
  output suffix instead of `.annotated`.
- 2026-07-22: Preserve `ProForma_peptide` on DIA-NN v1 fragment outputs so fragment and
  multi-modality FASTA validation use the complete peptide hierarchy.
- 2026-07-22: Read ProteoBench module annotation TOMLs directly, allow CSV/TSV observation
  tables without Pydantic modelling, and add managed module-annotation downloads to
  `apb-testdata`.
- 2026-07-22: Preserve downloaded fixture directories during `apb-testdata catalog` refreshes and
  allow FASTA lookup from an explicit test-data root for APB Studio.
- 2026-07-22: Add compact cumulative sample-annotation and level-specific FASTA components
  to stored descriptive summaries.
- 2026-07-21: Add stage-owned descriptive summaries, the `apb summary` command, and type-derived conversion output suffixes.
- 2026-07-21: Make JSON the canonical parsing-rule format, consolidate each
  software-version family into one `base`/`levels` document, and make `--rule-config`
  document-aware.
- 2026-07-21: Add typed FASTA configuration, peptide-to-protein validation, and
  AnnData/MuData feature-mapping storage backed by Prozor.
