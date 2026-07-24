# Changes

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
