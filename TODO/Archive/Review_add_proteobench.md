# Review — add ProteoBench: substantive gaps

Plan reviewed: `TODO_add_proteobench.md`. Formulas, `intermediate_hash`
(`sha1(intermediate.to_string())`), and JSON layout check out against the
ProteoBench source. The plan's input contract can be simplified — most of it is
already solved by existing infrastructure.

> Archived 2026-07-24. The reviewed HYE implementation and parity work are
> complete; PYE/plasma is tracked separately.

## Design (settled — not a gap)

- Roles are per-tool `uns['<app>']['column_roles']` (omics-bridge ADR). `convert`
  stays tool-agnostic and does **not** write any tool namespace — APB cannot
  enumerate every present/future consumer.
- APB gives as much uniformity as reasonable (uniform `X`/`obs`/`var` with vendor
  columns preserved, `software_name`/level in `uns`). Tool developers still own
  which columns their tool needs — same story coming for prolfquapp.
- ProteoBench's per-tool parse-settings TOMLs already carry the per-tool column
  knowledge (which column → `Proteins`, `decoy_flag`, `contaminant_flag`). Nothing
  to migrate into APB rules, no role for APB to declare.

So: **`apb proteobench <obj> <module_settings.toml> <per-tool.toml>`** — two
ProteoBench inputs, the split ProteoBench already uses (experiment design vs.
per-tool parsing). `apb proteobench` writes its own `uns['proteobench']`.

## Resolved by existing infrastructure

- **Name reconciliation.** ProteoBench's per-tool TOML names raw vendor columns
  (`Protein.Ids`); APB `var` holds sanitized names (`Protein_Ids`). `convert`
  already stores the full merged rule as a JSON string at
  `uns["anndata_proteomics"]["rule_json"]` (assemble.py:37-43). `apb proteobench`
  inverts its `columns.var.select` map (source → APB var name) to resolve the
  exact converted column. No sanitization heuristics.
- **ROC-AUC.** Use `sklearn.metrics.roc_auc_score` (ProteoBench does) for
  guaranteed parity: `y_score = |log2_A_vs_B|`, drop NaNs, `NaN` on single class,
  unchanged species = smallest `|log2_expectedRatio|`.

## The one substantive gap: column retention

`rule_json` can only map columns that were actually `select`ed into `var`.
MaxQuant's `rules.json` selects `Proteins`/`Leading_Razor_Protein` but **not**
`Reverse` or `Potential contaminant`, so the decoy/contaminant filter can't be
reproduced from the converted object. Fix = add the columns each vendor's
per-tool TOML references to that vendor's `rules.json` `select`, so they are
retained in `var` and appear in `rule_json`. Not a role, not a contract — just
retaining the referenced columns. Audit all six vendors before declaring each
scoreable.

## Parity must-checks before the first golden test

- **`proteobench_version`** in the compatible score doc is stamped from the
  installed package, not derivable from formulas. Decide what APB emits and
  whether it is compared against the golden JSON.
- **`annotations/*.toml`** (the downloaded module TOMLs `apb annotate` reads):
  verify they carry the full scoring subset (`species_mapper`,
  `species_expected_ratio`, `min_count_multispec`), not just `[[samples]]`.
