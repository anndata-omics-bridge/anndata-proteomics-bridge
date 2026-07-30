<!-- Managed by agent: keep commands and file references verified -->
<!-- Last updated: 2026-07-23 | Last verified: 2026-07-23 -->

# anndata_proteomics

Convert proteomics software output to AnnData format.

**Precedence:** the closest `AGENTS.md` to changed files wins. Explicit user
instructions override repository files.

## Terminology

- **APB** means this project, `anndata_proteomics_bridge`.
- The Python package remains `anndata_proteomics`; use **APB** only as the project-level
  shorthand in plans, architecture notes, and cross-repo migration discussions.

Design lives in the sibling docs repo [anndata_omics_bridge](../anndata_omics_bridge/):
- **[conventions.md](../anndata_omics_bridge/docs/conventions.md)** — column / layer name sanitisation rules (apply on `obs.columns`, `var.columns`, layer names; **not** on `obs_names`/`var_names`/`uns` keys)
- **[adr_tool_specific_views.md](../anndata_omics_bridge/docs/adr_tool_specific_views.md)** — per-tool `uns['<app_name>']['column_roles']` schema (authoritative ADR)
- **[proteomics_rationale.md](../anndata_omics_bridge/docs/proteomics_rationale.md)** — why AnnData for proteomics; ProteoBench / prolfquapp synergies

In-repo docs: [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md), [docs/json_schema.md](docs/json_schema.md).

JSON parsing-rule edits must follow [docs/json_schema.md](docs/json_schema.md).
In particular, `columns.*.select` is for original input-table columns only
(`"<sample>"` is the wide-file exception); APB-derived values such as
`proforma_sequence` and `stripped_sequence` must be declared via
`columns.var.compute`.

Every vendor/version group is one self-contained `rules.json` document with the same
shape: document metadata, a shared `base`, and a `levels` map. Put fields shared by every
level (`modifications`, `axis.obs_keys`, common scalars) in `base`; put level-specific
axis keys/layers/computes in `levels.<level>`. `rules.loader.load_rule` merges the selected
level over its same-document base before Pydantic validation. There are no `$extends`
paths or separate base files. Single-level vendors use the same document shape.

## Current Scope

Packaged parsing rules span **8 vendors** across the **ion / fragment / peptidoform / protein**
quantification levels:

- DIA-NN — ion, fragment, protein (version-specific rules under `diann/v1/`, `diann/v2/`)
- Spectronaut — ion, fragment, protein
- MaxQuant — ion (`evidence.txt`, 1.5.x through 2.7.x)
- FragPipe — ion
- PEAKS — ion
- Sage — ion, peptidoform (wide `lfq.tsv`; the level is parameter-gated, see below)
- AlphaPept — ion (long, comma-delimited PSM table)
- WOMBAT — peptidoform

Sage's `lfq_settings.combine_charge_states` (default `true`) collapses charge states and writes
`charge = -1`, so the same `lfq.tsv` schema is ion-level or peptidoform-level depending only on that
setting. The Sage document therefore gates each level on `requires_search_parameters` rather than a
version regex; see [docs/json_schema.md](docs/json_schema.md).

See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for the authoritative vendor/level/shape table.

## Status

The restart core (`vendor file + parsing JSON → AnnData`) is **complete**, and the package has
grown beyond it: vendor **parameter parsing** (`params/`), modified-sequence **normalisation**
(`modifications/`), and second-stage **annotation** (`annotation/` — obs joins and FASTA-driven
protein `varm['fasta']`). [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) is the current module map and the JSON rule
schema lives in [docs/json_schema.md](docs/json_schema.md). The pre-restart `src/` was deleted on
2026-05-01 and is recoverable from git history (last full commit before deletion: `f6bffda`).

**apb is a pure library + `apb` CLI — no GUI, no marimo dependency.** On 2026-06-28 the marimo
tooling (test-data browser, AnnData viewer, background-job runner, ProteoBench catalog) was
extracted to the sibling **`apb_studio`** package, which drives apb through the `apb` CLI. The
conversion core that `apb convert` orchestrates lives in `converters/pipeline.py` (`scripts/` now
holds only `cli.py`). Do not reintroduce marimo or a GUI here.

## Test Data

`test_data_download/` is the single canonical local ProteoBench cache. It is
generated with `apb-testdata catalog/select/download/fasta`, consumed by the
integration tests, and gitignored. Do not add parallel `benchmark_data/`,
`examples/`, or test-data download scripts.

**Start here instead of searching the cache** — these three CSVs index everything
that is catalogued, selected, and actually on disk:

| File | Contents |
| --- | --- |
| `test_data_download/raw_file_db_full.csv` | Catalog of every ProteoBench submission: `module`, `repo_name`, `intermediate_hash`, `software_name`, `software_version`, `nr_feature` |
| `test_data_download/raw_file_db_selected.csv` | The subset selected for download (same columns) |
| `test_data_download/raw_file_db_downloaded.csv` | Manifest of what is on disk: the same columns plus `input_file_path`, `input_file_size_bytes`, `status` |

Downloaded submissions live at
`test_data_download/json_dir/<Results_repo>/<intermediate_hash>/`, each holding the
vendor result (`input_file.*`), the vendor parameter file (`param_0..txt`),
`result_performance.csv`, and `comment.txt`. To find an example for one
vendor/module combination, query the manifest CSV — do not glob the tree.

## Coding Rules

- **Keep `__init__.py` files empty** (a single module docstring is acceptable). Put classes/functions in separate modules and import them directly from those modules.
- **APB owns reusable proteomics parsing infrastructure.** Modification cleanup/mapping
  rules currently duplicated in ProteoBench per-tool configs should migrate into APB parsing
  JSON/schema instead of being reimplemented downstream.
- **Parameter parsing belongs in APB.** ProteoBench parameter parsers should move into APB
  as shared code; ProteoBench should consume APB rather than remain the upstream owner of
  generic vendor parameter parsing.
- **Peptide-to-protein algorithms belong in Prozor.** APB may orchestrate Prozor over FASTA
  records and store results in AnnData/MuData, but must not vendor Aho--Corasick or protein
  inference implementations. The dependency direction is APB → Prozor, never Prozor → APB.

## Development

| Task | Command |
| --- | --- |
| Install | `uv sync --frozen --extra dev --group docs` |
| Fast checks | `uv run pre-commit run --hook-stage pre-commit --all-files` |
| Full gate | `uv run pre-commit run --hook-stage pre-push --all-files` |
| Single test | `uv run pytest tests/test_cli.py -q` |
| Security audit | `uv run pre-commit run dependency-audit --hook-stage manual --all-files` |

The pre-commit configuration is the command source of truth for CI. Do not
lower Ruff, strict Pyright, dependency, or coverage gates without explicit
approval.

## Scoped AGENTS.md

- [GitHub workflows](./.github/workflows/AGENTS.md)
