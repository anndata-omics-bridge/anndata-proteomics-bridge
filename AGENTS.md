# anndata_proteomics

Convert proteomics software output to AnnData format.

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

Packaged parsing rules span **6 vendors** across the **ion / fragment / peptidoform / protein**
quantification levels:

- DIA-NN — ion, fragment, protein (version-specific rules under `diann/v1/`, `diann/v2/`)
- Spectronaut — ion, fragment, protein
- MaxQuant — ion (`evidence.txt`)
- FragPipe — ion
- PEAKS — ion
- WOMBAT — peptidoform

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

```bash
uv venv && source .venv/bin/activate
uv pip install -e ".[dev]"   # the [dev] extra brings in pytest + ruff
pytest tests/
```
