# anndata_proteomics_bridge

Convert proteomics quantification outputs to AnnData using declarative TOML parsing rules.

> **Status: restart in progress.** The TOML rule schema, loader, registry, and validation are implemented. Vendor-file reading and AnnData conversion are still TODO. See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for what's implemented now and [docs/RESTART_PLAN.md](docs/RESTART_PLAN.md) for the roadmap.

## Packaged parsing rules

| Vendor | Format | Level |
|---|---|---|
| DIA-NN | long | ion |
| Spectronaut | long | ion |
| MaxQuant (`evidence.txt`) | long | ion |
| FragPipe (`combined_modified_peptide.tsv`) | wide | ion |
| PEAKS | wide | ion |
| WOMBAT | wide | peptidoform |

## Install

```bash
uv venv --python 3.13
source .venv/bin/activate
uv pip install -e '.[dev]'
```

## Quick start

The umbrella CLI is `anndata-proteomics`:

```bash
anndata-proteomics validate                                # walks all packaged rules
anndata-proteomics validate path/to/your_rule.toml         # validates a single TOML
anndata-proteomics list                                    # show packaged rules
anndata-proteomics export-schema                           # regenerate parse_rule.schema.json
anndata-proteomics convert data.tsv rule.toml              # STUB until step 5+ lands
```

Exit codes: `0` all pass, `1` validation failed, `2` not yet implemented (e.g. `convert`).

Load a packaged rule from Python:

```python
from anndata_proteomics.rules.loader import load_packaged_rule

rule = load_packaged_rule("diann", "ion")
print(rule.software_name, rule.input_shape, rule.axis.x_layer)
# DIA-NN long Precursor_Normalised
```

Validate a custom TOML rule programmatically:

```python
from anndata_proteomics.rules.validate import validate_file

result = validate_file("/path/to/parse_<software>_<level>_<v>.toml")
if not result.ok:
    print(result.error)
```


## Tests

```bash
pytest tests/
```

## Documentation

- **[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)** — current module map, public API, data flow
- **[docs/RESTART_PLAN.md](docs/RESTART_PLAN.md)** — restart roadmap and step-by-step plan
- **[docs/toml_schema.md](docs/toml_schema.md)** — TOML parsing-rule schema spec
- **[anndata_omics_bridge/docs/conventions.md](../anndata_omics_bridge/docs/conventions.md)** — column / layer name sanitisation rules
- **[anndata_omics_bridge/docs/adr_tool_specific_views.md](../anndata_omics_bridge/docs/adr_tool_specific_views.md)** — per-tool `uns` design (authoritative ADR)
- **[anndata_omics_bridge/docs/proteomics_rationale.md](../anndata_omics_bridge/docs/proteomics_rationale.md)** — why AnnData for proteomics
