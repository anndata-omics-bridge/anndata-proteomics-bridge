# anndata_proteomics_bridge

Convert proteomics quantification outputs to AnnData using declarative TOML parsing rules.

> **Status: restart core complete.** Read a vendor file, validate a TOML rule, convert to AnnData end-to-end. See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for the module map and [docs/RESTART_PLAN.md](docs/RESTART_PLAN.md) for the roadmap.

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
# Validate one or more TOML parsing rules. With no path: walks all packaged rules.
anndata-proteomics validate                                  # all packaged
anndata-proteomics validate path/to/your_rule.toml           # one (or several)

# List packaged rules.
anndata-proteomics list

# Convert a vendor file to AnnData (.h5ad). The rule is auto-recognized from the
# data's column headers; pass --rule-toml to override.
anndata-proteomics convert path/to/report.tsv                # writes report.h5ad next to it
anndata-proteomics convert report.tsv --output out.h5ad
anndata-proteomics convert report.tsv --rule-toml my_rule.toml --output out.h5ad

# Regenerate the JSON Schema after editing the pydantic models.
anndata-proteomics export-schema
```

Exit codes: `0` all pass, `1` validation / conversion failed.

## Report generation

Generate conversion reports for the packaged parsing rules with:

```bash
python tools/generate_report.py
```

By default, outputs are written under:

```text
examples/results/
```

The report index is therefore:

```text
examples/results/index.html
```

The index includes one row per packaged rule, with links to the input file, generated
`.h5ad`, rendered HTML report, per-rule log, and the input/`.h5ad` file sizes.

Use `--output-dir` to write somewhere else:

```bash
python tools/generate_report.py --output-dir path/to/outdir
```

That writes:

```text
path/to/outdir/index.html
```

Useful filters:

```bash
python tools/generate_report.py --rule DIA-NN --rule WOMBAT
python tools/generate_report.py --log-level DEBUG
```

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
