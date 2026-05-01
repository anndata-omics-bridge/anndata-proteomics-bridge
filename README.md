# anndata_proteomics_bridge

Convert proteomics software output to AnnData format.

Currently supported (ion / precursor level only):

- DIA-NN (`report.tsv`)
- MaxQuant (`evidence.txt`)
- Spectronaut (precursor exports)

Design rationale, role separation, naming conventions, and the per-tool `uns['<app_name>']['column_roles']` schema live in the sibling docs repo: [anndata_omics_bridge](../anndata_omics_bridge/).

## Install

```bash
uv venv
source .venv/bin/activate
uv pip install -e .
```

## Quick start

```python
from anndata_proteomics.builder import ConverterBuilder

# Auto-detect format
converter = ConverterBuilder.from_file("report.tsv")
adata = converter.convert("report.tsv", "annotation.csv")

# Explicit software
converter = ConverterBuilder.for_software("diann")
adata = converter.convert("report.tsv", "annotation.csv")
```

## Sample annotation file

A CSV with at least one column matching the software's run identifier (e.g. DIA-NN's `Run` or `File.Name`):

```csv
sample_id,condition,batch,replicate
sample_001,control,batch1,1
sample_002,control,batch1,2
sample_003,treated,batch1,1
sample_004,treated,batch1,2
```

## Tests

```bash
pytest tests/
```

## Documentation

- **[anndata_omics_bridge/docs/conventions.md](../anndata_omics_bridge/docs/conventions.md)** — column / layer name sanitisation rules
- **[anndata_omics_bridge/docs/adr_tool_specific_views.md](../anndata_omics_bridge/docs/adr_tool_specific_views.md)** — per-tool `uns` design (authoritative ADR)
- **[anndata_omics_bridge/docs/proteomics_rationale.md](../anndata_omics_bridge/docs/proteomics_rationale.md)** — why AnnData for proteomics
- **[docs/diann_mapping.md](docs/diann_mapping.md)** — DIA-NN-specific conversion details
- **[docs/toml_schema.md](docs/toml_schema.md)** — TOML rules schema for per-tool parsers
- **[docs/RESTART_PLAN.md](docs/RESTART_PLAN.md)** — current implementation roadmap
