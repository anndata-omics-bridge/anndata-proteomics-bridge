# anndata_proteomics

Convert proteomics software output to AnnData format.

Design lives in the sibling docs repo [anndata_omics_bridge](../anndata_omics_bridge/):
- **[conventions.md](../anndata_omics_bridge/docs/conventions.md)** — column / layer name sanitisation rules (apply on `obs.columns`, `var.columns`, layer names; **not** on `obs_names`/`var_names`/`uns` keys)
- **[adr_tool_specific_views.md](../anndata_omics_bridge/docs/adr_tool_specific_views.md)** — per-tool `uns['<app_name>']['column_roles']` schema (authoritative ADR)
- **[proteomics_rationale.md](../anndata_omics_bridge/docs/proteomics_rationale.md)** — why AnnData for proteomics; ProteoBench / prolfquapp synergies

In-repo docs: [docs/toml_schema.md](docs/toml_schema.md), [docs/RESTART_PLAN.md](docs/RESTART_PLAN.md).

## Current Scope

**Ion/precursor level quantification only:**
- DIA-NN (`report.tsv`)
- MaxQuant (`evidence.txt`)
- Spectronaut (precursor exports)

## Usage

```python
from anndata_proteomics.builder import ConverterBuilder

# Auto-detect format
converter = ConverterBuilder.from_file("report.tsv")
adata = converter.convert("report.tsv", "annotation.csv")

# Explicit software
converter = ConverterBuilder.for_software("diann")
adata = converter.convert("report.tsv", "annotation.csv")
```

## Project Structure

```
src/anndata_proteomics/
├── builder.py          # ConverterBuilder (auto-detect, for_software)
├── core.py             # Converter (pivot to AnnData)
├── annotation.py       # Sample annotation loading/matching
├── proforma.py         # ProForma sequence conversion
├── utils.py            # Utilities
├── strategies/         # One file per software
│   ├── diann.py
│   ├── maxquant.py
│   └── spectronaut.py
└── configs/            # ProForma TOML configs
    ├── diann.toml
    ├── maxquant.toml
    └── spectronaut.toml
```

## Strategy Interface

Each strategy defines:
- `name` - Software name (e.g., "DIA-NN")
- `obs_id` - Column identifying samples
- `var_id` - Column identifying precursors
- `VAR_COLUMNS` - Columns to include in var metadata
- `LAYER_COLUMNS` - Columns to include as layers (first = default X)
- `detect(path)` - Check if file matches this format
- `load(path)` - Load file to DataFrame
- `get_obs(df)` - Return obs DataFrame
- `get_var(df)` - Return var DataFrame
- `get_layers(df)` - Return DataFrame with obs_id, var_id, and layer columns

## Adding a New Strategy

1. Create `strategies/newsoftware.py`:
```python
class NewSoftwareStrategy:
    name = "NewSoftware"
    obs_id = "Run"
    var_id = "precursor_id"
    DETECTION_COLUMNS = ["RequiredCol1", "RequiredCol2"]
    VAR_COLUMNS = ["Sequence", "Charge", "Protein"]
    LAYER_COLUMNS = ["Intensity", "Score"]

    def detect(self, path): ...
    def load(self, path): ...
    def get_obs(self, df): ...
    def get_var(self, df): ...
    def get_layers(self, df): ...
```

2. Add import to `builder.py`:
```python
from .strategies.newsoftware import NewSoftwareStrategy

STRATEGY_REGISTRY = {
    ...
    _normalize_software_name(NewSoftwareStrategy.name): NewSoftwareStrategy,
}
```

## Test Data

ProteoBench test data:
- `/Users/wolski/projects/ProteoBench/test/data/quant/quant_lfq_ion_DIA_AIF/`
- `/Users/wolski/projects/ProteoBench/test/data/quant/quant_lfq_ion_DDA_QExactive/`

## Coding Rules

- **Keep `__init__.py` files empty** (a single module docstring is acceptable). Put classes/functions in separate modules and import them directly from those modules.

## Development

```bash
uv venv && source .venv/bin/activate
uv pip install -e .
pytest tests/
```
