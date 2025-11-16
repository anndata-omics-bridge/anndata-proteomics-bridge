# DIANN to AnnData Converter

## Project Overview

Convert DIANN (DIA-NN) proteomics output to AnnData format following the **AnnData Omics Bridge Specification**.

**Goal**: Create a Python package that converts DIANN parquet/TSV files to AnnData `.h5ad` files compatible with downstream analysis tools (exploreDE, prolfqua).

## Key Resources

### Specification Reference
- **Main Specification**: `/Users/wolski/projects/anndata_omics_bridge/docs/AnnData_Omics_Bridge_spec.qmd`
- **Mapping Guide**: `docs/DIANN_to_AnnData_mapping.md` (this project)

### Dependencies
- **omicsbridge**: Validation and column resolution utilities
  - Location: `/Users/wolski/projects/anndata_omics_bridge`
  - Import: `from omicsbridge import validate_anndata_omics, ColumnResolver`

## Project Structure

```
diann_to_anndata/
├── docs/
│   └── DIANN_to_AnnData_mapping.md   # Column mapping reference
├── src/diann_converter/
│   ├── __init__.py                   # Package exports
│   ├── converter.py                  # Main conversion logic
│   ├── reader.py                     # DIANN file reader
│   ├── annotation.py                 # Sample annotation handling
│   └── qc.py                         # Quality control and filtering
├── tests/
│   └── test_converter.py             # Unit tests
├── examples/
│   ├── sample_annotation.csv         # Example annotation file
│   └── convert_diann.py              # Example usage script
├── pyproject.toml                    # Package configuration
├── README.md                         # User documentation
└── CLAUDE.md                         # This file (AI assistant instructions)
```

## AnnData Requirements (exploreDE-compatible)

### Required Components

**var (features/proteins)**:
- `description`: Free-text protein description (REQUIRED)
- `label`: Gene symbols, protein IDs (optional)

**obs (samples)**:
- `factor`: Experimental factors like condition, batch (REQUIRED)
- `label`: Sample identifiers (optional)

**X (abundance matrix)**:
- Log2-transformed intensities (analysis-ready)

**layers**:
- `intensities`: Raw DIANN quantities

**uns (metadata)**:
```python
adata.uns['exploreDE'] = {
    'column_roles': {
        'var': {
            'description': ['description'],
            'label': ['gene_symbol', 'protein_id']
        },
        'obs': {
            'factor': ['condition', 'batch'],
            'label': ['run', 'file_name']
        }
    }
}
```

## Development Setup

### Using uv for Virtual Environment

This project uses `uv` for fast Python package management.

**Setup:**
```bash
# Create virtual environment
uv venv

# Activate virtual environment
source .venv/bin/activate  # On macOS/Linux
# .venv\Scripts\activate  # On Windows

# Install dependencies
uv pip install -e .

# Install optional dependencies for validation
uv pip install -e "/Users/wolski/projects/anndata_omics_bridge"
```

**Running examples:**
```bash
# Activate venv first
source .venv/bin/activate

# Run DIA-NN example
python examples/convert_diann_example.py

# Run Spectronaut example
python examples/convert_spectronaut_example.py
```

**Running tests:**
```bash
source .venv/bin/activate
pytest tests/
```

## Development Principles

1. **Self-contained code blocks**: Each example should import all necessary dependencies
2. **Validation**: Always validate output with `omicsbridge.validate_anndata_omics()`
3. **Error handling**: Provide clear error messages for common issues:
   - Missing annotation file
   - Sample ID mismatches
   - Invalid DIANN file format
4. **Flexible column mapping**: Allow users to specify which DIANN columns map to which roles

## Common Tasks

### Reading DIANN Files
- Support both parquet and TSV formats
- Handle different DIANN versions (column name variations)
- Extract protein-level, peptide-level, or precursor-level data

### Sample Annotation
- Match DIANN `Run` or `File.Name` to annotation file
- Handle filename variations (paths, extensions)
- Validate that all samples have annotations

### Data Reshaping
- Pivot long format → wide format (samples × proteins)
- Handle missing values appropriately
- Compute summary statistics (mean intensity, n_peptides)

### Quality Control
- Filter by protein Q-value
- Remove low-quality proteins
- Handle missing values (filtering vs imputation)

## Testing Strategy

1. **Unit tests**: Test individual functions (reshape, filter, validate)
2. **Integration tests**: Test full conversion pipeline
3. **Test data**: Use small synthetic DIANN output
4. **Validation**: All outputs must pass `validate_anndata_omics()`

## Code Style

- Follow PEP 8
- Use type hints
- Document functions with docstrings
- Keep functions focused and modular

## Version History

- v0.1.0: Initial project setup
