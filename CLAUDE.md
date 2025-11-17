# Proteomics to AnnData Converter

## Project Overview

Convert proteomics quantification software output to AnnData format following the **AnnData Omics Bridge Specification**.

**Goal**: Create a Python package that converts output from multiple proteomics software tools to AnnData `.h5ad` files compatible with downstream analysis tools (exploreDE, prolfqua).

### **IMPORTANT: Initial Scope - Ion/Precursor Level Quantification**

**In the first pass, this project focuses EXCLUSIVELY on ion/precursor level quantification data.**

This means:
- ✅ **Ion/precursor quantification** (e.g., DIA-NN report.tsv, MaxQuant evidence.txt)
- ❌ **NOT protein-level quantification** (defer to later phases)
- ❌ **NOT peptidoform-level quantification** (defer to later phases)

The initial implementation should follow ProteoBench's `parse_ion.py` approach, which handles precursor/ion-level data from various software tools. Protein-level aggregation and peptidoform-level data will be added in subsequent development phases.

### **Key Differences from ProteoBench**

While using ProteoBench as a reference implementation, this project has important distinctions:

1. **Comprehensive var (ion/precursor) Metadata Extraction**:
   - **ProteoBench**: Extracts minimal information needed for benchmarking
   - **This project**: Extract and preserve **ALL** ion/precursor-level information provided by the software
   - Examples of information to capture:
     - Modified sequence, stripped sequence
     - Charge state, m/z, retention time
     - All quality scores (Q-values, PEP, confidence scores)
     - Protein associations (protein IDs, gene symbols, protein names)
     - Modification details
     - Fragment ion information (if available)
     - Software-specific metadata (e.g., DIA-NN's `Q.Value`, `PG.Q.Value`, `GG.Q.Value`)

2. **Multiple Quantification Types in Layers**:
   - **ProteoBench**: Typically extracts single quantification value
   - **This project**: Extract **ALL** quantification outputs and store each in separate layers
   - Examples:
     - DIA-NN: Store both `Ms1.Area` and `Precursor.Quantity` if available
     - MaxQuant: Store `Intensity` and `Area` separately
     - Spectronaut: Store different normalization states if present
   - Naming convention: Normalize original names (see Column Naming Policy below)

3. **Column Naming Policy - Normalize, Don't Rename**:
   - **ProteoBench**: May standardize column names across software
   - **This project**: Preserve original column names, only normalize them
   - **Normalization rules** (minimal changes):
     - Replace whitespace with underscores: `"Raw file"` → `"raw_file"`
     - Remove or replace special characters: `"Modified.Sequence"` → `"modified_sequence"` (dots to underscores)
     - Convert to lowercase for consistency
     - Ensure valid Python identifiers (no leading numbers, reserved keywords)
   - **DO NOT rename** columns to different semantic names:
     - ❌ Don't rename `"Run"` to `"sample_id"`
     - ❌ Don't rename `"Modified.Sequence"` to `"peptide_sequence"`
     - ✅ Keep as `"run"` and `"modified_sequence"` (normalized only)
   - **Rationale**: Preserve traceability to source software documentation
   - **Implementation guidance**:
     ```python
     def normalize_column_name(col_name: str) -> str:
         """Normalize column name: lowercase, underscores, valid Python identifier."""
         # Convert to lowercase
         normalized = col_name.lower()
         # Replace whitespace and special chars with underscores
         normalized = re.sub(r'[\s\.\-]+', '_', normalized)
         # Remove other special characters
         normalized = re.sub(r'[^\w]', '', normalized)
         # Ensure doesn't start with number
         if normalized[0].isdigit():
             normalized = 'col_' + normalized
         return normalized

     # Examples:
     # "Modified.Sequence" → "modified_sequence"
     # "Raw file" → "raw_file"
     # "Q.Value" → "q_value"
     # "R.FileName" → "r_filename"
     ```

### Supported Software Tools

This project aims to support outputs from major proteomics quantification software, including:

**DIA (Data-Independent Acquisition)**:
- DIA-NN
- Spectronaut
- AlphaDIA
- FragPipe (DIA-NN quant)
- MSAID
- MaxDIA
- PEAKS

**DDA (Data-Dependent Acquisition)**:
- MaxQuant
- AlphaPept
- Sage
- FragPipe/MSFragger
- WOMBAT
- ProlineStudio
- MSAngel
- i2MassChroQ
- quantms
- MetaMorpheus

## Key Resources

### Specification Reference
- **Main Specification**: `/Users/wolski/projects/anndata_omics_bridge/docs/AnnData_Omics_Bridge_spec.qmd`
- **Mapping Guide**: `docs/DIANN_to_AnnData_mapping.md` (example for DIANN)

### Reference Implementation - ProteoBench
This project should base its implementation on the [ProteoBench](https://github.com/Proteobench/ProteoBench) project structure:

**Example Data**:
- Location: `/Users/wolski/projects/ProteoBench/test/data/quant`
- Contains test data for multiple software tools (DDA and DIA)
- **Use these subdirectories for initial development**:
  - **`quant_lfq_ion_DDA_QExactive/`** - DDA ion-level test data ✅
  - **`quant_lfq_ion_DIA_AIF/`** - DIA ion-level test data ✅
  - `quant_lfq_peptidoform_DDA/` - Peptidoform data (NOT for initial scope) ❌

**Parser Implementation**:
- Location: `/Users/wolski/projects/ProteoBench/proteobench/io`
- Key modules:
  - **`parsing/parse_ion.py`** - **PRIMARY REFERENCE**: Ion/precursor level parsers for different software formats
  - `parsing/parse_peptidoform.py` - Peptidoform parsers (NOT in initial scope)
  - `params/` - Parameter extraction for each software tool
    - Each software has its own param file (e.g., `diann.py`, `spectronaut.py`, `maxquant.py`)
- Design pattern:
  - Format-specific load functions in `parse_ion.py` (e.g., `_load_diann()`, `_load_maxquant()`)
  - `_LOAD_FUNCTIONS` dictionary maps software names to loader functions
  - Each loader standardizes output to a common DataFrame format with ion/precursor level data

**Key Insights from ProteoBench**:
1. Use a registry pattern (`_LOAD_FUNCTIONS`) to map software names to loaders
2. Each software format has its own parameter extraction module
3. Standardize column names across different software outputs
4. Handle software-specific quirks (e.g., MaxQuant fixed modifications, AlphaDIA secondary files)

### Dependencies
- **omicsbridge**: Validation and column resolution utilities
  - Location: `/Users/wolski/projects/anndata_omics_bridge`
  - Import: `from omicsbridge import validate_anndata_omics, ColumnResolver`

## Project Structure

```
proteomics_to_anndata/
├── docs/
│   └── software_mappings/            # Column mapping references for each software
│       ├── DIANN_mapping.md
│       ├── Spectronaut_mapping.md
│       └── MaxQuant_mapping.md
├── src/proteomics_converter/
│   ├── __init__.py                   # Package exports
│   ├── converter.py                  # Main conversion logic
│   ├── parsers/                      # Software-specific parsers
│   │   ├── __init__.py
│   │   ├── base.py                   # Base parser interface
│   │   ├── diann.py                  # DIA-NN parser
│   │   ├── spectronaut.py            # Spectronaut parser
│   │   ├── maxquant.py               # MaxQuant parser
│   │   └── ...                       # Other software parsers
│   ├── annotation.py                 # Sample annotation handling
│   └── qc.py                         # Quality control and filtering
├── tests/
│   ├── test_converter.py             # Unit tests
│   └── test_parsers/                 # Parser-specific tests
│       ├── test_diann.py
│       └── ...
├── examples/
│   ├── sample_annotation.csv         # Example annotation file
│   ├── convert_diann.py              # DIANN example
│   ├── convert_spectronaut.py        # Spectronaut example
│   └── convert_maxquant.py           # MaxQuant example
├── pyproject.toml                    # Package configuration
├── README.md                         # User documentation
└── CLAUDE.md                         # This file (AI assistant instructions)
```

## AnnData Requirements (exploreDE-compatible)

### Required Components

**var (features/ions/precursors)**:
- **Comprehensive ion/precursor metadata** - Extract ALL available information from software output
- Minimum required columns:
  - `description`: Peptide sequence, modifications, charge (REQUIRED for exploreDE)
  - `label`: Precursor ID, modified sequence (optional)
- **All other metadata columns from source** (normalized names):
  - Extract every non-quantification column from the software output
  - Normalize column names (lowercase, underscores, no special chars)
  - Preserve original column semantics - don't rename to different meanings
  - Examples of what you'll find (software-dependent):
    - DIA-NN: `modified_sequence`, `stripped_sequence`, `precursor_charge`, `precursor_mz`, `rt`, `q_value`, `pg_q_value`, `protein_ids`, `protein_names`, `genes`
    - MaxQuant: `modified_sequence`, `charge`, `mass`, `retention_time`, `pep`, `score`, `proteins`, `leading_proteins`, `gene_names`
    - Spectronaut: `eg_modifiedsequence`, `fg_charge`, `fg_precursormz`, `fg_rt`, `eg_qvalue`, `pg_proteingroups`, `pg_proteinnames`
- **Philosophy**: Keep all non-quantification columns in `var` to preserve maximum information with original naming

**obs (samples)**:
- `factor`: Experimental factors like condition, batch (REQUIRED)
- `label`: Sample identifiers (optional)

**X (abundance matrix)**:
- Primary quantification values (typically log2-transformed)
- Shape: (n_samples × n_precursors)
- Choose the most appropriate quantification type as default (e.g., normalized intensity)

**layers** (ALL available quantification types):
- Extract **every** quantification column from software output
- Store each quantification type in a separate layer
- Use descriptive naming convention
- Examples by software:
  - **DIA-NN**: `ms1_area`, `precursor_quantity`, `precursor_normalized`, `fragment_quant_corrected`
  - **MaxQuant**: `intensity`, `area`, `intensity_normalized`
  - **Spectronaut**: `fg_quantity`, `eg_quantity` (fragment/elution group), multiple normalization states
  - **AlphaDIA**: `ms1_height`, `ms1_area`, `ms2_area`
- Keep raw (non-log-transformed) values in layers for flexibility

**uns (metadata)**:
```python
adata.uns['exploreDE'] = {
    'column_roles': {
        'var': {
            'description': ['modified_sequence'],  # or 'peptide_sequence_charge'
            'label': ['precursor_id', 'protein_id', 'gene_symbol']
        },
        'obs': {
            'factor': ['condition', 'batch'],
            'label': ['run', 'file_name']
        }
    },
    'quantification_level': 'ion',  # Specify this is ion/precursor level
    'software': 'DIA-NN',  # Source software
    'software_version': '1.9.0',
    'layers': {
        # Document what each layer contains
        'ms1_area': 'Raw MS1 area values',
        'precursor_quantity': 'Precursor quantity (normalized)',
        'precursor_translated': 'Library-free precursor quantity',
        'fragment_quant_corrected': 'Fragment ion quantities (corrected)'
    },
    'primary_quantification': 'precursor_quantity'  # Which layer was used for X
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

# Run software-specific examples
python examples/convert_diann.py
python examples/convert_spectronaut.py
python examples/convert_maxquant.py
```

**Running tests:**
```bash
source .venv/bin/activate
pytest tests/
```

## Development Principles

1. **Self-contained code blocks**: Each example should import all necessary dependencies
2. **Validation**: Always validate output with `omicsbridge.validate_anndata_omics()`
3. **Comprehensive data extraction**:
   - Extract ALL metadata columns → store in `var`
   - Extract ALL quantification columns → store in separate `layers`
   - Preserve maximum information from source software
   - Document extracted columns in `uns['layers']`
4. **Error handling**: Provide clear error messages for common issues:
   - Missing annotation file
   - Sample ID mismatches
   - Invalid/unrecognized software format
   - Unsupported software version
   - Missing required columns in software output
5. **Flexible column mapping**: Allow users to specify custom column mappings
6. **Modular design**: Follow ProteoBench pattern with registry-based parsers

## Common Tasks

### Reading Software Output Files (Ion/Precursor Level)
- Support multiple file formats (parquet, TSV, CSV, Excel)
- Handle different software versions (column name variations)
- **Extract ion/precursor-level quantification data**:
  - DIA-NN: `report.tsv` (precursor-level)
  - MaxQuant: `evidence.txt` (precursor-level)
  - Spectronaut: precursor-level exports
  - etc.
- Auto-detect software format when possible

### Parser Development
- Follow ProteoBench design pattern with important enhancements:
  1. Create format-specific loader function (e.g., `_load_diann()`)
  2. Add to parser registry dictionary
  3. **Extract ALL columns from software output** (not just minimal set)
  4. **Normalize column names** (minimal transformation):
     - Convert to lowercase
     - Replace whitespace with underscores
     - Replace dots/special characters with underscores
     - Ensure valid Python identifiers
     - **DO NOT rename** to different semantic meanings
  5. Separate quantification columns from metadata:
     - Identify all quantification columns (intensity, area, quantity, etc.)
     - Keep all other columns as metadata for `var`
  6. Handle format-specific edge cases
- Parser output structure:
  ```python
  {
      'data': pd.DataFrame,  # Long format with normalized column names
      'quant_columns': List[str],  # Names of quantification columns to extract as layers
      'metadata_columns': List[str],  # Names of columns to include in var
      'sample_column': str,  # Column identifying samples (normalized name)
      'feature_column': str  # Column identifying precursors/ions (normalized name)
  }
  ```

### Sample Annotation
- Match software-specific run identifiers to annotation file (using normalized column names):
  - DIA-NN: `run` or `file_name` (normalized from `Run` or `File.Name`)
  - MaxQuant: `raw_file` (normalized from `Raw file`)
  - Spectronaut: `r_filename` (normalized from `R.FileName`)
  - etc.
- Handle filename variations (paths, extensions)
- Validate that all samples have annotations
- Preserve original column names in documentation/metadata

### Data Reshaping
- Pivot long format → wide format (samples × ions/precursors)
- **Create separate wide-format matrices for each quantification type**:
  1. Identify all quantification columns in source data
  2. Pivot each quantification column separately
  3. Store each pivoted matrix as a layer in AnnData
  4. Choose primary quantification for `X` (document choice in `uns`)
- Handle missing values appropriately
  - Preserve NaN/missing patterns from source data
  - Don't impute unless explicitly requested
- Extract and aggregate metadata columns into `var`:
  - Metadata should be consistent across samples (same for each precursor)
  - If conflicts exist, choose most informative value or concatenate
- Column names will already be normalized by parser (lowercase, underscores)
- **Do not rename** columns across software - preserve original semantics
- **Note**: Protein-level aggregation is NOT in initial scope

**Example reshaping workflow**:
```python
# Identify quantification columns
quant_cols = ['Ms1.Area', 'Precursor.Quantity', 'Precursor.Normalised']

# Create layers dictionary
layers = {}
for qcol in quant_cols:
    # Pivot each quantification type
    wide = df.pivot(index='Run', columns='precursor_id', values=qcol)
    layers[qcol.lower().replace('.', '_')] = wide

# Extract metadata (consistent per precursor)
var = df.groupby('precursor_id').first()[metadata_columns]

# Build AnnData
adata = AnnData(X=layers['precursor_quantity'], layers=layers, var=var)
```

### Quality Control
- Filter by protein/peptide Q-value or FDR
- Remove low-quality features
- Handle missing values (filtering vs imputation)
- Software-specific quality metrics

## Testing Strategy

1. **Unit tests**: Test individual functions (reshape, filter, validate)
2. **Parser tests**: Test each software-specific parser with real example data
   - **Use ion-level test data from ProteoBench**:
     - `/Users/wolski/projects/ProteoBench/test/data/quant/quant_lfq_ion_DDA_QExactive/`
     - `/Users/wolski/projects/ProteoBench/test/data/quant/quant_lfq_ion_DIA_AIF/`
   - Test both DDA and DIA formats
   - Verify correct extraction of ion/precursor level data
   - **Verify comprehensive extraction**:
     - Count columns in source file vs. columns in `var`
     - Ensure no metadata columns are lost
     - Verify all quantification columns are captured in layers
3. **Integration tests**: Test full conversion pipeline for each supported software
   - Verify output AnnData has correct shape (n_samples × n_precursors)
   - Confirm `uns['quantification_level'] == 'ion'`
   - **Test layer completeness**:
     - Verify number of layers matches expected quantification columns
     - Check layer names are descriptive and documented in `uns['layers']`
     - Ensure `uns['primary_quantification']` is set correctly
   - **Test var completeness**:
     - Verify all metadata columns from source are present
     - Check column names are normalized (lowercase, underscores, no special chars)
     - Verify column names preserve original semantics (not renamed to different meanings)
     - Ensure software-specific columns are preserved with their original names (normalized)
4. **Validation**: All outputs must pass `validate_anndata_omics()`
5. **Edge cases**: Test software version variations, missing values, malformed files
6. **Data integrity**: Verify no information loss during conversion
   - Compare precursor counts: source file vs. AnnData
   - Spot-check quantification values match between source and layers
   - Verify metadata values match between source and var

## Code Style

- Follow PEP 8
- Use type hints
- Document functions with docstrings
- Keep functions focused and modular

## Version History

- v0.1.0: Initial project setup
