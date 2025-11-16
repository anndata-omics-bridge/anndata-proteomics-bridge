# DIANN to AnnData Conversion Mapping Guide

## Overview

This guide maps DIANN (DIA-NN) parquet output to AnnData format following the **AnnData Omics Bridge Specification** for proteomics data compatible with **exploreDE** and **prolfqua**.

**Specification Location**: `/Users/wolski/projects/anndata_omics_bridge/docs/AnnData_Omics_Bridge_spec.qmd`

---

## DIANN Output Structure

DIANN produces parquet/TSV files with combined protein/peptide/precursor information across all samples. Typical structure:

```
File.Name | Run | Protein.Group | Protein.Ids | Protein.Names | Genes | PG.Quantity | PG.Normalised | ...
----------|-----|---------------|-------------|---------------|-------|-------------|---------------|----
sample1   | run1| 1             | Q9Y6K9      | NF-kB p65     | RELA  | 1234567     | 0.89          | ...
sample1   | run1| 2             | P04637      | p53 protein   | TP53  | 9876543     | 1.12          | ...
sample2   | run2| 1             | Q9Y6K9      | NF-kB p65     | RELA  | 1345678     | 0.91          | ...
```

**Key characteristic**: Long format (one row per protein per sample) → needs reshaping to wide format for AnnData

---

## Conversion Strategy

### Step 1: Identify Analysis Level

DIANN supports multiple aggregation levels:
- **Protein level** (`PG.Quantity`, `PG.Normalised`) ← **Recommended for exploreDE**
- **Peptide level** (`Precursor.Quantity`, `Precursor.Normalised`)
- **Precursor level** (for PTM analysis)

**Choose**: Protein-level for standard differential expression analysis.

### Step 2: Reshape Data

Transform from long format to wide format:
- **Rows** = samples (unique `Run` or `File.Name`)
- **Columns** = proteins (unique `Protein.Group` or `Protein.Ids`)
- **Values** = quantification (`PG.Quantity`, `PG.Normalised`)

---

## Column Mapping: DIANN → AnnData

### var (Feature Annotations)

| AnnData Role | DIANN Column(s) | Required? | Notes |
|--------------|----------------|-----------|-------|
| **Index** | `Protein.Group` or `Protein.Ids` | Yes | Unique protein identifiers |
| `description` | `Protein.Names` | **REQUIRED** | Free-text searchable (exploreDE requirement) |
| `label` | `Genes`, `Protein.Ids`, `First.Protein.Description` | No | Human-readable identifiers for display/filtering |
| (custom) | `PG.Q.Value`, `Global.PG.Q.Value` | No | Protein identification quality scores |
| (custom) | Mean intensity, # peptides, etc. | No | Computed summary statistics |

**Example**:
```python
var = pd.DataFrame({
    'protein_id': ['Q9Y6K9', 'P04637', ...],          # Index
    'description': ['NF-kB p65', 'p53 protein', ...], # REQUIRED
    'gene_symbol': ['RELA', 'TP53', ...],             # label
    'protein_names': ['NF-kB p65', 'p53 protein', ...], # label
    'global_qvalue': [0.001, 0.002, ...],             # custom
    'n_peptides': [12, 8, ...],                        # custom
    'mean_intensity': [1.2e6, 8.9e5, ...]              # custom
}, index=['Q9Y6K9', 'P04637', ...])
```

### obs (Sample Annotations)

| AnnData Role | DIANN Column(s) | Required? | Notes |
|--------------|----------------|-----------|-------|
| **Index** | `Run` or `File.Name` | Yes | Unique sample identifiers |
| `factor` | **External annotation file** | **REQUIRED** | Experimental design (condition, batch, etc.) |
| `label` | `Run`, `File.Name` | No | Sample identifiers for display |

**Important**: DIANN output typically **does NOT include experimental design**. You must provide a separate annotation file:

**Example annotation file** (`samples.csv`):
```csv
sample_id,condition,batch,replicate,file_name
sample_001,control,batch1,1,20231015_sample_001.raw
sample_002,control,batch1,2,20231015_sample_002.raw
sample_003,treated,batch1,1,20231015_sample_003.raw
sample_004,treated,batch1,2,20231015_sample_004.raw
```

**Example**:
```python
obs = pd.DataFrame({
    'run': ['sample_001', 'sample_002', ...],         # Index
    'condition': ['control', 'control', 'treated', ...], # factor (REQUIRED)
    'batch': ['batch1', 'batch1', 'batch1', ...],     # factor
    'replicate': [1, 2, 1, ...],                      # label
    'file_name': ['20231015_sample_001.raw', ...]     # label
}, index=['sample_001', 'sample_002', ...])
```

### X and layers (Abundance Matrices)

| Layer Name | DIANN Column | Transformation | Purpose |
|------------|-------------|----------------|---------|
| `intensities` | `PG.Quantity` | None (raw) | Raw protein intensities |
| `normalised` | `PG.Normalised` | DIANN internal normalization | DIANN-normalized intensities |
| `X` (default) | `PG.Quantity` | **log2(x + 1)** | **Analysis-ready** transformed data |

**Convention**: Store **log2-transformed** data in `X` (most DE methods expect this).

**Example**:
```python
# Reshape DIANN long format to wide matrix (samples × proteins)
intensity_matrix = diann_df.pivot(
    index='Run',           # Samples as rows
    columns='Protein.Ids', # Proteins as columns
    values='PG.Quantity'   # Quantification values
)

# Store in AnnData
adata.layers['intensities'] = intensity_matrix.values  # Raw
adata.X = np.log2(intensity_matrix.values + 1)          # Log2 transformed
adata.uns['X_layer_name'] = 'log2_intensity'            # Document transformation
```

### varm (DE Results)

If DIANN includes statistical testing results (or you add them later):

| AnnData Role | DIANN Column(s) | Required? | Notes |
|--------------|----------------|-----------|-------|
| `effect` | `Log2.Fold.Change`, `Effect.Size` | **REQUIRED** | Effect size for DE |
| `score` | `P.Value`, `Q.Value`, `Adjusted.P.Value` | **REQUIRED** | Significance scores |
| `label` | Other stats | No | Additional statistics |

**Note**: DIANN's main report typically doesn't include DE results. These are added by downstream analysis (DEqMS, limma, etc.).

---

## Metadata Structure (uns)

### exploreDE-compatible metadata:

```python
adata.uns['exploreDE'] = {
    'column_roles': {
        'var': {
            'description': ['description'],  # REQUIRED
            'label': ['gene_symbol', 'protein_names', 'protein_id', 'global_qvalue', 'n_peptides']
        },
        'obs': {
            'factor': ['condition', 'batch'],  # REQUIRED (at least one)
            'label': ['run', 'file_name', 'replicate']
        }
    }
}

# DE results added later when performing statistical analysis
# adata.uns['exploreDE']['de_tests'] = {
#     'DE_treated_vs_control': {
#         'layer_used': 'X',
#         'factor_used': ['condition'],
#         'contrast_formula': 'treated - control',
#         'model': 'limma'
#     }
# }
```

### prolfqua-compatible metadata (if needed):

```python
adata.uns['prolfqua'] = {
    'hierarchy': ['protein_id'],  # For protein-level data
    'column_roles': {
        'var': {
            'intensity': ['mean_intensity'],
            'qvalue': ['global_qvalue'],
            'label': ['gene_symbol', 'protein_names']
        },
        'obs': {
            'sample_id': ['run']
        }
    }
}
```

---

## Conversion Workflow Pseudocode

```python
import pandas as pd
import numpy as np
import anndata as ad

# 1. Load DIANN parquet output
diann_df = pd.read_parquet('diann_report.parquet')

# 2. Load sample annotation (experimental design)
sample_annot = pd.read_csv('samples.csv', index_col='sample_id')

# 3. Filter to protein-level quantification
protein_df = diann_df[['Run', 'Protein.Ids', 'Protein.Names', 'Genes',
                        'PG.Quantity', 'PG.Normalised', 'Global.PG.Q.Value']].copy()

# 4. Reshape to wide format (samples × proteins)
intensity_wide = protein_df.pivot(
    index='Run',
    columns='Protein.Ids',
    values='PG.Quantity'
).fillna(0)  # Handle missing values

# 5. Create var (protein annotations) - one row per protein
var_df = protein_df.groupby('Protein.Ids').first()[
    ['Protein.Names', 'Genes', 'Global.PG.Q.Value']
].rename(columns={
    'Protein.Names': 'description',    # REQUIRED role
    'Genes': 'gene_symbol',             # label role
    'Global.PG.Q.Value': 'global_qvalue'
})
var_df.index.name = 'protein_id'

# 6. Create obs (sample annotations)
obs_df = sample_annot.loc[intensity_wide.index].copy()

# 7. Build AnnData
adata = ad.AnnData(
    X=np.log2(intensity_wide.values + 1),  # Log2 transform
    obs=obs_df,
    var=var_df
)

# 8. Store raw intensities in layer
adata.layers['intensities'] = intensity_wide.values

# 9. Add metadata
adata.uns['X_layer_name'] = 'log2_intensity'
adata.uns['exploreDE'] = {
    'column_roles': {
        'var': {
            'description': ['description'],
            'label': ['gene_symbol', 'global_qvalue']
        },
        'obs': {
            'factor': ['condition', 'batch'],
            'label': ['run', 'file_name']
        }
    }
}

# 10. Save
adata.write('diann_proteomics.h5ad')
```

---

## Data Quality Considerations

### Missing Values
- **DIANN strategy**: Missing values often indicate protein not detected
- **Options**:
  - Set to 0 (treated as absent)
  - Use small imputation value (e.g., min_intensity / 2)
  - Keep as NaN and filter proteins with too many missing values

### Filtering Recommendations
1. **Protein quality**: Filter by `Global.PG.Q.Value` < 0.01
2. **Detection frequency**: Keep proteins detected in ≥50% of samples in at least one condition
3. **Intensity threshold**: Remove very low abundance proteins (background noise)

### Sample Matching
- **Critical**: Ensure DIANN `Run` or `File.Name` matches your annotation file
- Handle filename variations (`.raw`, `.mzML`, path prefixes)

---

## Validation

After conversion, validate with:

```python
from omicsbridge import validate_anndata_omics

result = validate_anndata_omics(adata, app_name='exploreDE')
print(f"Status: {result['overall_status']}")

if result['errors']:
    print("Errors:")
    for error in result['errors']:
        print(f"  - {error}")
```

---

## Next Steps

1. **Create converter script**: Implement the conversion workflow
2. **Test with sample data**: Use small DIANN output for testing
3. **Add DE analysis**: Integrate limma/DEqMS for differential expression
4. **Visualize in exploreDE**: Load `.h5ad` file for interactive exploration

---

## References

- **AnnData Specification**: `/Users/wolski/projects/anndata_omics_bridge/docs/AnnData_Omics_Bridge_spec.qmd`
- **DIANN Documentation**: https://github.com/vdemichev/DiaNN
- **Validator**: `src/omicsbridge/exploreDE_validator.py`
