# DIANN to AnnData Converter

Convert DIANN (DIA-NN) proteomics output to AnnData format for downstream analysis with exploreDE, prolfqua, and other tools.

## Features

- ✅ Convert DIANN parquet/TSV files to AnnData `.h5ad` format
- ✅ Support protein-level, peptide-level, and precursor-level data
- ✅ Follows AnnData Omics Bridge specification
- ✅ Compatible with exploreDE and prolfqua
- ✅ Built-in quality control and filtering
- ✅ Automatic validation of output files

## Installation

```bash
# Clone the repository
cd ~/projects/diann_to_anndata

# Install with uv (recommended)
uv venv
source .venv/bin/activate
uv pip install -e .

# Or with pip
pip install -e .
```

## Quick Start

```python
from diann_converter import convert_diann

# Convert DIANN output to AnnData
adata = convert_diann(
    diann_file='diann_report.parquet',
    annotation_file='samples.csv',
    output_file='proteomics.h5ad',
    level='protein'  # 'protein', 'peptide', or 'precursor'
)
```

## Sample Annotation File

Create a CSV file with your experimental design:

```csv
sample_id,condition,batch,replicate
sample_001,control,batch1,1
sample_002,control,batch1,2
sample_003,treated,batch1,1
sample_004,treated,batch1,2
```

**Important**: The `sample_id` column must match DIANN's `Run` or `File.Name` column.

## Documentation

- **Mapping Guide**: See `docs/DIANN_to_AnnData_mapping.md` for detailed column mapping
- **Specification**: Based on `/Users/wolski/projects/anndata_omics_bridge/docs/AnnData_Omics_Bridge_spec.qmd`

## Requirements

- Python ≥ 3.9
- pandas
- numpy
- anndata
- pyarrow (for parquet support)
- omicsbridge (for validation)

## License

MIT License
