"""
Example: Convert DIA-NN output to AnnData format.

This example uses test data from ProteoBench to demonstrate the conversion workflow.
"""

import sys
from pathlib import Path

# Add src to path
src_path = Path(__file__).parent.parent / "src"
sys.path.insert(0, str(src_path))

from diann_converter import proteomics_to_anndata

# Paths
DATA_PATH = "/Users/wolski/projects/ProteoBench/test/data/quant/quant_lfq_ion_DIA_AIF/DIANN_1.9_beta_sample_report.tsv"
ANNOTATION_PATH = Path(__file__).parent / "diann_annotation.csv"
OUTPUT_PATH = Path(__file__).parent / "diann_example_output.h5ad"

def main():
    """Convert DIA-NN data to AnnData."""
    print("=" * 60)
    print("DIA-NN to AnnData Conversion Example")
    print("=" * 60)

    print(f"\nInput data: {DATA_PATH}")
    print(f"Annotation: {ANNOTATION_PATH}")
    print(f"Output: {OUTPUT_PATH}\n")

    # Convert to AnnData
    adata = proteomics_to_anndata(
        data_path=DATA_PATH,
        annotation_path=str(ANNOTATION_PATH),
        software="diann",
        annotation_id_col="sample",
        factor_cols=["condition", "batch"],
        label_cols=["replicate"],
        qvalue_threshold=0.01,
        log2_transform=True,
    )

    print("\n" + "=" * 60)
    print("Conversion Complete!")
    print("=" * 60)

    print(f"\nAnnData object shape: {adata.shape} (samples × precursors)")
    print(f"Number of samples: {adata.n_obs}")
    print(f"Number of precursors: {adata.n_vars}")

    print("\nSample metadata (obs):")
    print(adata.obs)

    print("\nPrecursor metadata (var) - first 5:")
    print(adata.var.head())

    print("\nLayers:")
    print(f"  - X: Log2-transformed intensities, shape {adata.X.shape}")
    print(f"  - intensities: Raw intensities, shape {adata.layers['intensities'].shape}")

    print("\nexploreDE metadata:")
    print(f"  Factors: {adata.uns['exploreDE']['column_roles']['obs']['factor']}")
    print(f"  Labels: {adata.uns['exploreDE']['column_roles']['obs']['label']}")

    # Save to file
    adata.write_h5ad(OUTPUT_PATH)
    print(f"\nSaved AnnData object to: {OUTPUT_PATH}")

    # Optional: Validate with omicsbridge (if available)
    try:
        from omicsbridge import validate_anndata_omics
        print("\nValidating with omicsbridge...")
        is_valid = validate_anndata_omics(adata)
        print(f"Validation result: {'PASS' if is_valid else 'FAIL'}")
    except ImportError:
        print("\nomicsbridge not available - skipping validation")

    print("\n" + "=" * 60)


if __name__ == "__main__":
    main()
