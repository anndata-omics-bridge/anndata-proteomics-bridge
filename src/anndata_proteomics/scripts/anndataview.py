#!/usr/bin/env python3
"""
AnnData Viewer - A marimo-based interactive viewer for h5ad files.

Usage:
    marimo run src/anndata_proteomics/scripts/anndataview.py -- path/to/file.h5ad
    marimo edit src/anndata_proteomics/scripts/anndataview.py -- path/to/file.h5ad
"""

import marimo

__generated_with = "0.18.4"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    import sys
    import os
    return mo, os, sys


@app.cell
def _(mo, sys):
    # Get file path from CLI arguments
    cli_args = mo.cli_args()

    if cli_args and len(cli_args) > 0:
        file_path = cli_args[0]
    elif len(sys.argv) > 1:
        # Fallback for direct script execution
        file_path = sys.argv[-1]
    else:
        file_path = None

    # Create file browser for selecting h5ad files
    file_input = mo.ui.text(
        value=file_path or "",
        placeholder="Enter path to .h5ad file",
        label="H5AD File Path",
        full_width=True
    )
    return cli_args, file_input, file_path


@app.cell
def _(file_input, mo):
    mo.md(f"""
    # AnnData Viewer

    Interactive viewer for AnnData (.h5ad) files.

    {file_input}
    """)
    return


@app.cell
def _(file_input, mo, os):
    import anndata as ad

    current_path = file_input.value.strip()

    if not current_path:
        mo.stop(True, mo.md("**Please enter a path to an h5ad file above.**"))

    if not os.path.exists(current_path):
        mo.stop(True, mo.md(f"**File not found:** `{current_path}`"))

    if not current_path.endswith('.h5ad'):
        mo.stop(True, mo.md(f"**Not an h5ad file:** `{current_path}`"))

    # Load the AnnData object
    adata = ad.read_h5ad(current_path)

    mo.md(f"""
    ## File: `{os.path.basename(current_path)}`

    **Path:** `{current_path}`
    """)
    return ad, adata, current_path


@app.cell
def _(adata, mo):
    # Overview section
    n_obs, n_vars = adata.shape

    # Build structure info
    layers_list = list(adata.layers.keys()) if adata.layers else []
    obsm_list = list(adata.obsm.keys()) if adata.obsm else []
    varm_list = list(adata.varm.keys()) if adata.varm else []
    uns_keys = list(adata.uns.keys()) if adata.uns else []

    overview_md = f"""
    ## Overview

    | Property | Value |
    |----------|-------|
    | **Shape** | {n_obs} samples × {n_vars} features |
    | **obs columns** | {len(adata.obs.columns)} |
    | **var columns** | {len(adata.var.columns)} |
    | **Layers** | {', '.join(layers_list) if layers_list else 'None'} |
    | **obsm** | {', '.join(obsm_list) if obsm_list else 'None'} |
    | **varm** | {', '.join(varm_list) if varm_list else 'None'} |
    | **uns keys** | {', '.join(uns_keys) if uns_keys else 'None'} |
    """

    mo.md(overview_md)
    return layers_list, n_obs, n_vars, obsm_list, overview_md, uns_keys, varm_list


@app.cell
def _(mo):
    # Create tabs for different sections
    tab_selector = mo.ui.tabs({
        "obs": "Sample Metadata (obs)",
        "var": "Feature Metadata (var)",
        "X": "Expression Matrix (X)",
        "layers": "Layers",
        "uns": "Unstructured (uns)",
    })
    tab_selector
    return (tab_selector,)


@app.cell
def _(adata, mo, tab_selector):
    import pandas as pd
    import numpy as np

    selected_tab = tab_selector.value

    if selected_tab == "obs":
        # Sample metadata
        obs_df = adata.obs.copy()
        obs_df.insert(0, '_index', obs_df.index)

        content = mo.vstack([
            mo.md(f"### Sample Metadata ({len(adata.obs)} samples, {len(adata.obs.columns)} columns)"),
            mo.md(f"**Columns:** {', '.join(adata.obs.columns.tolist())}"),
            mo.ui.table(obs_df, selection=None, page_size=20)
        ])

    elif selected_tab == "var":
        # Feature metadata
        var_df = adata.var.copy()
        var_df.insert(0, '_index', var_df.index)

        # Add search/filter
        content = mo.vstack([
            mo.md(f"### Feature Metadata ({len(adata.var)} features, {len(adata.var.columns)} columns)"),
            mo.md(f"**Columns:** {', '.join(adata.var.columns.tolist())}"),
            mo.ui.table(var_df, selection=None, page_size=20)
        ])

    elif selected_tab == "X":
        # Expression matrix preview
        X_preview = pd.DataFrame(
            adata.X[:min(10, adata.n_obs), :min(20, adata.n_vars)],
            index=adata.obs_names[:min(10, adata.n_obs)],
            columns=adata.var_names[:min(20, adata.n_vars)]
        )

        # Basic statistics
        X_flat = adata.X.flatten()
        X_flat_valid = X_flat[~np.isnan(X_flat)]

        stats_md = f"""
        ### Expression Matrix (X)

        **Shape:** {adata.X.shape[0]} × {adata.X.shape[1]}

        | Statistic | Value |
        |-----------|-------|
        | Min | {np.min(X_flat_valid):.4f} |
        | Max | {np.max(X_flat_valid):.4f} |
        | Mean | {np.mean(X_flat_valid):.4f} |
        | Median | {np.median(X_flat_valid):.4f} |
        | Std | {np.std(X_flat_valid):.4f} |
        | NaN count | {np.sum(np.isnan(X_flat))} ({100*np.sum(np.isnan(X_flat))/len(X_flat):.1f}%) |

        **Preview (first 10 samples × first 20 features):**
        """

        content = mo.vstack([
            mo.md(stats_md),
            mo.ui.table(X_preview.round(4), selection=None)
        ])

    elif selected_tab == "layers":
        # Layers information
        if not adata.layers:
            content = mo.md("### Layers\n\nNo layers present in this AnnData object.")
        else:
            layers_info = []
            for layer_name, layer_data in adata.layers.items():
                layer_flat = layer_data.flatten()
                layer_valid = layer_flat[~np.isnan(layer_flat)]
                layers_info.append({
                    "Layer": layer_name,
                    "Shape": f"{layer_data.shape[0]} × {layer_data.shape[1]}",
                    "Min": f"{np.min(layer_valid):.4f}" if len(layer_valid) > 0 else "N/A",
                    "Max": f"{np.max(layer_valid):.4f}" if len(layer_valid) > 0 else "N/A",
                    "Mean": f"{np.mean(layer_valid):.4f}" if len(layer_valid) > 0 else "N/A",
                    "NaN %": f"{100*np.sum(np.isnan(layer_flat))/len(layer_flat):.1f}%"
                })

            content = mo.vstack([
                mo.md(f"### Layers ({len(adata.layers)} layers)"),
                mo.ui.table(pd.DataFrame(layers_info), selection=None)
            ])

    elif selected_tab == "uns":
        # Unstructured metadata
        if not adata.uns:
            content = mo.md("### Unstructured Metadata (uns)\n\nNo unstructured metadata present.")
        else:
            import json

            def format_uns_value(v, max_depth=3, current_depth=0):
                """Format uns values for display."""
                if current_depth >= max_depth:
                    return str(type(v).__name__)
                if isinstance(v, dict):
                    formatted = {k: format_uns_value(val, max_depth, current_depth+1)
                                for k, val in v.items()}
                    return formatted
                elif isinstance(v, (list, tuple)) and len(v) > 10:
                    return f"[{type(v).__name__} with {len(v)} items]"
                elif isinstance(v, np.ndarray):
                    return f"ndarray{v.shape}"
                elif isinstance(v, pd.DataFrame):
                    return f"DataFrame{v.shape}"
                else:
                    return v

            uns_formatted = {k: format_uns_value(v) for k, v in adata.uns.items()}

            try:
                uns_json = json.dumps(uns_formatted, indent=2, default=str)
            except (TypeError, ValueError):
                uns_json = str(uns_formatted)

            content = mo.vstack([
                mo.md("### Unstructured Metadata (uns)"),
                mo.md(f"**Keys:** {', '.join(adata.uns.keys())}"),
                mo.md(f"```json\n{uns_json}\n```")
            ])
    else:
        content = mo.md("Select a tab above to view data.")

    content
    return (
        X_flat,
        X_flat_valid,
        X_preview,
        content,
        np,
        obs_df,
        pd,
        selected_tab,
        stats_md,
        var_df,
    )


@app.cell
def _(adata, mo, np, pd):
    # Heatmap visualization (optional, shown below tabs)
    mo.md("---")

    show_heatmap = mo.ui.checkbox(label="Show expression heatmap (first 50 samples × 50 features)")
    show_heatmap
    return (show_heatmap,)


@app.cell
def _(adata, mo, np, pd, show_heatmap):
    import plotly.express as px

    heatmap_output = None

    if show_heatmap.value:
        # Subset for visualization
        n_samples = min(50, adata.n_obs)
        n_features = min(50, adata.n_vars)

        heatmap_data = pd.DataFrame(
            adata.X[:n_samples, :n_features],
            index=adata.obs_names[:n_samples],
            columns=adata.var_names[:n_features]
        )

        fig = px.imshow(
            heatmap_data,
            labels=dict(x="Features", y="Samples", color="Value"),
            aspect="auto",
            color_continuous_scale="RdBu_r"
        )
        fig.update_layout(height=600)

        heatmap_output = mo.ui.plotly(fig)

    heatmap_output
    return (heatmap_output,)


if __name__ == "__main__":
    app.run()
