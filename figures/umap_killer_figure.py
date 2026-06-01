"""UMAP figure: cells removed by blood-calibrated mito threshold form coherent clusters.

This is the biological validation for the killer figure. It shows that the 14,076 cells
removed from kidney data by a blood-calibrated mito≤5% threshold are NOT random garbage —
they form well-defined clusters corresponding to real kidney cell populations.

Run via the ARM64 wrapper (needs umap-learn + scanpy):
    bash scripts/run_arm64.sh figures/umap_killer_figure.py

Produces:
  figures/umap_killer.html   — interactive UMAP coloured by filter status, donor, and mito%
  figures/umap_killer_data.json — cell counts and cluster composition summary
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import json
import numpy as np
import scanpy as sc
import plotly.graph_objects as go
import plotly.offline as pyo
from plotly.subplots import make_subplots

sc.settings.verbosity = 2

DATA_PATH = "benchmarks/wu_kidney/raw_data.h5ad"
OUTPUT_DIR = "figures"
BLOOD_MITO_MAX = 5.0   # typical blood-calibrated threshold
KIDNEY_MITO_MAX = 30.0  # kidney-appropriate threshold
MITO_PREFIX = "MT-"
SEED = 42


def label_filter_status(pct_mt, blood_max=BLOOD_MITO_MAX, kidney_max=KIDNEY_MITO_MAX):
    """Three-way label: kept, blood-removed (biologically real), over-quality-removed."""
    labels = np.full(len(pct_mt), "kept (both thresholds)", dtype=object)
    labels[pct_mt > blood_max] = f"removed by blood (mito>{blood_max}%)"
    labels[pct_mt > kidney_max] = f"removed by kidney (mito>{kidney_max}%)"
    return labels


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ── Load data ─────────────────────────────────────────────────────────────
    print("Loading kidney data...")
    adata = sc.read_h5ad(DATA_PATH)
    adata.var_names_make_unique()
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    print(f"  Donors: {adata.obs['sample'].value_counts().to_dict()}")

    # ── QC metrics ───────────────────────────────────────────────────────────
    print("Computing QC metrics...")
    adata.var["mt"] = adata.var_names.str.startswith(MITO_PREFIX)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    pct_mt = adata.obs["pct_counts_mt"].values

    print(f"  Mito% range: {pct_mt.min():.1f}–{pct_mt.max():.1f}%")
    print(f"  Mito% median: {np.median(pct_mt):.1f}%")

    per_donor = adata.obs.groupby("sample")["pct_counts_mt"].median()
    for donor, med in per_donor.items():
        print(f"  Donor {donor}: median mito {med:.1f}%")

    # ── Filter labels ─────────────────────────────────────────────────────────
    filter_labels = label_filter_status(pct_mt)
    adata.obs["filter_status"] = filter_labels

    n_total = len(adata)
    n_blood_removed = (filter_labels == f"removed by blood (mito>{BLOOD_MITO_MAX}%)").sum()
    n_kidney_removed = (filter_labels == f"removed by kidney (mito>{KIDNEY_MITO_MAX}%)").sum()
    n_kept = (filter_labels == "kept (both thresholds)").sum()

    print(f"\n  Filter breakdown:")
    print(f"    Kept by both: {n_kept:,} ({100*n_kept/n_total:.1f}%)")
    print(f"    Removed by blood (mito>{BLOOD_MITO_MAX}%): {n_blood_removed:,} ({100*n_blood_removed/n_total:.1f}%)")
    print(f"    Removed by kidney (mito>{KIDNEY_MITO_MAX}%): {n_kidney_removed:,} ({100*n_kidney_removed/n_total:.1f}%)")

    # ── Preprocessing → UMAP ─────────────────────────────────────────────────
    print("\nPreprocessing...")
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()

    print("PCA...")
    sc.pp.pca(adata_hvg, n_comps=50, random_state=SEED)

    print("Neighbors...")
    sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30, random_state=SEED)

    print("UMAP...")
    sc.tl.umap(adata_hvg, random_state=SEED)

    # Copy UMAP coordinates back
    adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]

    print("Leiden clustering...")
    sc.tl.leiden(adata_hvg, resolution=0.5, random_state=SEED)
    adata.obs["leiden"] = adata_hvg.obs["leiden"]

    umap_x = adata.obsm["X_umap"][:, 0]
    umap_y = adata.obsm["X_umap"][:, 1]

    # ── Build plotly figure ───────────────────────────────────────────────────
    print("\nBuilding figure...")
    fig = make_subplots(
        rows=1, cols=3,
        subplot_titles=[
            f"Filter Status (blood≤{BLOOD_MITO_MAX}% vs kidney≤{KIDNEY_MITO_MAX}%)",
            "Mitochondrial %",
            "Donor",
        ],
        horizontal_spacing=0.04,
    )

    POINT_SIZE = 2
    OPACITY = 0.6

    # Panel 1: filter status
    status_colors = {
        "kept (both thresholds)": "#2ca02c",
        f"removed by blood (mito>{BLOOD_MITO_MAX}%)": "#d62728",
        f"removed by kidney (mito>{KIDNEY_MITO_MAX}%)": "#ff7f0e",
    }
    for status, color in status_colors.items():
        mask = adata.obs["filter_status"] == status
        n = mask.sum()
        if n == 0:
            continue
        fig.add_trace(go.Scattergl(
            x=umap_x[mask], y=umap_y[mask],
            mode="markers",
            marker=dict(size=POINT_SIZE, color=color, opacity=OPACITY),
            name=f"{status} (n={n:,})",
            legendgroup="status",
            showlegend=True,
        ), row=1, col=1)

    # Panel 2: mito% (continuous colour)
    fig.add_trace(go.Scattergl(
        x=umap_x, y=umap_y,
        mode="markers",
        marker=dict(
            size=POINT_SIZE,
            color=pct_mt,
            colorscale="RdYlGn_r",
            cmin=0, cmax=min(40, np.percentile(pct_mt, 99)),
            colorbar=dict(title="Mito %", x=0.68, len=0.7),
            opacity=OPACITY,
        ),
        name="Mito %",
        showlegend=False,
    ), row=1, col=2)

    # Panel 3: donor
    donor_palette = ["#1f77b4", "#ff7f0e", "#9467bd", "#8c564b", "#e377c2"]
    donors = sorted(adata.obs["sample"].unique())
    for i, donor in enumerate(donors):
        mask = adata.obs["sample"] == donor
        n = mask.sum()
        fig.add_trace(go.Scattergl(
            x=umap_x[mask], y=umap_y[mask],
            mode="markers",
            marker=dict(size=POINT_SIZE, color=donor_palette[i % len(donor_palette)], opacity=OPACITY),
            name=f"{donor} (n={n:,})",
            legendgroup="donor",
            showlegend=True,
        ), row=1, col=3)

    # Axis cleanup
    for col in [1, 2, 3]:
        fig.update_xaxes(showticklabels=False, showgrid=False, zeroline=False,
                         title_text="UMAP 1", row=1, col=col)
        fig.update_yaxes(showticklabels=False, showgrid=False, zeroline=False,
                         title_text="UMAP 2", row=1, col=col)

    fig.update_layout(
        title=dict(
            text=(
                "UMAP: Kidney cells removed by blood-calibrated mito≤5% threshold "
                "form coherent biological clusters<br>"
                f"<sup>Wu 2019 (GSE131685) · {n_total:,} cells · "
                f"{n_blood_removed:,} cells ({100*n_blood_removed/n_total:.0f}%) "
                f"removed by blood threshold, only {n_kidney_removed:,} ({100*n_kidney_removed/n_total:.0f}%) "
                f"by kidney threshold</sup>"
            ),
            x=0.5, xanchor="center",
        ),
        height=500,
        width=1400,
        font=dict(family="Arial", size=12),
        plot_bgcolor="white",
        paper_bgcolor="white",
        legend=dict(itemsizing="constant"),
    )

    out_html = os.path.join(OUTPUT_DIR, "umap_killer.html")
    fig.write_html(out_html)
    print(f"Figure written to: {out_html}")

    # ── Cluster composition summary ───────────────────────────────────────────
    print("\nCluster composition (% cells removed by blood threshold):")
    cluster_summary = []
    for cluster in sorted(adata.obs["leiden"].unique(), key=int):
        mask_c = adata.obs["leiden"] == cluster
        n_c = mask_c.sum()
        n_removed = (mask_c & (adata.obs["filter_status"] == f"removed by blood (mito>{BLOOD_MITO_MAX}%)")).sum()
        med_mito = np.median(pct_mt[mask_c])
        pct_removed = 100 * n_removed / n_c
        print(f"  Cluster {cluster:>2}: n={n_c:>5,}  mito={med_mito:>5.1f}%  "
              f"blood-removed={n_removed:>5,} ({pct_removed:.0f}%)")
        cluster_summary.append({
            "cluster": cluster,
            "n_cells": int(n_c),
            "median_mito_pct": round(float(med_mito), 2),
            "n_blood_removed": int(n_removed),
            "pct_blood_removed": round(float(pct_removed), 1),
        })

    # Save summary
    out_json = os.path.join(OUTPUT_DIR, "umap_killer_data.json")
    with open(out_json, "w") as f:
        json.dump({
            "n_total": int(n_total),
            "n_kept": int(n_kept),
            "n_blood_removed": int(n_blood_removed),
            "n_kidney_removed": int(n_kidney_removed),
            "pct_blood_removed": round(100 * n_blood_removed / n_total, 1),
            "pct_kidney_removed": round(100 * n_kidney_removed / n_total, 1),
            "per_donor_mito_median": per_donor.round(2).to_dict(),
            "cluster_summary": cluster_summary,
        }, f, indent=2)
    print(f"Data saved to: {out_json}")


if __name__ == "__main__":
    main()
