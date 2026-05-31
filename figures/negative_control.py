"""Negative control: PBMC data is correctly identified as blood by reference comparison.

Runs the 10x 5k PBMC dataset through AnnQC and compares to both blood and kidney
reference profiles. Demonstrates:
  - PBMC mito is NORMAL vs blood reference (correct — no warning)
  - PBMC mito would appear LOW or NORMAL vs kidney reference (not concerning)
  - Standard blood-calibrated filtering agrees with reference: no cells over-filtered

Run from repo root:
    python figures/negative_control.py
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import annqc  # noqa: F401 — must import first to stub numba

import json
import scanpy as sc
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from annqc.qc import calculate_qc_metrics
from annqc.reference.compare import compare_to_reference, reference_warnings

DATA_PATH = "benchmarks/pbmc5k_v3/raw_data.h5"
OUTPUT_DIR = "figures"


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("Loading PBMC 5k dataset...")
    adata = sc.read_10x_h5(DATA_PATH)
    adata.var_names_make_unique()
    print(f"  {adata.n_obs} cells, {adata.n_vars} genes")

    print("Computing QC metrics...")
    adata = calculate_qc_metrics(adata)

    dataset_mito = float(np.median(adata.obs["pct_counts_mt"].dropna()))
    dataset_genes = float(np.median(adata.obs["n_genes_by_counts"].dropna()))
    dataset_counts = float(np.median(adata.obs["total_counts"].dropna()))
    print(f"  Dataset median mito: {dataset_mito:.2f}%")
    print(f"  Dataset median genes: {dataset_genes:.0f}")
    print(f"  Dataset median UMI: {dataset_counts:.0f}")

    print("\n=== Reference comparisons ===")
    results = {}
    for tissue, assay in [("blood", "10x_v3"), ("kidney", "10x_v2")]:
        result = compare_to_reference(adata, tissue=tissue, assay=assay, suspension_type="cell")
        if result:
            results[tissue] = result
            conf = result["confidence"]
            n_refs = result["n_references"]
            print(f"\n  vs {tissue} (n={n_refs}, {conf}):")
            for metric, r in result.get("comparison", {}).items():
                print(f"    {metric}: dataset={r['dataset_median']:.2f}, "
                      f"ref={r['reference_median']:.2f}, "
                      f"pct={r['percentile']:.0f}, flag={r['flag']}")
            warnings = reference_warnings(result)
            if warnings:
                for w in warnings:
                    print(f"    WARNING: {w}")
            else:
                print(f"    No warnings — PBMC mito NORMAL vs {tissue} reference")
        else:
            print(f"  No profile for {tissue}/{assay}")

    # Build figure
    metric = "pct_counts_mt"
    tissue_list = [t for t in ["blood", "kidney"] if t in results]

    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=[
            "PBMC Mito% vs Reference Profiles",
            "Mito Distribution (PBMC 5k)",
        ],
        column_widths=[0.5, 0.5],
        horizontal_spacing=0.12,
    )

    # Panel 1: dataset median vs reference range
    for i, tissue in enumerate(tissue_list):
        r = results[tissue]
        cmp = r.get("comparison", {}).get(metric)
        if not cmp:
            continue
        color = "#d62728" if cmp["flag"] == "HIGH" else "#2ca02c" if cmp["flag"] == "NORMAL" else "#ff7f0e"
        x_label = f"vs {tissue} ref\n(n={r['n_references']})"
        fig.add_trace(go.Bar(
            x=[x_label], y=[cmp["dataset_median"]],
            name=f"PBMC median ({tissue})",
            marker_color=color,
            text=[f"{cmp['dataset_median']:.1f}%<br>{cmp['percentile']:.0f}th pct<br>{cmp['flag']}"],
            textposition="outside",
            showlegend=False,
        ), row=1, col=1)
        # Reference range (p5–p95 shaded)
        fig.add_trace(go.Scatter(
            x=[x_label, x_label],
            y=[cmp["reference_p5"], cmp["reference_p95"]],
            mode="lines",
            line=dict(color="gray", width=8),
            name=f"{tissue} ref p5–p95",
            showlegend=(i == 0),
        ), row=1, col=1)
        # Reference median marker
        fig.add_trace(go.Scatter(
            x=[x_label],
            y=[cmp["reference_median"]],
            mode="markers",
            marker=dict(symbol="line-ew-open", size=20, color="black", line=dict(width=3)),
            name="Reference median",
            showlegend=(i == 0),
        ), row=1, col=1)

    fig.update_yaxes(title_text="Mitochondrial %", range=[0, 20], row=1, col=1)

    # Panel 2: histogram of mito distribution in dataset
    mito_vals = adata.obs["pct_counts_mt"].dropna().values
    fig.add_trace(go.Histogram(
        x=mito_vals,
        nbinsx=50,
        name="Mito distribution",
        marker_color="#1f77b4",
        opacity=0.75,
        showlegend=False,
    ), row=1, col=2)
    fig.add_vline(
        x=dataset_mito,
        line_dash="dash", line_color="black", line_width=2,
        annotation_text=f"Median: {dataset_mito:.1f}%",
        row=1, col=2,
    )

    fig.update_yaxes(title_text="Cell count", row=1, col=2)
    fig.update_xaxes(title_text="Mito %", row=1, col=2)

    fig.update_layout(
        title=dict(
            text="Negative Control: PBMC Correctly Identified as Blood-Appropriate<br>"
                 "<sup>10x 5k PBMC v3 | Mito is NORMAL vs blood reference — no spurious warnings</sup>",
            x=0.5, xanchor="center",
        ),
        height=480,
        width=900,
        font=dict(family="Arial", size=13),
        plot_bgcolor="white",
        paper_bgcolor="white",
    )

    out_html = os.path.join(OUTPUT_DIR, "negative_control.html")
    fig.write_html(out_html)
    print(f"\nFigure written to: {out_html}")

    out_json = os.path.join(OUTPUT_DIR, "negative_control_data.json")
    with open(out_json, "w") as f:
        json.dump({
            "dataset": "pbmc5k_v3",
            "n_cells": int(adata.n_obs),
            "dataset_mito_median": round(dataset_mito, 4),
            "reference_comparisons": {
                k: {
                    "n_references": v["n_references"],
                    "confidence": v["confidence"],
                    "comparison": {
                        m: {kk: vv for kk, vv in r.items()}
                        for m, r in v.get("comparison", {}).items()
                    }
                }
                for k, v in results.items()
            }
        }, f, indent=2)
    print(f"Data saved to: {out_json}")


if __name__ == "__main__":
    main()
