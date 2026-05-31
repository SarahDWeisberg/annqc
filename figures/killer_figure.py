"""Killer figure: reference-aware QC catches tissue-inappropriate filtering.

Runs the Wu 2019 kidney dataset (GSE131685) through two AnnQC configurations:
  A) Naive blood-calibrated: mito_max=5% (typical PBMC QC threshold)
  B) Kidney reference-aware: tissue_preset=kidney + reference comparison

Shows:
  1. Cell retention waterfall for each approach
  2. Reference comparison table: kidney mito is HIGH vs blood but NORMAL vs kidney
  3. Bar chart showing cells lost under each regime

Run from repo root:
    python figures/killer_figure.py
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import annqc  # noqa: F401 — must import first to stub numba

import copy
import json
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from annqc.config import DEFAULT_CONFIG
from annqc.pipeline import run as annqc_run
from annqc.reference.compare import compare_to_reference, reference_warnings
from annqc.qc import calculate_qc_metrics

DATA_PATH = "benchmarks/wu_kidney/raw_data.h5ad"
OUTPUT_DIR = "figures"


def _make_config(mito_max: float) -> dict:
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    cfg["mito"]["max_pct"] = mito_max
    cfg["cells"]["min_genes"] = 200
    cfg["cells"]["min_counts"] = 500
    cfg["doublets"]["enabled"] = False
    return cfg


def run_pipeline_with_config(cfg: dict, label: str):
    """Run AnnQC pipeline on kidney data with the given config. Returns record dict."""
    import tempfile, scanpy as sc
    print(f"\n=== {label} ===")
    adata = sc.read_h5ad(DATA_PATH)
    with tempfile.NamedTemporaryFile(suffix=".yaml", mode="w", delete=False) as f:
        import yaml
        yaml.dump(cfg, f)
        config_path = f.name
    out_h5ad = f"/tmp/kidney_fig_{label.replace(' ', '_')}.h5ad"
    try:
        result_adata = annqc_run(
            adata_or_path=DATA_PATH,
            config=config_path,
            output=out_h5ad,
            input_file=DATA_PATH,
            dry_run=False,
            auto_thresholds=False,
            no_doublet_detection=True,
            no_cluster_qc=True,
        )
        record = result_adata.uns["annqc"]
        cc = record["cell_counts"]
        print(f"  Input: {cc['input']} cells")
        print(f"  After mito filter: {cc['after_mito_filter']} cells")
        print(f"  Output: {cc['output']} cells")
        pct_kept = 100 * cc["output"] / cc["input"]
        pct_removed = 100 - pct_kept
        print(f"  Removed: {cc['input'] - cc['output']} cells ({pct_removed:.1f}%)")
        return record, result_adata
    except Exception as e:
        print(f"  ERROR: {e}")
        raise
    finally:
        os.unlink(config_path)


def run_reference_comparison(data_path: str):
    """Compute reference comparison against both blood and kidney profiles."""
    import scanpy as sc
    adata = sc.read_h5ad(data_path)
    adata = calculate_qc_metrics(adata)

    print("\n=== Reference comparison ===")
    results = {}
    for tissue, assay in [("blood", "10x_v3"), ("kidney", "10x_v2")]:
        result = compare_to_reference(adata, tissue=tissue, assay=assay, suspension_type="cell")
        if result:
            results[tissue] = result
            print(f"\n  vs {tissue} reference (n={result['n_references']}, {result['confidence']}):")
            for metric, r in result.get("comparison", {}).items():
                print(f"    {metric}: median={r['dataset_median']:.2f}, "
                      f"ref_median={r['reference_median']:.2f}, "
                      f"percentile={r['percentile']:.0f}th, flag={r['flag']}")
            for w in reference_warnings(result):
                print(f"  WARNING: {w}")
        else:
            print(f"  No profile found for {tissue}/{assay}")
    return results


def make_figure(record_a: dict, record_b: dict, ref_results: dict, label_a: str, label_b: str):
    """Generate a 2-panel figure: cell retention + reference comparison."""
    cc_a = record_a["cell_counts"]
    cc_b = record_b["cell_counts"]

    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=["Cell Retention by QC Strategy", "Mito% vs Reference Atlas"],
        column_widths=[0.45, 0.55],
        horizontal_spacing=0.12,
    )

    # Panel 1: bar chart of cells kept vs removed
    strategies = [label_a, label_b]
    kept = [cc_a["output"], cc_b["output"]]
    removed = [cc_a["input"] - cc_a["output"], cc_b["input"] - cc_b["output"]]
    total = cc_a["input"]

    fig.add_trace(go.Bar(
        name="Kept", x=strategies, y=[k/total*100 for k in kept],
        marker_color="#2ca02c", text=[f"{k:,}" for k in kept],
        textposition="inside",
    ), row=1, col=1)
    fig.add_trace(go.Bar(
        name="Removed", x=strategies, y=[r/total*100 for r in removed],
        marker_color="#d62728", text=[f"{r:,}" for r in removed],
        textposition="inside",
    ), row=1, col=1)

    fig.update_yaxes(title_text="Cells (%)", range=[0, 100], row=1, col=1)
    fig.update_layout(barmode="stack")

    # Panel 2: mito percentile comparison
    metric = "pct_counts_mt"
    tissues = list(ref_results.keys())
    if not tissues:
        pass
    else:
        ref_labels = []
        dataset_medians = []
        ref_medians = []
        percentiles = []
        flags = []
        for tissue in tissues:
            r = ref_results[tissue].get("comparison", {}).get(metric)
            if r:
                ref_labels.append(f"vs {tissue}\nreference")
                dataset_medians.append(r["dataset_median"])
                ref_medians.append(r["reference_median"])
                percentiles.append(r["percentile"] or 0)
                flags.append(r["flag"])

        colors = ["#d62728" if f == "HIGH" else "#2ca02c" if f == "NORMAL" else "#1f77b4"
                  for f in flags]

        fig.add_trace(go.Bar(
            name="Dataset mito median",
            x=ref_labels,
            y=dataset_medians,
            marker_color=colors,
            text=[f"{v:.1f}% ({p:.0f}th pct)" for v, p in zip(dataset_medians, percentiles)],
            textposition="outside",
        ), row=1, col=2)

        # Reference median as horizontal reference lines (as scatter)
        for i, (label, ref_med) in enumerate(zip(ref_labels, ref_medians)):
            fig.add_trace(go.Scatter(
                x=[label],
                y=[ref_med],
                mode="markers",
                marker=dict(symbol="line-ew-open", size=20, color="black", line=dict(width=3)),
                name=f"Ref median ({label.split(chr(10))[1].strip()})",
                showlegend=(i == 0),
            ), row=1, col=2)

    fig.update_yaxes(title_text="Mitochondrial % (dataset median)", row=1, col=2)

    fig.update_layout(
        title=dict(
            text="AnnQC: Reference-Aware QC Catches Tissue-Inappropriate Filtering<br>"
                 "<sup>Wu 2019 kidney (GSE131685, n=25,404 cells) | Strategy A = blood-calibrated mito≤5% | "
                 "Strategy B = kidney preset (mito≤30%)</sup>",
            x=0.5, xanchor="center",
        ),
        height=500,
        width=1000,
        font=dict(family="Arial", size=13),
        plot_bgcolor="white",
        paper_bgcolor="white",
    )

    return fig


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Strategy A: blood-calibrated (5% mito, typical PBMC threshold)
    cfg_a = _make_config(mito_max=5.0)
    record_a, _ = run_pipeline_with_config(cfg_a, "A: blood-calibrated (mito≤5%)")

    # Strategy B: kidney preset (30% mito)
    cfg_b = _make_config(mito_max=30.0)
    record_b, _ = run_pipeline_with_config(cfg_b, "B: kidney preset (mito≤30%)")

    # Reference comparison on raw data
    ref_results = run_reference_comparison(DATA_PATH)

    # Summary table
    cc_a = record_a["cell_counts"]
    cc_b = record_b["cell_counts"]
    total = cc_a["input"]
    print("\n=== Summary ===")
    print(f"{'Strategy':<35} {'Input':>8} {'Output':>8} {'Removed':>8} {'% kept':>8}")
    print("-" * 75)
    print(f"{'A: blood-calibrated (mito≤5%)':<35} {total:>8,} {cc_a['output']:>8,} "
          f"{total - cc_a['output']:>8,} {100*cc_a['output']/total:>7.1f}%")
    print(f"{'B: kidney preset (mito≤30%)':<35} {total:>8,} {cc_b['output']:>8,} "
          f"{total - cc_b['output']:>8,} {100*cc_b['output']/total:>7.1f}%")

    # Generate figure
    fig = make_figure(record_a, record_b, ref_results,
                      "A: blood (mito≤5%)", "B: kidney (mito≤30%)")
    out_html = os.path.join(OUTPUT_DIR, "killer_figure.html")
    fig.write_html(out_html)
    print(f"\nFigure written to: {out_html}")

    # Save raw results
    out_json = os.path.join(OUTPUT_DIR, "killer_figure_data.json")
    summary = {
        "strategy_a": {
            "label": "blood-calibrated (mito≤5%)",
            "cell_counts": cc_a,
            "thresholds": record_a["thresholds"],
        },
        "strategy_b": {
            "label": "kidney preset (mito≤30%)",
            "cell_counts": cc_b,
            "thresholds": record_b["thresholds"],
        },
        "reference_comparison": {
            k: {
                "n_references": v["n_references"],
                "confidence": v["confidence"],
                "comparison": {
                    m: {kk: vv for kk, vv in r.items() if kk != "n_references"}
                    for m, r in v.get("comparison", {}).items()
                }
            }
            for k, v in ref_results.items()
        }
    }
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Data saved to: {out_json}")


if __name__ == "__main__":
    main()
