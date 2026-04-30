"""Threshold sensitivity analysis for AnnQC."""

import json
import logging
import os

import numpy as np

logger = logging.getLogger(__name__)


def run_sensitivity_analysis(
    adata_or_path,
    config=None,
    sample_key=None,
    seed: int = 0,
    output_path=None,
    report_path=None,
    profiles: bool = False,
    cluster_labels_path=None,
):
    """Test a range of QC thresholds and report their impact on cell retention and cluster composition.

    Parameters
    ----------
    adata_or_path : AnnData or str
        Input data: an AnnData, a path to .h5ad, a 10x directory, or a .h5 file.
    config : dict, str, or None
        Config dict, path to YAML, or None to use built-in defaults.
    sample_key : str or None
        obs column for per-sample breakdown.
    seed : int
        Random seed for reproducibility.
    output_path : str or None
        If provided, save results dict as JSON.
    report_path : str or None
        If provided, render and save an HTML sensitivity report.

    Returns
    -------
    dict
        Sensitivity results with keys: n_input, mito, min_genes, max_genes, min_counts.
    """
    import scanpy as sc
    from annqc.config import get_default_config
    from annqc.qc import calculate_qc_metrics
    from annqc.thresholds import suggest_thresholds
    import annqc

    # --- Step 1: Load data ---
    if isinstance(adata_or_path, (str, os.PathLike)):
        path = str(adata_or_path)
        if not os.path.exists(path):
            raise FileNotFoundError(f"Input not found: {path}")
        if os.path.isdir(path):
            adata = sc.read_10x_mtx(path, var_names="gene_symbols", cache=True)
        elif path.endswith(".h5"):
            adata = sc.read_10x_h5(path)
        else:
            adata = sc.read_h5ad(path)
    else:
        import anndata as ad
        adata = adata_or_path.copy()

    n_input = adata.n_obs
    logger.info("Sensitivity analysis: %d cells x %d genes", n_input, adata.n_vars)

    # --- Step 2: Compute baseline QC metrics ---
    cfg_defaults = get_default_config()
    adata = calculate_qc_metrics(
        adata,
        mito_prefix=cfg_defaults["mito"]["prefix"],
        ribo_prefix=cfg_defaults["ribo"]["prefix"],
    )

    # Raw metric arrays
    mt_vals = adata.obs["pct_counts_mt"].values if "pct_counts_mt" in adata.obs.columns else np.array([])
    gene_vals = adata.obs["n_genes_by_counts"].values if "n_genes_by_counts" in adata.obs.columns else np.array([])
    count_vals = adata.obs["total_counts"].values if "total_counts" in adata.obs.columns else np.array([])

    # --- Step 3: Get MAD-suggested thresholds ---
    suggested = suggest_thresholds(adata)["standard"]
    rec_mito = suggested.get("mito_max_pct")
    rec_min_genes = suggested.get("min_genes")
    rec_max_genes = suggested.get("max_genes")
    rec_min_counts = suggested.get("min_counts")

    # --- Step 4: Threshold ranges ---
    mito_range = [2, 5, 8, 10, 12, 15, 18, 20, 25, 30]
    min_genes_range = [100, 150, 200, 300, 400, 500, 600, 800, 1000]
    max_genes_range = [1500, 2000, 2500, 3000, 4000, 5000, 6000]
    min_counts_range = [200, 300, 500, 750, 1000, 1500, 2000]

    def _sweep_upper(vals, thresholds):
        removed = []
        pcts = []
        for t in thresholds:
            n_rem = int(np.sum(vals > t)) if len(vals) > 0 else 0
            removed.append(n_rem)
            pcts.append(100.0 * n_rem / n_input if n_input > 0 else 0.0)
        return removed, pcts

    def _sweep_lower(vals, thresholds):
        removed = []
        pcts = []
        for t in thresholds:
            n_rem = int(np.sum(vals < t)) if len(vals) > 0 else 0
            removed.append(n_rem)
            pcts.append(100.0 * n_rem / n_input if n_input > 0 else 0.0)
        return removed, pcts

    mito_removed, mito_pcts = _sweep_upper(mt_vals, mito_range)
    min_genes_removed, min_genes_pcts = _sweep_lower(gene_vals, min_genes_range)
    max_genes_removed, max_genes_pcts = _sweep_upper(gene_vals, max_genes_range)
    min_counts_removed, min_counts_pcts = _sweep_lower(count_vals, min_counts_range)

    # --- Step 5: Run full pipeline + clustering for cluster impact ---
    cluster_labels = None
    adata_clean = None

    try:
        adata_clean = annqc.run(
            adata_or_path,
            config=config,
            seed=seed,
        )
        logger.info("Full pipeline complete: %d cells retained", adata_clean.n_obs)

        # Use provided cluster labels if supplied, otherwise run Leiden
        if cluster_labels_path is not None:
            try:
                import pandas as pd
                cl_df = pd.read_csv(cluster_labels_path)
                cl_df = cl_df.set_index(cl_df.columns[0])
                label_col = cl_df.columns[0]
                cluster_labels = adata_clean.obs_names.map(
                    cl_df[label_col].to_dict()
                ).fillna("unknown").values
                logger.info("Loaded cluster labels from %s", cluster_labels_path)
            except Exception as exc:
                logger.warning("Failed to load cluster labels (%s) — falling back to Leiden", exc)
                cluster_labels_path = None

        if cluster_labels_path is None:
            try:
                n_comps = min(40, adata_clean.n_vars - 1, adata_clean.n_obs - 1)
                sc.pp.pca(adata_clean, n_comps=n_comps)
                sc.pp.neighbors(adata_clean, n_neighbors=10, n_pcs=n_comps)
                sc.tl.leiden(adata_clean, key_added="annqc_cluster", resolution=0.5)
                cluster_labels = adata_clean.obs["annqc_cluster"].values
                logger.info("Clustering complete: %d clusters", len(np.unique(cluster_labels)))
            except Exception as exc:
                logger.warning("Clustering failed (%s) — skipping cluster impact", exc)
                cluster_labels = None

    except Exception as exc:
        logger.warning("Full pipeline failed for cluster impact (%s) — skipping", exc)

    def _cluster_impact(metric_vals, thresholds, is_upper, cluster_labels, adata_clean):
        """For each threshold, compute fraction of each cluster's cells that would be removed."""
        impact = {}
        if cluster_labels is None or adata_clean is None or len(metric_vals) == 0:
            return impact

        # Map clean adata obs_names back to original indices
        # We use obs_names intersection if available, otherwise skip
        try:
            orig_names = set(adata.obs_names)
            clean_names = list(adata_clean.obs_names)
            clean_set = set(clean_names)

            for t in thresholds:
                key = f"threshold_{t}"
                if is_upper:
                    removed_mask = metric_vals > t
                else:
                    removed_mask = metric_vals < t

                # Map: for each cell in cleaned adata, was its original cell removed?
                removed_by_name = {
                    name: bool(removed_mask[i])
                    for i, name in enumerate(adata.obs_names)
                    if name in clean_set
                }

                cluster_fracs = {}
                for cluster in np.unique(cluster_labels):
                    c_str = str(cluster)
                    c_mask = cluster_labels == cluster
                    c_names = [n for n, m in zip(clean_names, c_mask) if m]
                    if not c_names:
                        continue
                    n_would_remove = sum(1 for n in c_names if removed_by_name.get(n, False))
                    frac = n_would_remove / len(c_names)
                    cluster_fracs[c_str] = round(frac, 4)

                impact[key] = cluster_fracs
        except Exception as exc:
            logger.warning("Cluster impact calculation failed: %s", exc)

        return impact

    mito_impact = _cluster_impact(mt_vals, mito_range, True, cluster_labels, adata_clean)
    min_genes_impact = _cluster_impact(gene_vals, min_genes_range, False, cluster_labels, adata_clean)
    max_genes_impact = _cluster_impact(gene_vals, max_genes_range, True, cluster_labels, adata_clean)
    min_counts_impact = _cluster_impact(count_vals, min_counts_range, False, cluster_labels, adata_clean)

    results = {
        "n_input": n_input,
        "mito": {
            "thresholds_tested": mito_range,
            "cells_removed": mito_removed,
            "pct_removed": mito_pcts,
            "cluster_impact": mito_impact,
            "mad_suggested": float(rec_mito) if rec_mito is not None else None,
            "mad_suggested_reason": "MAD-based: median + 5×MAD",
        },
        "min_genes": {
            "thresholds_tested": min_genes_range,
            "cells_removed": min_genes_removed,
            "pct_removed": min_genes_pcts,
            "cluster_impact": min_genes_impact,
            "mad_suggested": float(rec_min_genes) if rec_min_genes is not None else None,
            "mad_suggested_reason": "MAD-based: median − 5×MAD",
        },
        "max_genes": {
            "thresholds_tested": max_genes_range,
            "cells_removed": max_genes_removed,
            "pct_removed": max_genes_pcts,
            "cluster_impact": max_genes_impact,
            "mad_suggested": float(rec_max_genes) if rec_max_genes is not None else None,
            "mad_suggested_reason": "MAD-based: median + 5×MAD",
        },
        "min_counts": {
            "thresholds_tested": min_counts_range,
            "cells_removed": min_counts_removed,
            "pct_removed": min_counts_pcts,
            "cluster_impact": min_counts_impact,
            "mad_suggested": float(rec_min_counts) if rec_min_counts is not None else None,
            "mad_suggested_reason": "MAD-based: median − 5×MAD",
        },
    }

    # --- Step 6: Profile comparison (optional) ---
    if profiles:
        results["profiles"] = _compute_profiles(adata_or_path, adata, config, seed, sample_key)

    # --- Step 7: Save JSON ---
    if output_path is not None:
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as fh:
            json.dump(results, fh, indent=2)
        logger.info("Sensitivity results saved → %s", output_path)

    # --- Step 8: Build report ---
    if report_path is not None:
        try:
            _build_sensitivity_report(results, report_path)
        except Exception as exc:
            logger.warning("Sensitivity report generation failed: %s", exc)

    return results


def _compute_profiles(adata_or_path, adata_with_metrics, config, seed, sample_key):
    """Run strict / standard / permissive MAD profiles and return comparison data."""
    import annqc
    from annqc.thresholds import suggest_thresholds

    all_suggestions = suggest_thresholds(adata_with_metrics)
    profile_rows = []

    for level in ("strict", "standard", "permissive"):
        thr = all_suggestions[level]
        profile_cfg = {
            "mito": {"max_pct": thr["mito_max_pct"] if thr["mito_max_pct"] is not None else 20},
            "cells": {
                "min_genes": int(thr["min_genes"]) if thr["min_genes"] is not None else 200,
                "max_genes": int(thr["max_genes"]) if thr["max_genes"] is not None else 6000,
                "min_counts": int(thr["min_counts"]) if thr["min_counts"] is not None else 500,
                "max_counts": None,
            },
        }

        try:
            adata_p = annqc.run(
                adata_or_path,
                config=profile_cfg,
                seed=seed,
                dry_run=False,
            )
            cc = adata_p.uns["annqc"]["cell_counts"]
            n_in = cc["input"]
            n_out = cc["output"]

            # Per-step breakdown from waterfall counts
            row = {
                "profile": level,
                "mito_max_pct": round(float(thr["mito_max_pct"]), 2) if thr["mito_max_pct"] is not None else None,
                "min_genes": int(thr["min_genes"]) if thr["min_genes"] is not None else None,
                "max_genes": int(thr["max_genes"]) if thr["max_genes"] is not None else None,
                "min_counts": int(thr["min_counts"]) if thr["min_counts"] is not None else None,
                "n_input": n_in,
                "n_output": n_out,
                "pct_retained": round(100.0 * n_out / n_in, 1) if n_in > 0 else 0.0,
                "cells_removed_mito": n_in - cc.get("after_mito_filter", n_in),
                "cells_removed_genes": cc.get("after_mito_filter", n_in) - cc.get("after_gene_filter", cc.get("after_mito_filter", n_in)),
                "cells_removed_counts": cc.get("after_gene_filter", n_in) - cc.get("after_count_filter", cc.get("after_gene_filter", n_in)),
                "doublets_removed": cc.get("after_count_filter", n_in) - cc.get("after_doublet_filter", cc.get("after_count_filter", n_in)),
            }

            # Per-sample impact if sample_key provided
            if sample_key and sample_key in adata_p.obs.columns:
                row["per_sample"] = {
                    str(s): int((adata_p.obs[sample_key] == s).sum())
                    for s in adata_p.obs[sample_key].unique()
                }

        except Exception as exc:
            logger.warning("Profile '%s' run failed: %s", level, exc)
            row = {"profile": level, "error": str(exc)}

        profile_rows.append(row)
        logger.info("Profile '%s': %s cells retained", level, row.get("n_output", "?"))

    return profile_rows


def _build_sensitivity_report(results: dict, output_path: str) -> None:
    """Render the sensitivity HTML report."""
    import jinja2
    import plotly.graph_objects as go

    PALETTE = ["#4A90D9", "#E74C3C", "#27AE60", "#F39C12", "#9B59B6"]
    METRIC_LABELS = {
        "mito": "Mitochondrial % (upper bound)",
        "min_genes": "Min Genes per Cell (lower bound)",
        "max_genes": "Max Genes per Cell (upper bound)",
        "min_counts": "Min UMI Counts (lower bound)",
    }

    def _make_chart(metric_key, color):
        m = results.get(metric_key, {})
        xs = m.get("thresholds_tested", [])
        ys = m.get("pct_removed", [])
        rec = m.get("mad_suggested")
        label = METRIC_LABELS.get(metric_key, metric_key)

        traces = [
            go.Scatter(
                x=xs,
                y=ys,
                mode="lines+markers",
                name="% cells removed",
                line=dict(color=color, width=2),
                marker=dict(size=7),
                hovertemplate="Threshold: %{x}<br>Removed: %{y:.1f}%<extra></extra>",
            )
        ]

        shapes = []
        annotations = []

        # Color bands
        x_min = min(xs) if xs else 0
        x_max = max(xs) if xs else 1
        for y0, y1, fill in [(0, 10, "rgba(39,174,96,0.10)"),
                              (10, 20, "rgba(243,156,18,0.10)"),
                              (20, 100, "rgba(231,76,60,0.10)")]:
            shapes.append(dict(
                type="rect", xref="paper", yref="y",
                x0=0, x1=1, y0=y0, y1=y1,
                fillcolor=fill, line_width=0, layer="below",
            ))

        if rec is not None:
            shapes.append(dict(
                type="line", xref="x", yref="paper",
                x0=rec, x1=rec, y0=0, y1=1,
                line=dict(color="#2C3E50", width=2, dash="dash"),
            ))
            annotations.append(dict(
                x=rec, y=1.05, xref="x", yref="paper",
                text=f"MAD-suggested: {rec:.3g}",
                showarrow=False,
                font=dict(size=11, color="#2C3E50"),
                xanchor="center",
            ))

        fig = go.Figure(data=traces)
        fig.update_layout(
            title=dict(text=label, x=0.5, xanchor="center", font=dict(size=13)),
            xaxis_title="Threshold Value",
            yaxis_title="% Cells Removed",
            yaxis=dict(range=[0, min(100, max(ys or [0]) * 1.3 + 5)]),
            shapes=shapes,
            annotations=annotations,
            margin=dict(l=50, r=30, t=60, b=50),
            height=320,
            plot_bgcolor="white",
            paper_bgcolor="white",
            font=dict(family="system-ui, sans-serif", size=12, color="#2C3E50"),
            showlegend=False,
        )
        fig.update_xaxes(showgrid=True, gridcolor="#E8ECF0", zeroline=False)
        fig.update_yaxes(showgrid=True, gridcolor="#E8ECF0", zeroline=False)
        return fig.to_json()

    charts = {
        "mito_chart": _make_chart("mito", PALETTE[0]),
        "min_genes_chart": _make_chart("min_genes", PALETTE[2]),
        "max_genes_chart": _make_chart("max_genes", PALETTE[3]),
        "min_counts_chart": _make_chart("min_counts", PALETTE[1]),
    }

    # Pre-compute MAD-suggested summaries (avoids complex Jinja2 list-building)
    rec_summaries = {}
    for metric_key in ("mito", "min_genes", "max_genes", "min_counts"):
        m = results.get(metric_key, {})
        rec = m.get("mad_suggested")
        if rec is None:
            continue
        disproportionate = []
        for thr_key, cluster_fracs in m.get("cluster_impact", {}).items():
            thr_val = float(thr_key.replace("threshold_", ""))
            if thr_val < rec:
                for cluster, frac in cluster_fracs.items():
                    if frac > 0.10:
                        disproportionate.append(f"Cluster {cluster}")
        rec_summaries[metric_key] = {
            "rec": rec,
            "reason": m.get("mad_suggested_reason", ""),
            "label": METRIC_LABELS.get(metric_key, metric_key),
            "disproportionate": list(dict.fromkeys(disproportionate)),
        }

    # Build profile summary sentence if profiles exist
    profiles_data = results.get("profiles", [])
    profiles_summary = ""
    if len(profiles_data) >= 2:
        strict_row = next((r for r in profiles_data if r.get("profile") == "strict"), None)
        standard_row = next((r for r in profiles_data if r.get("profile") == "standard"), None)
        if strict_row and standard_row and "pct_retained" in strict_row and "pct_retained" in standard_row:
            diff = round(standard_row["pct_retained"] - strict_row["pct_retained"], 1)
            # Find which metric causes the biggest difference between strict and standard
            metric_diffs = {}
            for key in ("cells_removed_mito", "cells_removed_genes", "cells_removed_counts", "doublets_removed"):
                s_val = strict_row.get(key, 0) or 0
                st_val = standard_row.get(key, 0) or 0
                metric_diffs[key] = abs(s_val - st_val)
            biggest = max(metric_diffs, key=metric_diffs.get)
            label_map = {
                "cells_removed_mito": "mitochondrial",
                "cells_removed_genes": "gene count",
                "cells_removed_counts": "UMI count",
                "doublets_removed": "doublet",
            }
            profiles_summary = (
                f"Profile comparison shows that strict thresholds remove {diff}% "
                f"more cells than standard. "
                f"The largest difference is in {label_map.get(biggest, biggest)} filtering."
            )

    template_dir = os.path.join(os.path.dirname(__file__), "report", "templates")
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(template_dir),
        autoescape=False,
    )
    template = env.get_template("sensitivity.html")
    rendered = template.render(
        results=results,
        metric_labels=METRIC_LABELS,
        rec_summaries=rec_summaries,
        n_input_fmt=f"{results.get('n_input', 0):,}",
        profiles_data=profiles_data,
        profiles_summary=profiles_summary,
        **charts,
    )

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(rendered)
    logger.info("Sensitivity report written → %s", output_path)
