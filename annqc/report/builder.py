import logging
import math
import os

import jinja2
import yaml

from annqc import utils
from annqc.report import plots

logger = logging.getLogger(__name__)


def build_report(adata, output_path: str) -> None:
    """Build an HTML QC report and write it to output_path."""
    record = adata.uns["annqc"]
    cfg = record["config"]

    def _round_thr(v):
        """Round a threshold float to 3 dp; leave None/NaN unchanged."""
        if v is None:
            return None
        if isinstance(v, float) and math.isnan(v):
            return v
        try:
            return round(float(v), 3)
        except (TypeError, ValueError):
            return v

    # Round ALL threshold values to 3 dp for display and for plot annotations
    thr = {
        k: (_round_thr(v) if k != "method" else v)
        for k, v in record["thresholds"].items()
    }
    doublet_threshold = thr.get("doublet_threshold")

    # Build plots — use rounded thresholds so plot annotations match displayed values
    mito_violin = plots.plot_qc_violin(
        adata,
        "pct_counts_mt",
        threshold_high=thr["mito_max_pct"],
    )
    genes_violin = plots.plot_qc_violin(
        adata,
        "n_genes_by_counts",
        threshold_low=thr["min_genes"],
        threshold_high=thr["max_genes"],
    )
    counts_violin = plots.plot_qc_violin(
        adata,
        "total_counts",
        threshold_low=thr["min_counts"],
    )
    ribo_violin = plots.plot_qc_violin(adata, "pct_counts_ribo")
    doublet_scores = plots.plot_doublet_scores(
        adata,
        threshold=doublet_threshold,
    )
    waterfall = plots.plot_filtering_waterfall(record["cell_counts"])
    per_sample_plot = plots.plot_per_sample_summary(record.get("per_sample", {}))

    # Before/after plots (uses raw_obs_metrics snapshot from pipeline)
    raw = record.get("raw_obs_metrics", {})
    before_after_plots = {}
    for metric in ("pct_counts_mt", "n_genes_by_counts", "total_counts", "pct_counts_ribo"):
        before_vals = raw.get(metric, [])
        after_vals = adata.obs[metric].tolist() if metric in adata.obs.columns else []
        n_removed = len(before_vals) - len(after_vals)
        if before_vals and after_vals:
            before_after_plots[metric] = plots.plot_before_after(
                before_vals, after_vals, metric, n_removed=n_removed
            )
        else:
            before_after_plots[metric] = None

    # Generate methods text
    try:
        from annqc.methods_text import generate_methods_text, generate_methods_short
        methods_text_full = generate_methods_text(adata)
        methods_text_short = generate_methods_short(adata)
    except Exception:
        methods_text_full = ""
        methods_text_short = ""

    # Build template data
    data = {
        "title": cfg.get("report", {}).get("title", "Single-Cell QC Report"),
        "author": cfg.get("report", {}).get("author", ""),
        "run_date": record["date"],
        "version": record["version"],
        "seed": record["seed"],
        "input_file": record["input_file"],
        "status": record["status"],
        "cell_counts": record["cell_counts"],
        "thresholds": thr,
        "warnings": record["warnings"],
        "per_sample": record.get("per_sample", {}),
        "software_versions": utils.get_software_versions(),
        "config_yaml": yaml.dump(record["config"], default_flow_style=False, sort_keys=False),
        "dry_run": record.get("dry_run", False),
        "doublet_status": record.get("doublet_status"),
        "threshold_method": record.get("threshold_method", "manual"),
        "before_after_plots": before_after_plots,
        "qc_arrays": record.get("raw_obs_metrics", {}),
        "threshold_explanations": record.get("threshold_explanations", {}),
        "methods_text_full": methods_text_full,
        "methods_text_short": methods_text_short,
        "plots": {
            "mito_violin": mito_violin,
            "genes_violin": genes_violin,
            "counts_violin": counts_violin,
            "ribo_violin": ribo_violin,
            "doublet_scores": doublet_scores,
            "waterfall": waterfall,
            "per_sample": per_sample_plot,
        },
    }

    # Load and render Jinja2 template
    template_dir = os.path.join(os.path.dirname(__file__), "templates")
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(template_dir),
        autoescape=False,
    )
    template = env.get_template("report.html")
    rendered = template.render(**data)

    # Write output file, creating parent directories if needed
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(rendered)

    logger.info("Report written to %s", output_path)


def _plot_comparison_violin(vals_a: list, vals_b: list, metric: str, label_a: str, label_b: str) -> str:
    """Two-violin comparison plot for the comparison report."""
    from annqc.report.plots import _metric_label, _PALETTE, _LAYOUT_DEFAULTS, _clean_axes
    import plotly.graph_objects as go
    label = _metric_label(metric)
    traces = []
    if vals_a:
        traces.append(go.Violin(
            y=vals_a, name=label_a, box_visible=True, meanline_visible=True,
            fillcolor=_PALETTE[0], line_color=_PALETTE[0], opacity=0.7,
        ))
    if vals_b:
        traces.append(go.Violin(
            y=vals_b, name=label_b, box_visible=True, meanline_visible=True,
            fillcolor=_PALETTE[1], line_color=_PALETTE[1], opacity=0.7,
        ))
    if not traces:
        return None
    fig = go.Figure(data=traces)
    fig.update_layout(
        **_LAYOUT_DEFAULTS,
        title=dict(text=f"{label}: {label_a} vs {label_b}", x=0.5, xanchor="center"),
        yaxis_title=label,
        showlegend=True,
    )
    _clean_axes(fig)
    return fig.to_json()


def build_comparison_report(adata_a, adata_b, label_a: str, label_b: str, output_path: str) -> None:
    """Build a side-by-side comparison HTML report for two AnnQC runs."""
    rec_a = adata_a.uns["annqc"]
    rec_b = adata_b.uns["annqc"]

    # Build comparison violin plots for QC metrics
    comparison_plots = {}
    for metric in ("pct_counts_mt", "n_genes_by_counts", "total_counts", "pct_counts_ribo"):
        vals_a = adata_a.obs[metric].dropna().tolist() if metric in adata_a.obs.columns else []
        vals_b = adata_b.obs[metric].dropna().tolist() if metric in adata_b.obs.columns else []
        if vals_a or vals_b:
            comparison_plots[metric] = _plot_comparison_violin(vals_a, vals_b, metric, label_a, label_b)
        else:
            comparison_plots[metric] = None

    data = {
        "title": f"AnnQC Comparison: {label_a} vs {label_b}",
        "label_a": label_a,
        "label_b": label_b,
        "record_a": rec_a,
        "record_b": rec_b,
        "comparison_plots": comparison_plots,
        "software_versions": utils.get_software_versions(),
    }

    template_dir = os.path.join(os.path.dirname(__file__), "templates")
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(template_dir),
        autoescape=False,
    )
    template = env.get_template("comparison.html")
    rendered = template.render(**data)

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(rendered)
    logger.info("Comparison report written to %s", output_path)
