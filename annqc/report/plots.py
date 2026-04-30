"""
annqc.report.plots
~~~~~~~~~~~~~~~~~~
Plotly-based visualisation helpers for the annqc QC report.
Every public function returns a JSON string (fig.to_json()) or None.
"""

import math
import numpy as np
import plotly.graph_objects as go

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

_METRIC_LABELS = {
    "n_genes_by_counts": "Genes per Cell",
    "total_counts": "Total UMI Counts",
    "pct_counts_mt": "Mitochondrial %",
    "pct_counts_ribo": "Ribosomal %",
}

# Colour palette used when multiple samples are needed
_PALETTE = [
    "#4A90D9", "#E74C3C", "#2ECC71", "#F39C12", "#9B59B6",
    "#1ABC9C", "#E67E22", "#3498DB", "#E91E63", "#00BCD4",
]

_LAYOUT_DEFAULTS = dict(
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
    font=dict(family="Arial, sans-serif", size=13, color="#333333"),
    margin=dict(l=60, r=40, t=60, b=60),
)


def _metric_label(metric: str) -> str:
    return _METRIC_LABELS.get(metric, metric)


def _clean_axes(fig, xgrid: bool = False, ygrid: bool = False):
    """Apply minimal axis styling to all axes in the figure."""
    axis_style = dict(
        showgrid=False,
        zeroline=False,
        showline=True,
        linecolor="#CCCCCC",
        ticks="outside",
        tickcolor="#CCCCCC",
    )
    for ax in ("xaxis", "yaxis"):
        current = getattr(fig.layout, ax).to_plotly_json()
        fig.update_layout(
            **{ax: {**axis_style, **{k: v for k, v in current.items() if v is not None}}}
        )
    # Re-apply grid settings explicitly
    fig.update_xaxes(showgrid=xgrid)
    fig.update_yaxes(showgrid=ygrid)
    return fig


# ---------------------------------------------------------------------------
# Function 1: plot_qc_violin
# ---------------------------------------------------------------------------

def plot_qc_violin(adata, metric, threshold_low=None, threshold_high=None, sample_key=None):
    """Violin plot of adata.obs[metric], optionally split by sample."""
    if metric not in adata.obs.columns:
        return None

    label = _metric_label(metric)
    traces = []

    if sample_key is not None and sample_key in adata.obs.columns:
        samples = adata.obs[sample_key].unique().tolist()
        for i, sample in enumerate(samples):
            values = adata.obs.loc[adata.obs[sample_key] == sample, metric].dropna().tolist()
            color = _PALETTE[i % len(_PALETTE)]
            traces.append(
                go.Violin(
                    y=values,
                    name=str(sample),
                    box_visible=True,
                    meanline_visible=True,
                    fillcolor=color,
                    line_color=color,
                    opacity=0.7,
                )
            )
    else:
        values = adata.obs[metric].dropna().tolist()
        traces.append(
            go.Violin(
                y=values,
                name=label,
                box_visible=True,
                meanline_visible=True,
                fillcolor="#4A90D9",
                line_color="#4A90D9",
                opacity=0.7,
            )
        )

    fig = go.Figure(data=traces)

    # Threshold lines + annotations
    shapes = []
    annotations = []
    for tval, side in ((threshold_low, "low"), (threshold_high, "high")):
        if tval is None:
            continue
        shapes.append(
            dict(
                type="line",
                xref="paper",
                x0=0,
                x1=1,
                yref="y",
                y0=tval,
                y1=tval,
                line=dict(color="#E74C3C", width=1.5, dash="dash"),
            )
        )
        annotations.append(
            dict(
                xref="paper",
                x=1.01,
                yref="y",
                y=tval,
                text=f"{tval:g}",
                showarrow=False,
                font=dict(color="#E74C3C", size=11),
                xanchor="left",
            )
        )

    fig.update_layout(
        **_LAYOUT_DEFAULTS,
        title=dict(text=label, x=0.5, xanchor="center"),
        yaxis_title=label,
        shapes=shapes,
        annotations=annotations,
        showlegend=(sample_key is not None and sample_key in adata.obs.columns),
    )
    _clean_axes(fig)
    return fig.to_json()


# ---------------------------------------------------------------------------
# Function 2: plot_before_after
# ---------------------------------------------------------------------------

def plot_before_after(before_values: list, after_values: list, metric: str, n_removed: int = 0):
    """Side-by-side violin plots comparing a QC metric before and after filtering.

    Parameters
    ----------
    before_values : list of float
        Raw metric values for ALL cells before filtering.
    after_values : list of float
        Metric values for cells that PASSED filtering.
    metric : str
        Column name key (used for label lookup).
    n_removed : int
        Number of cells removed (shown in title).
    """
    if not before_values or not after_values:
        return None

    label = _metric_label(metric)

    if n_removed > 0:
        title_text = f"Before vs After — {label} ({n_removed:,} cells removed)"
    else:
        title_text = f"Before vs After — {label}"

    traces = [
        go.Violin(
            y=before_values,
            name="Before Filtering",
            box_visible=True,
            meanline_visible=True,
            fillcolor="#AAAAAA",
            line_color="#AAAAAA",
            opacity=0.7,
        ),
        go.Violin(
            y=after_values,
            name="After Filtering",
            box_visible=True,
            meanline_visible=True,
            fillcolor="#4A90D9",
            line_color="#4A90D9",
            opacity=0.7,
        ),
    ]

    fig = go.Figure(data=traces)
    fig.update_layout(
        **_LAYOUT_DEFAULTS,
        title=dict(
            text=title_text,
            x=0.5,
            xanchor="center",
        ),
        yaxis_title=label,
        showlegend=True,
    )
    _clean_axes(fig)
    return fig.to_json()


# ---------------------------------------------------------------------------
# Function 3: plot_doublet_scores
# ---------------------------------------------------------------------------

def plot_doublet_scores(adata, threshold=None):
    """Histogram of doublet scores, coloured by whether they exceed threshold."""
    col = "annqc_doublet_score"
    if col not in adata.obs.columns:
        return None

    scores = adata.obs[col].dropna().values
    print(f"[plot_doublet_scores] {len(scores)} non-NaN doublet scores available")
    if len(scores) == 0:
        return None

    # Convert to plain Python list — required for reliable Plotly 6 JSON serialization
    scores_list = scores.tolist()

    # Determine a sensible bin count (Sturges-ish, capped)
    n_bins = min(max(int(math.ceil(math.log2(len(scores_list)))) + 1, 20), 100)
    bin_min = float(np.min(scores))
    bin_max = float(np.max(scores))

    # Guard against degenerate range
    if bin_min == bin_max:
        bin_max = bin_min + 1e-6

    bin_width = (bin_max - bin_min) / n_bins
    # Shared xbins covering the full score range — both traces must use identical
    # bin boundaries so the overlay aligns correctly.
    xbins = dict(start=bin_min, end=bin_max + bin_width, size=bin_width)

    thr_valid = (
        threshold is not None
        and not (isinstance(threshold, float) and math.isnan(threshold))
    )

    traces = []
    if thr_valid:
        below = [s for s in scores_list if s <= threshold]
        above = [s for s in scores_list if s > threshold]
        if below:
            traces.append(
                go.Histogram(
                    x=below,
                    xbins=xbins,
                    name="Below Threshold",
                    marker_color="#4A90D9",
                    opacity=0.85,
                )
            )
        if above:
            traces.append(
                go.Histogram(
                    x=above,
                    xbins=xbins,
                    name="Above Threshold",
                    marker_color="#E74C3C",
                    opacity=0.85,
                )
            )
    else:
        traces.append(
            go.Histogram(
                x=scores_list,
                xbins=xbins,
                name="Doublet Score",
                marker_color="#4A90D9",
                opacity=0.85,
            )
        )

    fig = go.Figure(data=traces)

    shapes = []
    annotations = []
    if thr_valid:
        shapes.append(
            dict(
                type="line",
                xref="x",
                x0=threshold,
                x1=threshold,
                yref="paper",
                y0=0,
                y1=1,
                line=dict(color="#E74C3C", width=1.5, dash="dash"),
            )
        )
        annotations.append(
            dict(
                xref="x",
                x=threshold,
                yref="paper",
                y=1.02,
                text=f"Threshold: {threshold:.3f}",
                showarrow=False,
                font=dict(color="#E74C3C", size=11),
                xanchor="center",
            )
        )

    fig.update_layout(
        **_LAYOUT_DEFAULTS,
        title=dict(text="Doublet Score Distribution", x=0.5, xanchor="center"),
        xaxis_title="Doublet Score",
        yaxis_title="Number of Cells",
        barmode="overlay",
        shapes=shapes,
        annotations=annotations,
        showlegend=thr_valid,
    )
    _clean_axes(fig, ygrid=True)
    fig.update_yaxes(showgrid=True, gridcolor="#EEEEEE")
    return fig.to_json()


# ---------------------------------------------------------------------------
# Function 4: plot_filtering_waterfall
# ---------------------------------------------------------------------------

def plot_filtering_waterfall(cell_counts_dict):
    """Horizontal bar chart of cell counts at each filtering step."""
    step_keys = [
        ("input", "Input"),
        ("after_mito_filter", "After Mito Filter"),
        ("after_gene_filter", "After Gene Filter"),
        ("after_count_filter", "After Count Filter"),
        ("after_doublet_filter", "After Doublet Filter"),
        ("output", "Output"),
    ]

    labels = []
    counts = []
    for key, display in step_keys:
        val = cell_counts_dict.get(key)
        if val is not None:
            labels.append(display)
            counts.append(int(val))

    if not counts:
        return None

    input_count = counts[0] if counts[0] > 0 else 1  # avoid ZeroDivisionError
    text_labels = [
        f"{c:,} ({c / input_count * 100:.1f}%)"
        for c in counts
    ]

    fig = go.Figure(
        go.Bar(
            x=counts,
            y=labels,
            orientation="h",
            marker_color="#4A90D9",
            text=text_labels,
            textposition="outside",
            cliponaxis=False,
        )
    )

    # Extend x range slightly so outside text labels don't get clipped
    max_count = max(counts)
    layout = {**_LAYOUT_DEFAULTS, "margin": dict(l=160, r=120, t=60, b=60)}
    fig.update_layout(
        **layout,
        title=dict(
            text="Cell Counts at Each Filtering Step",
            x=0.5,
            xanchor="center",
        ),
        xaxis_title="Cells Remaining",
        xaxis_range=[0, max_count * 1.25],
        yaxis=dict(autorange="reversed"),
        showlegend=False,
    )
    _clean_axes(fig, xgrid=True)
    fig.update_xaxes(showgrid=True, gridcolor="#EEEEEE")
    return fig.to_json()


# ---------------------------------------------------------------------------
# Function 5: plot_per_sample_summary
# ---------------------------------------------------------------------------

def plot_per_sample_summary(per_sample_dict):
    """Grouped bar chart of per-sample cell retention (% of input)."""
    if not per_sample_dict:
        return None

    samples = list(per_sample_dict.keys())
    pct_retained = []
    for sample in samples:
        info = per_sample_dict[sample]
        inp = info.get("input", 0)
        out = info.get("output", 0)
        pct = (out / inp * 100) if inp and inp > 0 else 0.0
        pct_retained.append(round(pct, 2))

    fig = go.Figure(
        go.Bar(
            x=samples,
            y=pct_retained,
            marker_color="#4A90D9",
            text=[f"{p:.1f}%" for p in pct_retained],
            textposition="outside",
            cliponaxis=False,
        )
    )

    max_pct = max(pct_retained) if pct_retained else 100
    fig.update_layout(
        **_LAYOUT_DEFAULTS,
        title=dict(text="Per-Sample Cell Retention", x=0.5, xanchor="center"),
        xaxis_title=None,
        yaxis_title="% Cells Retained",
        yaxis_range=[0, min(max_pct * 1.2, 110)],
        showlegend=False,
    )
    _clean_axes(fig, ygrid=True)
    fig.update_yaxes(showgrid=True, gridcolor="#EEEEEE")
    return fig.to_json()
