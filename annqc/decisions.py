"""Analyzes QC metric distributions and generates plain-English threshold explanations."""

import math

import numpy as np


def _isnan(v):
    try:
        return math.isnan(float(v))
    except (TypeError, ValueError):
        return False


_METRIC_LABELS = {
    "pct_counts_mt": "mitochondrial content",
    "n_genes_by_counts": "gene count per cell",
    "total_counts": "total UMI count per cell",
    "pct_counts_ribo": "ribosomal content",
}
_METRIC_UNITS = {
    "pct_counts_mt": "%",
    "pct_counts_ribo": "%",
    "n_genes_by_counts": " genes",
    "total_counts": " UMIs",
}


def analyze_distributions(adata) -> None:
    """Populate record['threshold_explanations'].

    Called from pipeline.py after thresholds are finalized. Reads raw_obs_metrics
    from the record (pre-filtering values). Results stored in-place in adata.uns['annqc'].
    """
    record = adata.uns["annqc"]
    raw = record.get("raw_obs_metrics", {})
    thr = record.get("thresholds", {})
    method = record.get("threshold_method", "manual")

    threshold_explanations = {}

    METRICS = [
        ("pct_counts_mt",      None,                    thr.get("mito_max_pct")),
        ("n_genes_by_counts",  thr.get("min_genes"),    thr.get("max_genes")),
        ("total_counts",       thr.get("min_counts"),   thr.get("max_counts")),
        ("pct_counts_ribo",    None,                    None),
    ]

    for metric, thr_low, thr_high in METRICS:
        values = raw.get(metric, [])
        if not values:
            continue
        threshold_explanations[metric] = explain_threshold(
            metric, values, thr_low, thr_high, method=method
        )

    # Doublet explanation
    dbl_thr = thr.get("doublet_threshold")
    cc = record.get("cell_counts", {})
    n_input = cc.get("input", 0)
    n_after_doublet = cc.get("after_doublet_filter", n_input)
    if dbl_thr is not None and not (isinstance(dbl_thr, float) and _isnan(dbl_thr)):
        n_doublets = n_input - n_after_doublet
        doublet_pct = 100.0 * n_doublets / n_input if n_input > 0 else 0.0
        threshold_explanations["doublet"] = (
            f"Doublet threshold set to {dbl_thr:.3f} by Scrublet (automatic). "
            f"{n_doublets} doublets identified ({doublet_pct:.1f}% of input cells)."
        )

    record["threshold_explanations"] = threshold_explanations


def explain_threshold(metric: str, values: list, threshold_low, threshold_high,
                      method: str = "manual") -> str:
    """Return a plain English explanation for the threshold applied to this metric."""
    arr = np.array(values, dtype=float)
    n = len(arr)
    if n == 0:
        return "No data available."

    median = float(np.median(arr))
    mad = float(np.median(np.abs(arr - median)))
    label = _METRIC_LABELS.get(metric, metric)
    unit = _METRIC_UNITS.get(metric, "")

    method_phrase = "MAD-based adaptive method" if method == "auto_mad" else "manual threshold"

    parts = []

    if threshold_high is not None:
        pct_below = 100.0 * np.sum(arr <= threshold_high) / n
        if method == "auto_mad":
            n_mads = (threshold_high - median) / mad if mad > 0 else 0
            parts.append(
                f"Upper threshold {threshold_high:.3g}{unit} set using {method_phrase} "
                f"(median {median:.3g}{unit} + {n_mads:.1f}×MAD {mad:.3g}{unit}). "
                f"{pct_below:.1f}% of cells fall below this threshold."
            )
        else:
            parts.append(
                f"Upper threshold {threshold_high:.3g}{unit} set manually. "
                f"{pct_below:.1f}% of cells fall below this threshold "
                f"(median {median:.3g}{unit})."
            )

    if threshold_low is not None and threshold_low > 0:
        pct_above = 100.0 * np.sum(arr >= threshold_low) / n
        if method == "auto_mad":
            n_mads_lo = (median - threshold_low) / mad if mad > 0 else 0
            parts.append(
                f"Lower threshold {threshold_low:.3g}{unit} set using {method_phrase} "
                f"(median {median:.3g}{unit} − {n_mads_lo:.1f}×MAD {mad:.3g}{unit}). "
                f"{pct_above:.1f}% of cells exceed this threshold."
            )
        else:
            parts.append(
                f"Lower threshold {threshold_low:.3g}{unit} set manually. "
                f"{pct_above:.1f}% of cells exceed this threshold "
                f"(median {median:.3g}{unit})."
            )

    if not parts:
        return f"No threshold applied for {label} (median {median:.3g}{unit})."

    return " ".join(parts)


