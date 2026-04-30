"""MAD-based threshold suggestion for AnnQC."""

import numpy as np


def suggest_thresholds(adata, n_mads: int = 5) -> dict:
    """Suggest QC thresholds based on median ± N × MAD for each QC metric.

    Parameters
    ----------
    adata : AnnData
        AnnData object with QC metrics already computed in adata.obs.
    n_mads : int
        Default MAD multiplier (used for the ``"standard"`` level).
        The three preset levels always use 3 / 5 / 7.

    Returns
    -------
    dict
        Keys: ``"strict"`` (3-MAD), ``"standard"`` (5-MAD), ``"permissive"``
        (7-MAD), and ``"raw"`` (median/MAD values per metric).
        Each threshold level contains:
        ``mito_max_pct``, ``min_genes``, ``max_genes``, ``min_counts``,
        ``max_counts`` — all as float or None when the metric is absent.
    """
    _METRICS = (
        "pct_counts_mt",
        "n_genes_by_counts",
        "total_counts",
        "pct_counts_ribo",
    )
    _PCT_METRICS = {"pct_counts_mt", "pct_counts_ribo"}

    # ---------- compute raw median / MAD for each metric ----------
    raw: dict = {}
    for col in _METRICS:
        if col in adata.obs.columns:
            vals = adata.obs[col].values.astype(float)
            median = float(np.median(vals))
            mad = float(np.median(np.abs(vals - median)))
            raw[col] = {"median": median, "mad": mad}

    def _bounds(col, n):
        """Return (lower, upper) bounds clamped appropriately."""
        if col not in raw:
            return None, None
        median = raw[col]["median"]
        mad = raw[col]["mad"]
        lo = max(0.0, median - n * mad)
        hi = median + n * mad
        if col in _PCT_METRICS:
            hi = min(100.0, hi)
        return lo, hi

    def _level(n):
        lo_mt, hi_mt = _bounds("pct_counts_mt", n)
        lo_genes, hi_genes = _bounds("n_genes_by_counts", n)
        lo_counts, hi_counts = _bounds("total_counts", n)
        return {
            "mito_max_pct": hi_mt,
            "min_genes": lo_genes,
            "max_genes": hi_genes,
            "min_counts": lo_counts,
            "max_counts": hi_counts,
        }

    return {
        "strict": _level(3),
        "standard": _level(5),
        "permissive": _level(7),
        "raw": raw,
    }
