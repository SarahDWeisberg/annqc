"""MAD-based threshold suggestion for AnnQC."""

import numpy as np


def suggest_thresholds(
    adata,
    n_mads: int = 5,
    level: str = "standard",
    mito_min_pct: float = 10.0,
    min_genes_floor: int = 200,
) -> dict:
    """Suggest QC thresholds using MAD (upper) and percentile (lower) bounds.

    Parameters
    ----------
    adata : AnnData
        AnnData object with QC metrics already computed in adata.obs.
    n_mads : int
        Unused; kept for backward compatibility.
    level : str
        Convenience selector: ``"strict"``, ``"standard"``, or ``"permissive"``.
        Controls which level is also returned under the ``"selected"`` key.
    mito_min_pct : float
        Floor applied to the auto mito upper bound. Prevents over-filtering in
        tissues with naturally elevated mitochondrial content (e.g. gut, heart).
        Default 10.0.
    min_genes_floor : int
        Floor applied to the auto min_genes lower bound. Prevents accepting
        cells with near-zero gene counts as valid. Default 200.

    Returns
    -------
    dict
        Keys: ``"strict"`` (p95 mito / p2.5 lower), ``"standard"``
        (p97.5 mito / p5 lower), ``"permissive"`` (p99 mito / p10 lower),
        ``"selected"`` (same dict as the chosen level), and ``"raw"``
        (median / MAD / vals per metric).
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

    # ---------- compute raw median / MAD / vals for each metric ----------
    raw: dict = {}
    for col in _METRICS:
        if col in adata.obs.columns:
            vals = adata.obs[col].values.astype(float)
            median = float(np.median(vals))
            mad = float(np.median(np.abs(vals - median)))
            raw[col] = {"median": median, "mad": mad, "vals": vals}

    def _bounds(col, n, lo_pct):
        """Return (lower, upper) MAD bounds.

        Upper bound: median + n * MAD (clamped to 100 for pct metrics;
        clamped to mito_min_pct floor for pct_counts_mt).
        Lower bound: lo_pct percentile, floored at min_genes_floor for genes.

        Returns (None, None) when the metric is absent or all-zero.
        Lower bound returns None when the computed percentile is 0.0.
        """
        if col not in raw:
            return None, None
        median = raw[col]["median"]
        mad = raw[col]["mad"]
        vals = raw[col]["vals"]
        if median == 0.0 and mad == 0.0:
            return None, None
        lo = float(np.percentile(vals, lo_pct))
        if lo == 0.0:
            lo = None
        elif col == "n_genes_by_counts":
            lo = max(lo, float(min_genes_floor))
        hi = median + n * mad
        if col in _PCT_METRICS:
            hi = min(100.0, hi)
        if col == "pct_counts_mt":
            hi = max(hi, mito_min_pct)
        return lo, hi

    def _level(n, lo_pct):
        lo_mt, hi_mt = _bounds("pct_counts_mt", n, lo_pct)
        lo_genes, hi_genes = _bounds("n_genes_by_counts", n, lo_pct)
        lo_counts, hi_counts = _bounds("total_counts", n, lo_pct)
        return {
            "mito_max_pct": hi_mt,
            "min_genes": lo_genes,
            "max_genes": hi_genes,
            "min_counts": lo_counts,
            "max_counts": hi_counts,
        }

    # Strip vals from raw before returning (keep only median/mad for public API)
    raw_public = {
        col: {"median": v["median"], "mad": v["mad"]}
        for col, v in raw.items()
    }

    strict     = _level(3, 2.5)
    standard   = _level(5, 5.0)
    permissive = _level(7, 10.0)

    levels = {
        "strict": strict,
        "standard": standard,
        "permissive": permissive,
    }

    return {
        "strict": strict,
        "standard": standard,
        "permissive": permissive,
        "selected": levels.get(level, standard),
        "raw": raw_public,
    }
