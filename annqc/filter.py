"""Cell and gene filtering logic for AnnQC."""

import logging

import scanpy as sc

logger = logging.getLogger(__name__)


def flag_cells(adata, config: dict):
    """Flag cells that fail QC thresholds without removing them.

    Adds to adata.obs:
        annqc_filter_reason (str) — first reason for failure, "" if kept
        annqc_pass          (bool) — True if cell passed all filters

    Priority order: mito > min_genes > max_genes > min_counts > max_counts > doublet.
    Does not modify cell count.
    """
    n_cells = adata.n_obs
    reason = [""] * n_cells

    mito_max = config["mito"].get("max_pct")
    min_genes = config["cells"].get("min_genes")
    max_genes = config["cells"].get("max_genes")
    min_counts = config["cells"].get("min_counts")
    max_counts = config["cells"].get("max_counts")

    def _apply(col, op, threshold, label):
        """Flag cells where obs[col] op threshold, recording label as first reason."""
        if threshold is None or col not in adata.obs.columns:
            return
        vals = adata.obs[col].values
        for i, v in enumerate(vals):
            if reason[i] == "" and op(v, threshold):
                reason[i] = label

    import operator as _op

    if mito_max is not None:
        _apply("pct_counts_mt", _op.gt, mito_max, "mito")
        logger.debug(
            "Mito filter (>%.1f%%): %d cells flagged",
            mito_max,
            sum(1 for r in reason if r == "mito"),
        )

    if min_genes is not None:
        _apply("n_genes_by_counts", _op.lt, min_genes, "min_genes")

    if max_genes is not None:
        _apply("n_genes_by_counts", _op.gt, max_genes, "max_genes")

    if min_counts is not None:
        _apply("total_counts", _op.lt, min_counts, "min_counts")

    if max_counts is not None:
        _apply("total_counts", _op.gt, max_counts, "max_counts")

    if "annqc_is_doublet" in adata.obs.columns:
        doublets = adata.obs["annqc_is_doublet"].values.astype(bool)
        for i, is_d in enumerate(doublets):
            if is_d and reason[i] == "":
                reason[i] = "doublet"

    adata.obs["annqc_filter_reason"] = reason
    adata.obs["annqc_pass"] = [r == "" for r in reason]

    n_fail = sum(1 for r in reason if r != "")
    logger.info("flag_cells: %d/%d cells flagged for removal", n_fail, n_cells)
    return adata


def apply_filters(adata):
    """Remove cells where annqc_pass is False.

    Raises ValueError if all cells would be removed or if annqc_pass is absent.
    """
    if "annqc_pass" not in adata.obs.columns:
        raise ValueError("annqc_pass column not found in obs — run flag_cells first")

    mask = adata.obs["annqc_pass"].values.astype(bool)
    n_keep = int(mask.sum())

    if n_keep == 0:
        raise ValueError(
            "All cells were removed — check your thresholds. "
            f"Input had {adata.n_obs} cells."
        )

    n_remove = adata.n_obs - n_keep
    logger.info("apply_filters: removing %d cells, keeping %d", n_remove, n_keep)
    return adata[mask].copy()


def filter_genes(adata, min_cells: int = 3):
    """Remove genes expressed in fewer than min_cells cells.

    Uses scanpy.pp.filter_genes internally.
    """
    n_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=min_cells)
    n_after = adata.n_vars
    logger.info(
        "filter_genes (min_cells=%d): removed %d genes, kept %d",
        min_cells,
        n_before - n_after,
        n_after,
    )
    return adata
