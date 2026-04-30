"""Scrublet-based doublet detection for AnnQC."""

import logging

import numpy as np

logger = logging.getLogger(__name__)


def detect_doublets(adata, threshold="auto", simulate_doublet_ratio: float = 2.0, seed: int = 0):
    """Run Scrublet doublet detection on adata.

    Adds to adata.obs:
        annqc_doublet_score (float) — scrublet doublet probability score
        annqc_is_doublet    (bool)  — True if cell is predicted doublet

    If scrublet fails for any reason, logs a warning and falls back to
    NaN scores and False predictions so the pipeline can continue.

    The resolved doublet threshold (float) is stored in
    adata.uns['annqc']['thresholds']['doublet_threshold'] if the key exists.
    """
    n_cells = adata.n_obs
    logger.info("Running Scrublet on %d cells (seed=%d, threshold=%s)", n_cells, seed, threshold)

    _elements = adata.n_obs * adata.n_vars
    if _elements > 500_000_000:
        _estimated_gb = (_elements * 4) / 1e9
        logger.warning(
            "Large matrix detected (%d cells × %d genes = %d elements). "
            "Dense conversion for Scrublet may require %.1f GB of RAM. "
            "Options: "
            "(1) use --no-doublet-detection if memory is limited (results will be incomplete), "
            "(2) pre-filter to fewer genes before running AnnQC, "
            "(3) subsample to <50,000 cells for doublet detection. "
            "Future versions of AnnQC will support memory-efficient doublet detection.",
            adata.n_obs, adata.n_vars, _elements, _estimated_gb,
        )

    resolved_threshold: float = float("nan")

    scrublet_failed = False
    failure_reason = ""

    try:
        import scipy.sparse as sp
        import scrublet as scr

        counts = adata.X
        if sp.issparse(counts):
            counts = counts.toarray()
        counts = counts.astype(np.float32)

        n_components = min(30, n_cells - 1, counts.shape[1] - 1)

        scrub = scr.Scrublet(
            counts_matrix=counts,
            random_state=seed,
            sim_doublet_ratio=simulate_doublet_ratio,
        )

        doublet_scores, predicted_doublets = scrub.scrub_doublets(
            min_counts=2,
            min_cells=3,
            min_gene_variability_pctl=85,
            n_prin_comps=n_components,
            verbose=False,
        )

        if threshold != "auto":
            thresh_val = float(threshold)
            predicted_doublets = doublet_scores >= thresh_val
            resolved_threshold = thresh_val
        else:
            resolved_threshold = float(scrub.threshold_)

        adata.obs["annqc_doublet_score"] = doublet_scores.astype(float)
        adata.obs["annqc_is_doublet"] = predicted_doublets.astype(bool)

        n_doublets = int(predicted_doublets.sum())
        logger.info(
            "Scrublet complete: %d doublets (%.1f%%), threshold=%.3f",
            n_doublets,
            100.0 * n_doublets / max(n_cells, 1),
            resolved_threshold,
        )

        if "annqc" in adata.uns:
            adata.uns["annqc"]["doublet_status"] = "PASS"

    except Exception as exc:
        scrublet_failed = True
        failure_reason = f"{type(exc).__name__}: {exc}"
        logger.warning(
            "Scrublet failed (%s) — all cells marked as non-doublets",
            failure_reason,
        )
        adata.obs["annqc_doublet_score"] = np.nan
        adata.obs["annqc_is_doublet"] = False

        if "annqc" in adata.uns:
            adata.uns["annqc"]["doublet_status"] = "FAILED"
            adata.uns["annqc"]["doublet_failure_reason"] = failure_reason
            adata.uns["annqc"].setdefault("warnings", []).append(
                f"⚠️ CRITICAL: Doublet detection failed — {failure_reason}. "
                "All cells are marked as non-doublets. "
                "Status set to INCOMPLETE. "
                "Do not publish without manual doublet inspection."
            )
            # Downgrade PASS → INCOMPLETE; never let a failed doublet run be PASS
            if adata.uns["annqc"].get("status") == "PASS":
                adata.uns["annqc"]["status"] = "INCOMPLETE"

    if "annqc" in adata.uns and "thresholds" in adata.uns["annqc"]:
        adata.uns["annqc"]["thresholds"]["doublet_threshold"] = resolved_threshold

    return adata
