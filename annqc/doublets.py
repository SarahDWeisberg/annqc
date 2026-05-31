"""Scrublet-based doublet detection for AnnQC."""

import logging

import numpy as np

logger = logging.getLogger(__name__)


def detect_doublets(
    adata,
    threshold="auto",
    simulate_doublet_ratio: float = 2.0,
    seed: int = 0,
    subsample_n: int = 0,
):
    """Run Scrublet doublet detection on adata.

    Adds to adata.obs:
        annqc_doublet_score (float) — scrublet doublet probability score
        annqc_is_doublet    (bool)  — True if cell is predicted doublet

    If scrublet fails for any reason, logs a warning and falls back to
    NaN scores and False predictions so the pipeline can continue.

    The resolved doublet threshold (float) is stored in
    adata.uns['annqc']['thresholds']['doublet_threshold'] if the key exists.

    Parameters
    ----------
    subsample_n : int
        If > 0 and n_cells > subsample_n, run Scrublet on a random subsample
        of subsample_n cells. Unsampled cells receive NaN scores and are
        marked as non-doublets. Use to avoid out-of-memory errors on large
        datasets. 0 means no subsampling.
    """
    n_cells = adata.n_obs
    logger.info("Running Scrublet on %d cells (seed=%d, threshold=%s)", n_cells, seed, threshold)

    # --- Subsampling ---
    _subsampled = False
    _subsample_indices = None
    if subsample_n > 0 and n_cells > subsample_n:
        rng = np.random.default_rng(seed)
        _subsample_indices = rng.choice(n_cells, size=subsample_n, replace=False)
        _subsampled = True
        logger.warning(
            "Subsampling %d → %d cells for doublet detection (memory-safe mode). "
            "%d cells will not be scored and will be marked as non-doublets.",
            n_cells, subsample_n, n_cells - subsample_n,
        )
        if "annqc" in adata.uns:
            adata.uns["annqc"]["doublet_subsampled"] = True
            adata.uns["annqc"]["doublet_subsample_n"] = subsample_n

    _elements = adata.n_obs * adata.n_vars
    if _elements > 500_000_000 and not _subsampled:
        _estimated_gb = (_elements * 4) / 1e9
        logger.warning(
            "Large matrix detected (%d cells × %d genes = %d elements). "
            "Dense conversion for Scrublet may require %.1f GB of RAM. "
            "Use --subsample-doublets 50000 to score a representative subset.",
            adata.n_obs, adata.n_vars, _elements, _estimated_gb,
        )

    resolved_threshold: float = float("nan")

    scrublet_failed = False
    failure_reason = ""

    try:
        import scipy.sparse as sp
        import scrublet as scr

        if _subsampled:
            counts_full = adata.X
            if sp.issparse(counts_full):
                counts_input = counts_full[_subsample_indices].toarray()
            else:
                counts_input = counts_full[_subsample_indices]
        else:
            counts_input = adata.X
            if sp.issparse(counts_input):
                counts_input = counts_input.toarray()
        counts_input = counts_input.astype(np.float32)

        n_scored = counts_input.shape[0]
        n_components = min(30, n_scored - 1, counts_input.shape[1] - 1)

        scrub = scr.Scrublet(
            counts_matrix=counts_input,
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

        # Expand scores to full cell array (unsampled cells get NaN / False)
        if _subsampled:
            full_scores = np.full(n_cells, np.nan)
            full_doublets = np.zeros(n_cells, dtype=bool)
            full_scores[_subsample_indices] = doublet_scores.astype(float)
            full_doublets[_subsample_indices] = predicted_doublets.astype(bool)
        else:
            full_scores = doublet_scores.astype(float)
            full_doublets = predicted_doublets.astype(bool)

        adata.obs["annqc_doublet_score"] = full_scores
        adata.obs["annqc_is_doublet"] = full_doublets

        n_doublets = int(full_doublets.sum())
        logger.info(
            "Scrublet complete: %d doublets (%.1f%% of %d scored cells), threshold=%.3f",
            n_doublets,
            100.0 * n_doublets / max(n_scored, 1),
            n_scored,
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
                "Do not publish without manual doublet inspection."
            )

    if "annqc" in adata.uns and "thresholds" in adata.uns["annqc"]:
        adata.uns["annqc"]["thresholds"]["doublet_threshold"] = resolved_threshold

    return adata
