"""Cluster-aware mitochondrial QC warnings for AnnQC."""

import logging

logger = logging.getLogger(__name__)


def run_cluster_qc(
    adata,
    mito_threshold: float,
    resolution: float = 0.3,
    n_comps: int = 20,
):
    """Cluster cells and flag clusters with biology-driven high mito content.

    Runs PCA + Leiden clustering on a normalized copy of adata, then copies
    cluster assignments back. Raises warnings when a cluster's median mito%
    exceeds half the filtering threshold, or when >50% of cells in a cluster
    would be removed.

    Parameters
    ----------
    adata : AnnData
        Raw (unfiltered) AnnData. Must have pct_counts_mt in obs.
    mito_threshold : float
        The mito max_pct that will be applied to filter cells.
    resolution : float
        Leiden resolution. 0.3 produces coarse clusters (5-15 typically).
    n_comps : int
        Number of PCA components.

    Returns
    -------
    adata : AnnData
        Input adata with annqc_cluster added to obs (int string labels).
    cluster_stats : dict
        {cluster_id: {n_cells, median_mito, pct_would_remove}}.
    warnings : list of str
        Human-readable warnings for flagged clusters.
    """
    import scanpy as sc

    if "pct_counts_mt" not in adata.obs.columns:
        logger.warning("pct_counts_mt not in adata.obs — cluster QC skipped")
        return adata, {}, []

    if adata.n_obs < 50:
        logger.warning("Too few cells (%d) for cluster QC — skipped", adata.n_obs)
        return adata, {}, []

    logger.info(
        "Running cluster QC on %d cells (resolution=%.2f, n_comps=%d)",
        adata.n_obs, resolution, n_comps,
    )

    copy_ = adata.copy()

    sc.pp.normalize_total(copy_, target_sum=10_000)
    sc.pp.log1p(copy_)

    n_comps_actual = min(n_comps, copy_.n_obs - 1, copy_.n_vars - 1)
    sc.pp.pca(copy_, n_comps=n_comps_actual)

    # Try Leiden first (requires igraph or leidenalg). Fall back to k-means.
    _cluster_method = "leiden"
    try:
        n_pcs = min(n_comps_actual, copy_.obsm["X_pca"].shape[1])
        sc.pp.neighbors(copy_, n_pcs=n_pcs, random_state=0)
        sc.tl.leiden(
            copy_,
            resolution=resolution,
            key_added="annqc_cluster",
            random_state=0,
            flavor="leidenalg",
        )
        cluster_labels = copy_.obs["annqc_cluster"].values
    except (ImportError, TypeError):
        try:
            sc.tl.leiden(
                copy_,
                resolution=resolution,
                key_added="annqc_cluster",
                random_state=0,
            )
            cluster_labels = copy_.obs["annqc_cluster"].values
        except (ImportError, TypeError):
            from sklearn.cluster import KMeans
            n_clusters = max(5, min(15, copy_.n_obs // 200))
            km = KMeans(n_clusters=n_clusters, random_state=0, n_init=10)
            cluster_labels = km.fit_predict(copy_.obsm["X_pca"][:, :n_comps_actual]).astype(str)
            _cluster_method = "kmeans"
            logger.info(
                "Leiden not available (igraph/leidenalg not installed). "
                "Falling back to KMeans (%d clusters) for cluster QC.",
                n_clusters,
            )

    adata.obs["annqc_cluster"] = cluster_labels
    del copy_
    logger.debug("Cluster QC method: %s", _cluster_method)

    cluster_stats = {}
    warnings = []

    clusters = sorted(adata.obs["annqc_cluster"].unique(), key=lambda x: int(x))
    for cluster in clusters:
        mask = adata.obs["annqc_cluster"] == cluster
        cluster_mito = adata.obs.loc[mask, "pct_counts_mt"]
        n_cells = int(mask.sum())
        median_mito = float(cluster_mito.median())
        pct_would_remove = float((cluster_mito > mito_threshold).mean()) * 100

        cluster_stats[str(cluster)] = {
            "n_cells": n_cells,
            "median_mito": round(median_mito, 2),
            "pct_would_remove": round(pct_would_remove, 1),
        }
        logger.debug(
            "Cluster %s: %d cells, median_mito=%.1f%%, pct_removed=%.0f%%",
            cluster, n_cells, median_mito, pct_would_remove,
        )

        if median_mito > mito_threshold / 2:
            warnings.append(
                f"Cluster {cluster} ({n_cells:,} cells): median mito {median_mito:.1f}% "
                f"exceeds half the filtering threshold ({mito_threshold / 2:.1f}%). "
                f"Verify this is biology, not damage, before accepting the mito threshold."
            )

        if pct_would_remove > 50:
            warnings.append(
                f"Cluster {cluster} ({n_cells:,} cells): {pct_would_remove:.0f}% of cells "
                f"would be removed by the mito threshold ({mito_threshold}%). "
                f"If this is a known high-mito cell type (e.g. cardiomyocytes), "
                f"consider raising mito.max_pct or using a tissue preset."
            )

    logger.info(
        "Cluster QC complete: %d clusters, %d warning(s)", len(cluster_stats), len(warnings),
    )
    return adata, cluster_stats, warnings
