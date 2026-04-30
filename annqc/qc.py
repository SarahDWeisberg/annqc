"""QC metric calculation for AnnQC (wraps scanpy)."""

import logging

import scanpy as sc

logger = logging.getLogger(__name__)


def calculate_qc_metrics(adata, mito_prefix: str = "MT-", ribo_prefix: str = "RPS|RPL"):
    """Calculate QC metrics and annotate mitochondrial / ribosomal genes.

    Adds to adata.var:
        mt   (bool) — True if gene name starts with mito_prefix
        ribo (bool) — True if gene name matches ribo_prefix regex

    Adds to adata.obs:
        pct_counts_mt     — % of counts from MT genes
        pct_counts_ribo   — % of counts from ribosomal genes
        n_genes_by_counts — number of genes with at least 1 count
        total_counts      — total UMI counts per cell

    Uses scanpy.pp.calculate_qc_metrics internally.
    """
    logger.info("Calculating QC metrics (mito_prefix=%r, ribo_prefix=%r)", mito_prefix, ribo_prefix)

    adata.var["mt"] = adata.var_names.str.startswith(mito_prefix)
    adata.var["ribo"] = adata.var_names.str.contains(ribo_prefix, regex=True)

    n_mito = int(adata.var["mt"].sum())
    n_ribo = int(adata.var["ribo"].sum())
    logger.debug("Annotated %d mito genes, %d ribo genes", n_mito, n_ribo)

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    logger.info(
        "QC metrics calculated for %d cells: "
        "pct_counts_mt, pct_counts_ribo, n_genes_by_counts, total_counts",
        adata.n_obs,
    )
    return adata
