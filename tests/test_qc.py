"""Tests for annqc.qc.calculate_qc_metrics."""

import numpy as np
import pytest

from annqc.qc import calculate_qc_metrics


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _run_qc(adata):
    """Run calculate_qc_metrics with default mito/ribo prefixes and return result."""
    return calculate_qc_metrics(adata, mito_prefix="MT-", ribo_prefix="RPS|RPL")


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_qc_metrics_columns_added(small_adata):
    """Test that calculate_qc_metrics adds all required obs columns.

    The function must write at minimum: n_genes_by_counts, total_counts,
    and pct_counts_mt into adata.obs.
    """
    adata = _run_qc(small_adata)

    required_obs_cols = [
        "n_genes_by_counts",
        "total_counts",
        "pct_counts_mt",
    ]
    for col in required_obs_cols:
        assert col in adata.obs.columns, f"Missing obs column: {col}"


def test_qc_metrics_mito_annotation(small_adata):
    """Test that MT- genes are annotated as mito in var['mt'].

    After calculate_qc_metrics, adata.var['mt'] must exist and be True for
    every gene whose name starts with 'MT-'.
    """
    adata = _run_qc(small_adata)

    assert "mt" in adata.var.columns, "adata.var['mt'] column is missing"

    mt_genes = [g for g in adata.var_names if g.startswith("MT-")]
    assert len(mt_genes) > 0, "No MT- genes found — fixture problem"

    for gene in mt_genes:
        assert adata.var.loc[gene, "mt"], f"Gene {gene} should be flagged as mito"

    non_mt_genes = [g for g in adata.var_names if not g.startswith("MT-")]
    for gene in non_mt_genes[:5]:  # spot-check first 5
        assert not adata.var.loc[gene, "mt"], f"Gene {gene} should NOT be flagged as mito"


def test_qc_metrics_pct_counts_mt_range(small_adata):
    """Test that pct_counts_mt values are between 0 and 100.

    A percentage outside [0, 100] indicates a calculation error.
    """
    adata = _run_qc(small_adata)

    pct = adata.obs["pct_counts_mt"]
    assert (pct >= 0).all(), "pct_counts_mt has values below 0"
    assert (pct <= 100).all(), "pct_counts_mt has values above 100"


def test_qc_metrics_n_genes_positive(small_adata):
    """Test that n_genes_by_counts is positive for all cells.

    Every cell in the fixture has non-zero counts, so n_genes_by_counts must
    be >= 1 for all cells.
    """
    adata = _run_qc(small_adata)

    n_genes = adata.obs["n_genes_by_counts"]
    assert (n_genes > 0).all(), "Some cells have n_genes_by_counts == 0"


def test_qc_metrics_ribo_annotation(small_adata_with_ribo):
    """Test that RPS/RPL genes are annotated as ribo in var['ribo'].

    The fixture includes genes prefixed with 'RPS' and 'RPL'.  After running
    calculate_qc_metrics these must be flagged in adata.var['ribo'] and the
    corresponding pct_counts_ribo column must be present in adata.obs.
    """
    adata = calculate_qc_metrics(
        small_adata_with_ribo, mito_prefix="MT-", ribo_prefix="RPS|RPL"
    )

    assert "ribo" in adata.var.columns, "adata.var['ribo'] column is missing"

    ribo_genes = [g for g in adata.var_names if g.startswith("RPS") or g.startswith("RPL")]
    assert len(ribo_genes) > 0, "No RPS/RPL genes found — fixture problem"

    for gene in ribo_genes:
        assert adata.var.loc[gene, "ribo"], f"Gene {gene} should be flagged as ribo"

    assert "pct_counts_ribo" in adata.obs.columns, (
        "pct_counts_ribo should be added to obs when ribo genes are present"
    )
    pct = adata.obs["pct_counts_ribo"]
    assert (pct >= 0).all()
    assert (pct <= 100).all()


def test_qc_metrics_total_counts_equals_row_sum(small_adata):
    """Test that total_counts matches the row-wise sum of the count matrix.

    This is a correctness sanity-check: total_counts must equal numpy's sum
    of each cell's raw count vector.
    """
    adata = _run_qc(small_adata)

    expected = np.asarray(small_adata.X).sum(axis=1)
    actual = adata.obs["total_counts"].values
    np.testing.assert_allclose(actual, expected, rtol=1e-5)


def test_qc_metrics_does_not_modify_X(small_adata):
    """Test that calculate_qc_metrics does not alter the count matrix.

    The raw count matrix X must be identical before and after the call.
    """
    import copy

    original_X = small_adata.X.copy()
    _run_qc(small_adata)

    np.testing.assert_array_equal(
        np.asarray(small_adata.X),
        np.asarray(original_X),
        err_msg="calculate_qc_metrics must not modify adata.X",
    )


def test_qc_metrics_high_mito_cells_have_higher_pct(small_adata):
    """Test that cells with artificially boosted MT counts rank highest in pct_counts_mt.

    Cells 0-9 in the fixture had their mito counts multiplied by 10.  After
    calculate_qc_metrics these cells should have higher pct_counts_mt than the
    median of the rest.
    """
    adata = _run_qc(small_adata)

    high_mito_pct = adata.obs["pct_counts_mt"].iloc[:10].mean()
    rest_pct = adata.obs["pct_counts_mt"].iloc[10:].median()

    assert high_mito_pct > rest_pct, (
        f"High-mito cells (mean {high_mito_pct:.1f}%) should exceed "
        f"rest median ({rest_pct:.1f}%)"
    )
