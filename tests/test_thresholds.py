"""Tests for annqc.thresholds.suggest_thresholds."""

import copy

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp

from annqc.qc import calculate_qc_metrics
from annqc.thresholds import suggest_thresholds


# ---------------------------------------------------------------------------
# Fixture
# ---------------------------------------------------------------------------


@pytest.fixture
def adata_with_qc():
    """Create a small synthetic AnnData with QC metrics already calculated."""
    np.random.seed(7)
    n_cells, n_genes = 200, 500

    X = np.random.negative_binomial(5, 0.5, size=(n_cells, n_genes)).astype(float)

    var_names = [f"Gene{i}" for i in range(n_genes)]
    # Make first 20 genes mitochondrial and genes 20-29 ribosomal
    var_names[:20] = [f"MT-Gene{i}" for i in range(20)]
    var_names[20:30] = [f"RPS{i}" for i in range(10)]

    obs_names = [f"Cell{i}" for i in range(n_cells)]

    adata = ad.AnnData(
        X=X,
        obs={"cell_id": obs_names},
        var={"gene_id": var_names},
    )
    adata.obs_names = obs_names
    adata.var_names = var_names

    adata = calculate_qc_metrics(adata, mito_prefix="MT-", ribo_prefix="RPS|RPL")
    return adata


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_suggest_thresholds_returns_all_levels(adata_with_qc):
    """suggest_thresholds must return a dict with all four top-level keys."""
    result = suggest_thresholds(adata_with_qc)
    for key in ("strict", "standard", "permissive", "raw"):
        assert key in result, f"Missing top-level key in result: '{key}'"


def test_suggest_thresholds_standard_keys(adata_with_qc):
    """The 'standard' level must contain all five expected threshold keys."""
    result = suggest_thresholds(adata_with_qc)
    standard = result["standard"]
    for key in ("mito_max_pct", "min_genes", "max_genes", "min_counts", "max_counts"):
        assert key in standard, f"Missing key in result['standard']: '{key}'"


def test_suggest_thresholds_strict_lt_permissive(adata_with_qc):
    """For mito_max_pct, strict < standard < permissive (stricter = lower cap)."""
    result = suggest_thresholds(adata_with_qc)
    strict_mito = result["strict"]["mito_max_pct"]
    standard_mito = result["standard"]["mito_max_pct"]
    permissive_mito = result["permissive"]["mito_max_pct"]

    # All three should be non-None since MT- genes are present
    assert strict_mito is not None
    assert standard_mito is not None
    assert permissive_mito is not None

    assert strict_mito < standard_mito, (
        f"strict mito_max_pct ({strict_mito}) should be < standard ({standard_mito})"
    )
    assert standard_mito < permissive_mito, (
        f"standard mito_max_pct ({standard_mito}) should be < permissive ({permissive_mito})"
    )


def test_suggest_thresholds_lower_bounds_nonnegative(adata_with_qc):
    """All threshold values must be >= 0 (or None when metric is absent)."""
    result = suggest_thresholds(adata_with_qc)
    for level_name in ("strict", "standard", "permissive"):
        level = result[level_name]
        for key, value in level.items():
            if value is not None:
                assert value >= 0.0, (
                    f"result['{level_name}']['{key}'] = {value} is negative"
                )


def test_suggest_thresholds_mad_formula(adata_with_qc):
    """Verify the MAD threshold formula on a known-value column.

    For n_genes_by_counts:
        lower = max(0, median - 5 * MAD)
        upper = median + 5 * MAD
    """
    vals = adata_with_qc.obs["n_genes_by_counts"].values.astype(float)
    median = float(np.median(vals))
    mad = float(np.median(np.abs(vals - median)))

    expected_lower = max(0.0, median - 5 * mad)
    expected_upper = median + 5 * mad

    result = suggest_thresholds(adata_with_qc)
    standard = result["standard"]

    assert abs(standard["min_genes"] - expected_lower) < 1e-9, (
        f"min_genes expected {expected_lower}, got {standard['min_genes']}"
    )
    assert abs(standard["max_genes"] - expected_upper) < 1e-9, (
        f"max_genes expected {expected_upper}, got {standard['max_genes']}"
    )


def test_suggest_thresholds_raw_has_median_mad(adata_with_qc):
    """result['raw']['pct_counts_mt'] must have 'median' and 'mad' with sensible values."""
    result = suggest_thresholds(adata_with_qc)
    raw = result["raw"]

    assert "pct_counts_mt" in raw, "pct_counts_mt missing from result['raw']"
    entry = raw["pct_counts_mt"]
    assert "median" in entry, "Missing 'median' in raw['pct_counts_mt']"
    assert "mad" in entry, "Missing 'mad' in raw['pct_counts_mt']"

    # Median mito% must be a non-negative float
    assert isinstance(entry["median"], float), (
        f"raw['pct_counts_mt']['median'] should be float, got {type(entry['median'])}"
    )
    assert entry["median"] >= 0.0, (
        f"raw['pct_counts_mt']['median'] should be >= 0, got {entry['median']}"
    )
    # MAD must be non-negative
    assert entry["mad"] >= 0.0, (
        f"raw['pct_counts_mt']['mad'] should be >= 0, got {entry['mad']}"
    )


def test_suggest_thresholds_missing_metric():
    """suggest_thresholds must not crash when QC columns are absent; values should be None."""
    np.random.seed(1)
    n_cells, n_genes = 50, 100
    X = np.random.negative_binomial(5, 0.5, size=(n_cells, n_genes)).astype(float)

    # Create adata with no QC metrics in obs
    adata = ad.AnnData(X=X)
    adata.obs_names = [f"Cell{i}" for i in range(n_cells)]
    adata.var_names = [f"Gene{i}" for i in range(n_genes)]

    # Should not raise
    result = suggest_thresholds(adata)

    # Since no metrics exist, all threshold values should be None
    for level_name in ("strict", "standard", "permissive"):
        level = result[level_name]
        for key, value in level.items():
            assert value is None, (
                f"result['{level_name}']['{key}'] should be None when metric is absent, "
                f"got {value}"
            )

    # raw should be empty
    assert result["raw"] == {}, f"raw should be empty dict, got {result['raw']}"
