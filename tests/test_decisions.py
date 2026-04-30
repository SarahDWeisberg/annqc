"""Tests for annqc.decisions — explain_threshold, analyze_distributions."""

import math

import anndata as ad
import numpy as np
import pytest

from annqc.decisions import analyze_distributions, explain_threshold


# ── explain_threshold ────────────────────────────────────────────────────────

def test_explain_threshold_no_data():
    result = explain_threshold("pct_counts_mt", [], None, None)
    assert result == "No data available."


def test_explain_threshold_manual_high_only():
    vals = list(np.random.default_rng(0).uniform(0, 30, 200))
    result = explain_threshold("pct_counts_mt", vals, None, 20.0, method="manual")
    assert "20" in result
    assert "manual" in result
    assert "%" in result


def test_explain_threshold_manual_low_only():
    vals = list(np.arange(100, 5100, 25, dtype=float))  # 200 values
    result = explain_threshold("n_genes_by_counts", vals, 200.0, None, method="manual")
    assert "200" in result
    assert "manual" in result


def test_explain_threshold_auto_mad_high():
    rng = np.random.default_rng(1)
    vals = list(rng.normal(2000, 500, 300))
    result = explain_threshold("total_counts", vals, None, 3500.0, method="auto_mad")
    assert "MAD-based" in result
    assert "3500" in result or "3.5e+03" in result or "3.5" in result


def test_explain_threshold_auto_mad_low():
    rng = np.random.default_rng(2)
    vals = list(rng.normal(2000, 500, 300))
    result = explain_threshold("n_genes_by_counts", vals, 500.0, 4000.0, method="auto_mad")
    assert "MAD-based" in result
    assert "500" in result


def test_explain_threshold_no_threshold():
    vals = list(np.arange(10, 60, dtype=float))
    result = explain_threshold("pct_counts_ribo", vals, 0.0, None)
    assert "No threshold applied" in result


def test_explain_threshold_pct_below_reported():
    # 80 values below 20, 20 values above → 80% below
    vals = list(np.linspace(5, 19.9, 80)) + list(np.linspace(21, 40, 20))
    result = explain_threshold("pct_counts_mt", vals, None, 20.0, method="manual")
    assert "%" in result
    assert "80.0" in result


# ── analyze_distributions ────────────────────────────────────────────────────

def _make_adata_with_record(n=300):
    """Build a minimal AnnData with adata.uns['annqc'] populated."""
    rng = np.random.default_rng(0)
    adata = ad.AnnData(X=rng.poisson(5, size=(n, 50)).astype(float))
    mt_vals = list(rng.uniform(0, 25, n))
    gene_vals = list(rng.integers(200, 4000, n).astype(float))
    count_vals = list(rng.integers(500, 8000, n).astype(float))
    adata.obs["pct_counts_mt"] = mt_vals
    adata.obs["n_genes_by_counts"] = gene_vals
    adata.obs["total_counts"] = count_vals
    adata.uns["annqc"] = {
        "raw_obs_metrics": {
            "pct_counts_mt": mt_vals,
            "n_genes_by_counts": gene_vals,
            "total_counts": count_vals,
        },
        "thresholds": {
            "mito_max_pct": 20.0,
            "min_genes": 300.0,
            "max_genes": 3500.0,
            "min_counts": 600.0,
            "max_counts": None,
            "doublet_threshold": 0.25,
        },
        "threshold_method": "manual",
        "threshold_explanations": {},
        "cell_counts": {"input": n, "after_doublet_filter": n - 10},
    }
    return adata


def test_analyze_distributions_populates_record():
    adata = _make_adata_with_record()
    analyze_distributions(adata)
    record = adata.uns["annqc"]
    assert "pct_counts_mt" in record["threshold_explanations"]
    assert "n_genes_by_counts" in record["threshold_explanations"]
    assert "total_counts" in record["threshold_explanations"]


def test_analyze_distributions_doublet_explanation():
    adata = _make_adata_with_record()
    analyze_distributions(adata)
    record = adata.uns["annqc"]
    assert "doublet" in record["threshold_explanations"]
    assert "Scrublet" in record["threshold_explanations"]["doublet"]


def test_analyze_distributions_missing_metrics():
    rng = np.random.default_rng(0)
    adata = ad.AnnData(X=rng.poisson(5, size=(50, 20)).astype(float))
    adata.uns["annqc"] = {
        "raw_obs_metrics": {},
        "thresholds": {"mito_max_pct": 20.0},
        "threshold_method": "manual",
        "threshold_explanations": {},
        "cell_counts": {"input": 50, "after_doublet_filter": 48},
    }
    analyze_distributions(adata)
    record = adata.uns["annqc"]
    assert record["threshold_explanations"] == {} or "doublet" in record["threshold_explanations"]
