"""Tests for the _generate_warnings function and per-sample status logic."""

import copy

import anndata as ad
import numpy as np
import pytest

import annqc
from annqc.config import DEFAULT_CONFIG
from annqc.pipeline import _generate_warnings


# ---------------------------------------------------------------------------
# Minimal config for full-pipeline tests
# ---------------------------------------------------------------------------

_MINIMAL_CFG = {
    "mito": {"prefix": "MT-", "max_pct": 20},
    "ribo": {"prefix": "RPS|RPL", "max_pct": 50},
    "cells": {"min_genes": 100, "max_genes": 5000, "min_counts": 100, "max_counts": None},
    "genes": {"min_cells": 3},
    "doublets": {"method": "scrublet", "threshold": "auto", "simulate_doublet_ratio": 2.0},
    "normalization": {"method": "log1p", "target_sum": 10000},
    "report": {"title": "Test", "author": ""},
}


# ---------------------------------------------------------------------------
# Helper: build a minimal fake record for unit-testing _generate_warnings
# ---------------------------------------------------------------------------


def _make_record(
    n_input=2000,
    n_after_mito=2000,
    n_after_gene=2000,
    n_after_count=2000,
    n_after_doublet=2000,
    n_output=2000,
    per_sample=None,
):
    """Return a minimal fake record dict accepted by _generate_warnings."""
    return {
        "cell_counts": {
            "input": n_input,
            "after_mito_filter": n_after_mito,
            "after_gene_filter": n_after_gene,
            "after_count_filter": n_after_count,
            "after_doublet_filter": n_after_doublet,
            "output": n_output,
        },
        "warnings": [],
        "per_sample": per_sample or {},
    }


def _make_small_adata(n_cells=300, n_genes=500, seed=42, high_mito_cells=0):
    """Create a small synthetic AnnData for full-pipeline tests."""
    np.random.seed(seed)
    X = np.random.negative_binomial(5, 0.5, size=(n_cells, n_genes)).astype(float)

    var_names = [f"Gene{i}" for i in range(n_genes)]
    var_names[:20] = [f"MT-Gene{i}" for i in range(20)]

    obs_names = [f"Cell{i}" for i in range(n_cells)]

    adata = ad.AnnData(
        X=X,
        obs={"cell_id": obs_names},
        var={"gene_id": var_names},
    )
    adata.obs_names = obs_names
    adata.var_names = var_names

    if high_mito_cells > 0:
        # Inflate MT counts for first `high_mito_cells` cells
        adata.X[:high_mito_cells, :20] = adata.X[:high_mito_cells, :20] * 100

    return adata


# ---------------------------------------------------------------------------
# Tests: _generate_warnings in isolation
# ---------------------------------------------------------------------------


def test_warning_few_cells():
    """When output < 500 cells, a warning about remaining cells must be added."""
    cfg = {"mito": {"max_pct": 20}}
    record = _make_record(n_input=1000, n_output=300)

    _generate_warnings(record, cfg)

    warnings = record["warnings"]
    assert len(warnings) > 0, "Expected at least one warning for few remaining cells"
    assert any("cells remain" in w for w in warnings), (
        f"Expected a warning containing 'cells remain', got: {warnings}"
    )


def test_warning_high_doublet_rate():
    """When doublet rate > 8%, a doublet-rate warning must be added."""
    cfg = {"mito": {"max_pct": 20}}
    # 1000 cells after count filter, only 800 after doublet filter → 20% rate
    record = _make_record(
        n_input=1500,
        n_after_count=1000,
        n_after_doublet=800,
        n_output=800,
    )

    _generate_warnings(record, cfg)

    warnings = record["warnings"]
    assert any(
        "doublet" in w.lower() or "Doublet" in w for w in warnings
    ), f"Expected a doublet-rate warning, got: {warnings}"


def test_warning_mito_filter_aggressive():
    """When mito filter removes > 15% of cells, a mito warning must be added."""
    cfg = {"mito": {"max_pct": 20}}
    # 2000 input, only 1500 after mito → 25% removed
    record = _make_record(
        n_input=2000,
        n_after_mito=1500,
        n_output=1500,
    )

    _generate_warnings(record, cfg)

    warnings = record["warnings"]
    assert any(
        "mito" in w.lower() or "Mito" in w for w in warnings
    ), f"Expected a mito-filter warning, got: {warnings}"


def test_no_warnings_normal_run():
    """A clean dataset with healthy cell counts should produce no filter warnings."""
    cfg = {"mito": {"max_pct": 20}}
    # Healthy run: 2000 input, all filters remove very few cells
    record = _make_record(
        n_input=2000,
        n_after_mito=1950,
        n_after_gene=1940,
        n_after_count=1930,
        n_after_doublet=1920,
        n_output=1920,
    )

    _generate_warnings(record, cfg)

    warnings = record["warnings"]
    # Should have zero warnings (or none about aggressive filtering)
    filter_warnings = [
        w for w in warnings
        if any(kw in w for kw in ("cells remain", "Mito filter", "Doublet rate"))
    ]
    assert len(filter_warnings) == 0, (
        f"Expected no filter warnings for a clean run, got: {filter_warnings}"
    )


# ---------------------------------------------------------------------------
# Test: per-sample WARN status via full pipeline
# ---------------------------------------------------------------------------


def test_per_sample_warn_status():
    """A sample with high median mito% should have status == 'WARN' in per_sample."""
    np.random.seed(0)
    n_cells, n_genes = 400, 500

    X = np.random.negative_binomial(5, 0.5, size=(n_cells, n_genes)).astype(float)

    var_names = [f"Gene{i}" for i in range(n_genes)]
    var_names[:20] = [f"MT-Gene{i}" for i in range(20)]

    obs_names = [f"Cell{i}" for i in range(n_cells)]

    # Sample A: cells 0-199 are normal
    # Sample B: cells 200-399 have very high MT counts → high median mito%
    sample_labels = ["A"] * 200 + ["B"] * 200

    adata = ad.AnnData(
        X=X,
        obs={"cell_id": obs_names, "sample": sample_labels},
        var={"gene_id": var_names},
    )
    adata.obs_names = obs_names
    adata.var_names = var_names

    # Inflate MT counts for Sample B cells to drive mito% well above threshold
    adata.X[200:, :20] = adata.X[200:, :20] * 500

    cfg = copy.deepcopy(_MINIMAL_CFG)
    # Use a permissive mito filter so Sample B cells aren't all removed (just flagged)
    cfg["mito"]["max_pct"] = 100
    cfg["cells"]["min_counts"] = 10
    cfg["cells"]["min_genes"] = 10

    result = annqc.run(
        adata,
        config=cfg,
        seed=42,
        sample_key="sample",
    )

    record = result.uns["annqc"]
    per_sample = record.get("per_sample", {})

    assert "B" in per_sample, (
        f"Sample 'B' missing from per_sample. Keys: {list(per_sample.keys())}"
    )
    assert per_sample["B"]["status"] == "WARN", (
        f"Expected sample 'B' to have status='WARN' due to high mito%, "
        f"got {per_sample['B']['status']!r}. median_mito_pct={per_sample['B'].get('median_mito_pct')}"
    )
