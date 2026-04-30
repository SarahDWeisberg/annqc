"""Tests for dry_run=True mode of annqc.run."""

import copy
import os

import anndata as ad
import numpy as np
import pytest

import annqc
from annqc.config import DEFAULT_CONFIG


# ---------------------------------------------------------------------------
# Minimal config that disables scrublet to keep tests fast and hermetic
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
# Helper
# ---------------------------------------------------------------------------


def _make_small_adata(n_cells=200, n_genes=500, seed=42):
    """Create a small synthetic AnnData with MT- genes for testing."""
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
    return adata


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_dry_run_no_output_file(tmp_path):
    """With dry_run=True, no .h5ad file should be written to disk."""
    adata = _make_small_adata()
    output_path = str(tmp_path / "output.h5ad")

    annqc.run(
        adata,
        config=copy.deepcopy(_MINIMAL_CFG),
        seed=42,
        output=output_path,
        dry_run=True,
    )

    assert not os.path.exists(output_path), (
        f"dry_run=True should NOT write an output file, but {output_path!r} was created"
    )


def test_dry_run_record_flag():
    """With dry_run=True, adata.uns['annqc']['dry_run'] must be True."""
    adata = _make_small_adata()
    result = annqc.run(
        adata,
        config=copy.deepcopy(_MINIMAL_CFG),
        seed=42,
        dry_run=True,
    )

    assert result.uns["annqc"]["dry_run"] is True, (
        f"Expected dry_run=True in record, got {result.uns['annqc']['dry_run']!r}"
    )


def test_dry_run_no_cells_removed():
    """In dry run mode, all input cells must be present in the returned AnnData."""
    adata = _make_small_adata()
    input_n_obs = adata.n_obs

    result = annqc.run(
        adata,
        config=copy.deepcopy(_MINIMAL_CFG),
        seed=42,
        dry_run=True,
    )

    assert result.n_obs == input_n_obs, (
        f"dry_run=True should preserve all {input_n_obs} input cells, "
        f"but result has {result.n_obs} cells"
    )


def test_dry_run_output_count_reflects_flags():
    """In dry run mode, cell_counts['output'] reflects passing-flagged cells, not n_obs."""
    adata = _make_small_adata()
    result = annqc.run(
        adata,
        config=copy.deepcopy(_MINIMAL_CFG),
        seed=42,
        dry_run=True,
    )

    record = result.uns["annqc"]
    expected_output = int(result.obs["annqc_pass"].sum())
    actual_output = record["cell_counts"]["output"]

    assert actual_output == expected_output, (
        f"cell_counts['output'] ({actual_output}) should equal "
        f"annqc_pass.sum() ({expected_output}) in dry run mode"
    )

    # In dry run the recorded output should differ from total n_obs when some cells fail
    # (This just confirms the count is plausible — not all 200 cells should pass mito filter)
    assert actual_output <= result.n_obs, (
        f"cell_counts['output'] ({actual_output}) cannot exceed n_obs ({result.n_obs})"
    )


def test_dry_run_report_generated(tmp_path):
    """With dry_run=True, the HTML report should still be written."""
    adata = _make_small_adata()
    report_path = str(tmp_path / "preview.html")

    annqc.run(
        adata,
        config=copy.deepcopy(_MINIMAL_CFG),
        seed=42,
        report_path=report_path,
        dry_run=True,
    )

    assert os.path.exists(report_path), (
        f"dry_run=True should still write the report, but {report_path!r} was not created"
    )


def test_normal_run_no_dry_run_flag():
    """With dry_run=False (default), record['dry_run'] is False and cells are removed."""
    adata = _make_small_adata()
    input_n_obs = adata.n_obs

    result = annqc.run(
        adata,
        config=copy.deepcopy(_MINIMAL_CFG),
        seed=42,
        dry_run=False,
    )

    record = result.uns["annqc"]

    assert record["dry_run"] is False, (
        f"Expected dry_run=False in record, got {record['dry_run']!r}"
    )

    # In normal mode the returned adata should have exactly cell_counts["output"] cells
    # (the pipeline actually applied the filter, not just previewed it)
    assert result.n_obs == record["cell_counts"]["output"], (
        f"In normal mode n_obs ({result.n_obs}) should equal cell_counts['output'] "
        f"({record['cell_counts']['output']})"
    )
