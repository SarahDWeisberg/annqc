"""Tests for annqc.doublets.detect_doublets."""

import copy
from unittest.mock import patch

import numpy as np
import pytest

import annqc
from annqc.config import DEFAULT_CONFIG
from annqc.doublets import detect_doublets
from annqc.qc import calculate_qc_metrics


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _with_qc(adata):
    """Return adata with QC metrics pre-calculated (some implementations need this)."""
    return calculate_qc_metrics(adata, mito_prefix="MT-", ribo_prefix="RPS|RPL")


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_detect_doublets_adds_obs_columns(small_adata):
    """Test that detect_doublets adds annqc_doublet_score and annqc_is_doublet.

    Both columns are required by the per-cell contract regardless of whether
    any doublets are actually detected.
    """
    adata = detect_doublets(small_adata.copy(), seed=42)

    assert "annqc_doublet_score" in adata.obs.columns, (
        "annqc_doublet_score missing from obs after detect_doublets"
    )
    assert "annqc_is_doublet" in adata.obs.columns, (
        "annqc_is_doublet missing from obs after detect_doublets"
    )


def test_detect_doublets_score_range(small_adata):
    """Test that doublet scores are between 0 and 1 (or NaN).

    Scrublet scores represent a probability-like value in [0, 1].  Any score
    outside this range indicates a normalisation error.  NaN is allowed for
    edge-case cells where a score cannot be computed.
    """
    adata = detect_doublets(small_adata.copy(), seed=42)

    scores = adata.obs["annqc_doublet_score"]
    valid_mask = scores.notna()

    assert (scores[valid_mask] >= 0.0).all(), "Some doublet scores are below 0"
    assert (scores[valid_mask] <= 1.0).all(), "Some doublet scores are above 1"


def test_detect_doublets_is_doublet_is_bool(small_adata):
    """Test that annqc_is_doublet contains only bool values.

    Downstream logic uses boolean indexing on this column; a non-bool dtype
    would silently corrupt cell removal steps.
    """
    adata = detect_doublets(small_adata.copy(), seed=42)

    dtype = adata.obs["annqc_is_doublet"].dtype
    assert dtype == bool or str(dtype) in ("bool", "boolean"), (
        f"annqc_is_doublet dtype should be bool, got {dtype}"
    )


def test_detect_doublets_seed_deterministic(small_adata):
    """Test that the same seed produces the same doublet scores.

    Reproducibility is a core requirement; two runs with identical seeds on
    the same data must yield bit-identical doublet score vectors.
    """
    adata1 = detect_doublets(small_adata.copy(), seed=42)
    adata2 = detect_doublets(small_adata.copy(), seed=42)

    scores1 = adata1.obs["annqc_doublet_score"].values
    scores2 = adata2.obs["annqc_doublet_score"].values

    np.testing.assert_array_equal(
        scores1,
        scores2,
        err_msg="detect_doublets with the same seed produced different scores",
    )


def test_detect_doublets_different_seeds_differ(small_adata):
    """Test that different seeds produce different doublet score distributions.

    While not guaranteed in the strictest sense, two different seeds should
    almost always produce at least one differing score in a 200-cell dataset.
    """
    adata1 = detect_doublets(small_adata.copy(), seed=0)
    adata2 = detect_doublets(small_adata.copy(), seed=99)

    scores1 = adata1.obs["annqc_doublet_score"].values
    scores2 = adata2.obs["annqc_doublet_score"].values

    # If all scores happen to be identical with different seeds something is wrong
    assert not np.array_equal(scores1, scores2), (
        "detect_doublets produced identical scores with seed=0 and seed=99 — "
        "the seed parameter may not be respected"
    )


def test_detect_doublets_does_not_change_n_obs(small_adata):
    """Test that detect_doublets does not remove cells.

    detect_doublets is an annotation step; it must not alter the cell count.
    Cell removal happens later in apply_filters.
    """
    n_before = small_adata.n_obs
    adata = detect_doublets(small_adata.copy(), seed=42)

    assert adata.n_obs == n_before, (
        f"detect_doublets changed n_obs from {n_before} to {adata.n_obs}"
    )


def test_detect_doublets_does_not_change_n_vars(small_adata):
    """Test that detect_doublets does not alter the gene set.

    The gene dimension must be unchanged by doublet detection.
    """
    n_before = small_adata.n_vars
    adata = detect_doublets(small_adata.copy(), seed=42)

    assert adata.n_vars == n_before, (
        f"detect_doublets changed n_vars from {n_before} to {adata.n_vars}"
    )


def test_detect_doublets_with_explicit_threshold(small_adata):
    """Test that detect_doublets respects an explicit numeric threshold.

    When a float threshold (e.g. 0.3) is passed, every cell with a doublet
    score >= threshold must be marked as a doublet and vice versa.  This
    validates that the threshold parameter is actually applied.
    """
    threshold = 0.3
    adata = detect_doublets(small_adata.copy(), seed=42, threshold=threshold)

    scores = adata.obs["annqc_doublet_score"]
    is_doublet = adata.obs["annqc_is_doublet"]

    # For cells with valid scores: score >= threshold iff is_doublet is True
    valid_mask = scores.notna()
    for idx in adata.obs_names[valid_mask]:
        expected = scores[idx] >= threshold
        actual = is_doublet[idx]
        assert actual == expected, (
            f"Cell {idx}: score={scores[idx]:.4f}, threshold={threshold}, "
            f"is_doublet={actual} (expected {expected})"
        )


def test_detect_doublets_auto_threshold_produces_calls(small_adata):
    """Test that the 'auto' threshold still produces boolean calls for every cell.

    When threshold='auto' (the default), scrublet chooses a threshold
    automatically.  The result must still be a fully-populated bool column.
    """
    adata = detect_doublets(small_adata.copy(), seed=42, threshold="auto")

    is_doublet = adata.obs["annqc_is_doublet"]

    # No NaN in is_doublet — must be bool for every cell
    assert is_doublet.notna().all(), (
        "annqc_is_doublet has NaN values even when threshold='auto'"
    )


def test_detect_doublets_score_count_matches_n_obs(small_adata):
    """Test that annqc_doublet_score has exactly one entry per cell.

    The length of the score column must equal adata.n_obs so every cell
    is accounted for.
    """
    adata = detect_doublets(small_adata.copy(), seed=42)

    assert len(adata.obs["annqc_doublet_score"]) == adata.n_obs, (
        "Length of annqc_doublet_score does not match n_obs"
    )


def test_scrublet_failure_sets_incomplete_status(small_adata):
    """When Scrublet raises an exception, status must be INCOMPLETE."""
    import sys

    adata = small_adata.copy()
    adata.uns["annqc"] = {
        "status": "PASS",
        "warnings": [],
        "thresholds": {"doublet_threshold": None},
    }

    # Scrublet is imported lazily inside the try block; force ImportError via sys.modules
    with patch.dict(sys.modules, {"scrublet": None}):
        result = detect_doublets(adata, seed=42)

    assert result.uns["annqc"]["doublet_status"] == "FAILED"
    assert result.uns["annqc"]["status"] == "INCOMPLETE"
    assert result.uns["annqc"].get("doublet_failure_reason") is not None
    assert any("CRITICAL" in w for w in result.uns["annqc"]["warnings"])


def test_no_doublet_detection_flag(pbmc3k):
    """Running with no_doublet_detection=True must set doublet_status to SKIPPED."""
    import copy

    adata = pbmc3k.copy()
    result = annqc.run(
        adata,
        config=copy.deepcopy(DEFAULT_CONFIG),
        seed=42,
        input_file="pbmc3k",
        no_doublet_detection=True,
    )

    assert result.uns["annqc"]["doublet_status"] == "SKIPPED"


def test_raw_counts_layer_exists(pbmc3k):
    """After a full pipeline run, raw counts must be saved before normalization."""
    import copy

    import numpy as np

    adata = pbmc3k.copy()
    result = annqc.run(
        adata,
        config=copy.deepcopy(DEFAULT_CONFIG),
        seed=42,
        input_file="pbmc3k",
    )

    assert "counts" in result.layers, "adata.layers['counts'] was not created"
    # counts layer (raw integers) must differ from adata.X (log-normalized floats)
    x_vals = result.X if not hasattr(result.X, "toarray") else result.X.toarray()
    counts_vals = result.layers["counts"] if not hasattr(result.layers["counts"], "toarray") else result.layers["counts"].toarray()
    assert not np.allclose(x_vals, counts_vals), (
        "adata.layers['counts'] appears identical to adata.X — "
        "raw counts should differ from log-normalized values"
    )
