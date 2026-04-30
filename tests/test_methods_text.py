"""Tests for annqc.methods_text — generate_methods_text, generate_methods_short."""

import math

import anndata as ad
import numpy as np
import pytest

from annqc.methods_text import generate_methods_text, generate_methods_short


def _make_adata(threshold_method="manual", n_input=2638, n_output=2480, doublet_thr=0.25):
    """Return a minimal AnnData with adata.uns['annqc'] set."""
    adata = ad.AnnData(X=np.zeros((10, 10)))
    adata.uns["annqc"] = {
        "version": "0.3.0",
        "threshold_method": threshold_method,
        "config": {
            "mito": {"prefix": "MT-"},
            "ribo": {"prefix": "RPS|RPL"},
            "genes": {"min_cells": 3},
        },
        "thresholds": {
            "mito_max_pct": 20.0,
            "min_genes": 200.0,
            "max_genes": 6000.0,
            "min_counts": 500.0,
            "doublet_threshold": doublet_thr,
        },
        "cell_counts": {
            "input": n_input,
            "output": n_output,
            "after_count_filter": n_input - 100,
            "after_doublet_filter": n_input - 158,
        },
    }
    return adata


# ── generate_methods_text ────────────────────────────────────────────────────

def test_generate_methods_text_contains_version():
    adata = _make_adata()
    text = generate_methods_text(adata)
    assert "0.3.0" in text


def test_generate_methods_text_contains_cell_counts():
    adata = _make_adata(n_input=2638, n_output=2480)
    text = generate_methods_text(adata)
    assert "2,638" in text
    assert "2,480" in text


def test_generate_methods_text_manual_phrase():
    adata = _make_adata(threshold_method="manual")
    text = generate_methods_text(adata)
    assert "manually specified thresholds" in text


def test_generate_methods_text_auto_mad_phrase():
    adata = _make_adata(threshold_method="auto_mad")
    text = generate_methods_text(adata)
    assert "MAD-based adaptive thresholds" in text


def test_generate_methods_text_mito_prefix():
    adata = _make_adata()
    text = generate_methods_text(adata)
    assert "MT-" in text


def test_generate_methods_text_doublet_sentence():
    adata = _make_adata(doublet_thr=0.25)
    text = generate_methods_text(adata)
    assert "Scrublet" in text
    assert "0.250" in text


def test_generate_methods_text_nan_doublet_skipped():
    adata = _make_adata(doublet_thr=float("nan"))
    text = generate_methods_text(adata)
    assert "not performed" in text or "skipped" in text


def test_generate_methods_text_retention_pct():
    adata = _make_adata(n_input=2000, n_output=1800)
    text = generate_methods_text(adata)
    assert "90.0%" in text


# ── generate_methods_short ───────────────────────────────────────────────────

def test_generate_methods_short_basic():
    adata = _make_adata(n_input=2638, n_output=2480)
    text = generate_methods_short(adata)
    assert "AnnQC" in text
    assert "2,638" in text
    assert "2,480" in text


def test_generate_methods_short_auto_mad():
    adata = _make_adata(threshold_method="auto_mad")
    text = generate_methods_short(adata)
    assert "MAD-based" in text


def test_generate_methods_short_manual():
    adata = _make_adata(threshold_method="manual")
    text = generate_methods_short(adata)
    assert "standard thresholds" in text
