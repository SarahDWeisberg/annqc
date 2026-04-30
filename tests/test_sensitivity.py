"""Tests for annqc.sensitivity.run_sensitivity_analysis."""

import json
import os

import pytest

from annqc.sensitivity import run_sensitivity_analysis


def test_sensitivity_runs_on_small_adata(small_adata):
    result = run_sensitivity_analysis(small_adata, seed=0)
    assert isinstance(result, dict)


def test_sensitivity_runs_on_pbmc3k(pbmc3k):
    result = run_sensitivity_analysis(pbmc3k, seed=0)
    assert isinstance(result, dict)
    assert "n_input" in result
    assert result["n_input"] > 0


def test_sensitivity_output_contains_all_metrics(small_adata):
    result = run_sensitivity_analysis(small_adata, seed=0)
    for metric in ("mito", "min_genes", "max_genes", "min_counts"):
        assert metric in result, f"Missing metric: {metric}"


def test_sensitivity_each_metric_has_thresholds_tested(small_adata):
    result = run_sensitivity_analysis(small_adata, seed=0)
    for metric in ("mito", "min_genes", "max_genes", "min_counts"):
        assert "thresholds_tested" in result[metric]
        assert isinstance(result[metric]["thresholds_tested"], list)
        assert len(result[metric]["thresholds_tested"]) > 0


def test_sensitivity_cells_removed_matches_thresholds_length(small_adata):
    result = run_sensitivity_analysis(small_adata, seed=0)
    for metric in ("mito", "min_genes", "max_genes", "min_counts"):
        n_thr = len(result[metric]["thresholds_tested"])
        assert len(result[metric]["cells_removed"]) == n_thr


def test_sensitivity_cells_removed_are_integers(small_adata):
    result = run_sensitivity_analysis(small_adata, seed=0)
    for metric in ("mito", "min_genes", "max_genes", "min_counts"):
        for v in result[metric]["cells_removed"]:
            assert isinstance(v, int), f"{metric} cells_removed value {v!r} is not int"


def test_sensitivity_pct_removed_in_range(small_adata):
    result = run_sensitivity_analysis(small_adata, seed=0)
    for metric in ("mito", "min_genes", "max_genes", "min_counts"):
        for pct in result[metric]["pct_removed"]:
            assert 0.0 <= pct <= 100.0, f"{metric} pct_removed {pct} out of range"


def test_sensitivity_mad_suggested_present(small_adata):
    result = run_sensitivity_analysis(small_adata, seed=0)
    for metric in ("mito", "min_genes", "max_genes", "min_counts"):
        assert "mad_suggested" in result[metric]


def test_sensitivity_report_generated(small_adata, tmp_path):
    report_path = str(tmp_path / "sensitivity_test.html")
    run_sensitivity_analysis(small_adata, seed=0, report_path=report_path)
    assert os.path.exists(report_path)
    assert os.path.getsize(report_path) > 0


def test_sensitivity_json_output(small_adata, tmp_path):
    output_path = str(tmp_path / "sensitivity_test.json")
    run_sensitivity_analysis(small_adata, seed=0, output_path=output_path)
    assert os.path.exists(output_path)
    with open(output_path) as fh:
        data = json.load(fh)
    assert isinstance(data, dict)


def test_profiles_returns_three_profiles(small_adata):
    """When profiles=True, result must contain strict, standard, and permissive entries."""
    result = run_sensitivity_analysis(small_adata, seed=0, profiles=True)
    assert "profiles" in result, "profiles key missing from result when profiles=True"
    profile_names = {p["profile"] for p in result["profiles"] if "profile" in p}
    for expected in ("strict", "standard", "permissive"):
        assert expected in profile_names, f"Profile '{expected}' missing from profiles result"


def test_profiles_strict_retains_fewer_than_permissive(small_adata):
    """Strict profile should retain fewer or equal cells than permissive profile."""
    result = run_sensitivity_analysis(small_adata, seed=0, profiles=True)
    profiles = {p["profile"]: p for p in result["profiles"] if "profile" in p and "error" not in p}
    if "strict" in profiles and "permissive" in profiles:
        assert profiles["strict"]["n_output"] <= profiles["permissive"]["n_output"], (
            f"strict n_output ({profiles['strict']['n_output']}) > "
            f"permissive n_output ({profiles['permissive']['n_output']})"
        )


def test_cluster_labels_used_when_provided(small_adata, tmp_path):
    """When a cluster labels CSV is provided, it should be used without error."""
    import pandas as pd

    labels_path = str(tmp_path / "labels.csv")
    labels = pd.Series(
        ["A"] * (small_adata.n_obs // 2) + ["B"] * (small_adata.n_obs - small_adata.n_obs // 2),
        index=small_adata.obs_names,
        name="cluster",
    )
    labels.to_csv(labels_path, header=True)

    result = run_sensitivity_analysis(small_adata, seed=0, cluster_labels_path=labels_path)
    assert isinstance(result, dict)
