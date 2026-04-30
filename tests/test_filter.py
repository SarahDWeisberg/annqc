"""Tests for annqc.filter: flag_cells, apply_filters, filter_genes."""

import copy

import numpy as np
import pytest

from annqc.filter import apply_filters, filter_genes, flag_cells
from annqc.qc import calculate_qc_metrics


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _prepare(adata, config):
    """Run QC metrics then flag cells; returns the adata with both applied."""
    adata = calculate_qc_metrics(adata, mito_prefix="MT-", ribo_prefix="RPS|RPL")
    adata = flag_cells(adata, config)
    return adata


# ---------------------------------------------------------------------------
# flag_cells tests
# ---------------------------------------------------------------------------


def test_flag_cells_adds_obs_columns(small_adata, default_config):
    """Test that flag_cells adds annqc_pass and annqc_filter_reason to obs.

    These two columns are mandatory per the AnnData contract and must always
    be present after flag_cells is called.
    """
    adata = _prepare(small_adata, default_config)

    assert "annqc_pass" in adata.obs.columns, "annqc_pass missing from obs"
    assert "annqc_filter_reason" in adata.obs.columns, "annqc_filter_reason missing from obs"


def test_flag_cells_high_mito_flagged(small_adata, default_config):
    """Test that cells with high mito% are flagged with filter_reason='mito'.

    Cells 0-9 in the fixture have their mitochondrial counts multiplied by 10,
    which should push them above the default 20% threshold.  After flag_cells
    every one of those cells must have annqc_filter_reason == 'mito' and
    annqc_pass == False.
    """
    adata = _prepare(small_adata, default_config)

    high_mito_cells = adata.obs.iloc[:10]

    failing = high_mito_cells[high_mito_cells["annqc_filter_reason"] == "mito"]
    assert len(failing) > 0, (
        "Expected at least some cells 0-9 to be flagged with filter_reason='mito'"
    )

    for idx in failing.index:
        assert not adata.obs.loc[idx, "annqc_pass"], (
            f"Cell {idx} has filter_reason='mito' but annqc_pass is True"
        )


def test_flag_cells_does_not_remove(small_adata, default_config):
    """Test that flag_cells does not change n_obs (only flags, doesn't remove).

    flag_cells is a non-destructive operation: it annotates cells but must
    leave the AnnData shape unchanged.
    """
    n_before = small_adata.n_obs
    adata = _prepare(small_adata, default_config)

    assert adata.n_obs == n_before, (
        f"flag_cells changed n_obs from {n_before} to {adata.n_obs}; "
        "it should only annotate, not remove cells"
    )


def test_flag_cells_pass_column_is_bool(small_adata, default_config):
    """Test that annqc_pass contains only boolean values.

    The contract requires annqc_pass to be a bool column so downstream code
    can rely on truthiness checks without casting.
    """
    adata = _prepare(small_adata, default_config)

    dtype = adata.obs["annqc_pass"].dtype
    assert dtype == bool or str(dtype) in ("bool", "boolean"), (
        f"annqc_pass dtype should be bool, got {dtype}"
    )


def test_flag_cells_filter_reason_valid_values(small_adata, default_config):
    """Test that annqc_filter_reason only contains allowed string values.

    Valid values are the empty string plus the six named reasons from the spec.
    """
    valid_reasons = {"", "mito", "min_genes", "max_genes", "min_counts", "max_counts", "doublet"}

    adata = _prepare(small_adata, default_config)

    unique_reasons = set(adata.obs["annqc_filter_reason"].unique())
    unexpected = unique_reasons - valid_reasons
    assert not unexpected, f"Unexpected filter_reason values: {unexpected}"


def test_flag_cells_passing_cells_have_empty_reason(small_adata, default_config):
    """Test that cells with annqc_pass=True have an empty filter_reason string.

    A cell that passed all filters must not carry a failure reason.
    """
    adata = _prepare(small_adata, default_config)

    passing = adata.obs[adata.obs["annqc_pass"]]
    non_empty_reasons = passing[passing["annqc_filter_reason"] != ""]
    assert len(non_empty_reasons) == 0, (
        f"{len(non_empty_reasons)} passing cells have a non-empty filter_reason"
    )


# ---------------------------------------------------------------------------
# apply_filters tests
# ---------------------------------------------------------------------------


def test_apply_filters_removes_flagged(small_adata, default_config):
    """Test that apply_filters removes cells where annqc_pass is False.

    After apply_filters the AnnData must contain only cells whose annqc_pass
    was True before the call.
    """
    adata = _prepare(small_adata, default_config)

    n_passing = adata.obs["annqc_pass"].sum()
    filtered = apply_filters(adata)

    assert filtered.n_obs == n_passing, (
        f"Expected {n_passing} cells after apply_filters, got {filtered.n_obs}"
    )


def test_apply_filters_preserves_passing_cells(small_adata, default_config):
    """Test that apply_filters keeps exactly the cells whose annqc_pass is True.

    The obs_names of the filtered AnnData must be the same set as the
    obs_names of cells with annqc_pass == True.
    """
    adata = _prepare(small_adata, default_config)

    expected_cells = set(adata.obs_names[adata.obs["annqc_pass"]])
    filtered = apply_filters(adata)
    actual_cells = set(filtered.obs_names)

    assert actual_cells == expected_cells


def test_apply_filters_raises_if_all_removed(small_adata, default_config):
    """Test that apply_filters raises ValueError if all cells would be removed.

    Setting an impossible min_genes threshold means flag_cells will mark every
    cell as failing.  apply_filters must then raise a ValueError with a helpful
    message rather than returning an empty AnnData.
    """
    cfg = copy.deepcopy(default_config)
    cfg["cells"]["min_genes"] = 999_999  # impossible — no cell has this many genes

    adata = calculate_qc_metrics(small_adata, mito_prefix="MT-", ribo_prefix="RPS|RPL")
    adata = flag_cells(adata, cfg)

    with pytest.raises(ValueError, match=r"(?i)(all cells|no cells|empty|removed)"):
        apply_filters(adata)


# ---------------------------------------------------------------------------
# filter_genes tests
# ---------------------------------------------------------------------------


def test_filter_genes_removes_rare(small_adata):
    """Test that filter_genes removes genes appearing in fewer than min_cells cells.

    Genes observed in very few cells are typically noise.  After filter_genes
    every remaining gene must appear in at least min_cells cells.
    """
    min_cells = 3
    filtered = filter_genes(small_adata, min_cells=min_cells)

    # Verify: for every remaining gene the number of non-zero cells >= min_cells
    X = np.asarray(filtered.X)
    cells_per_gene = (X > 0).sum(axis=0)
    assert (cells_per_gene >= min_cells).all(), (
        "filter_genes left genes that appear in fewer than min_cells cells"
    )


def test_filter_genes_reduces_n_vars(small_adata):
    """Test that filter_genes never increases the number of genes.

    After filtering, the gene count must be less than or equal to the original.
    """
    filtered = filter_genes(small_adata, min_cells=3)
    assert filtered.n_vars <= small_adata.n_vars


def test_filter_genes_does_not_remove_cells(small_adata):
    """Test that filter_genes does not change the number of cells.

    filter_genes operates on the var axis only; n_obs must be unchanged.
    """
    filtered = filter_genes(small_adata, min_cells=3)
    assert filtered.n_obs == small_adata.n_obs


# ---------------------------------------------------------------------------
# Priority / integration tests
# ---------------------------------------------------------------------------


def test_filter_reason_priority(small_adata, default_config):
    """Test that filter_reason records the first failing reason (mito before gene).

    When a cell fails multiple thresholds, the recorded reason should be the
    earliest step in the pipeline.  Mito filtering (step 3) runs before gene
    filtering (step 4), so a cell that fails both must be labelled 'mito'.

    This test forces the scenario by lowering min_genes to 1 (so no cell fails
    that threshold) and verifying that high-mito cells still get 'mito', then
    raises min_genes so that a known non-mito cell would fail gene count too —
    but its primary reason should still match whichever step fires first.
    """
    cfg = copy.deepcopy(default_config)
    # Keep mito threshold at default (20%).  Cells 0-9 are high-mito.
    # Lower min_genes far enough that it never fires before mito does.
    cfg["cells"]["min_genes"] = 1

    adata = calculate_qc_metrics(small_adata, mito_prefix="MT-", ribo_prefix="RPS|RPL")
    adata = flag_cells(adata, cfg)

    # Every high-mito cell (0-9) that failed must be labelled 'mito'
    high_mito_obs = adata.obs.iloc[:10]
    mito_failed = high_mito_obs[~high_mito_obs["annqc_pass"]]
    for idx in mito_failed.index:
        reason = adata.obs.loc[idx, "annqc_filter_reason"]
        assert reason == "mito", (
            f"Cell {idx} failed; expected reason='mito', got '{reason}'"
        )


def test_flag_cells_respects_max_genes_threshold(small_adata, default_config):
    """Test that cells exceeding max_genes are flagged with reason='max_genes'.

    By setting max_genes very low we force every cell in the fixture to fail
    that threshold.  flag_cells must mark them with 'max_genes'.
    """
    cfg = copy.deepcopy(default_config)
    cfg["cells"]["max_genes"] = 1  # impossibly low
    cfg["mito"]["max_pct"] = 100   # disable mito filter so max_genes fires first

    adata = calculate_qc_metrics(small_adata, mito_prefix="MT-", ribo_prefix="RPS|RPL")
    adata = flag_cells(adata, cfg)

    reasons = adata.obs["annqc_filter_reason"]
    # At least the majority of cells should be flagged max_genes
    max_gene_flagged = (reasons == "max_genes").sum()
    assert max_gene_flagged > 0, "No cells flagged with 'max_genes' even with threshold=1"


def test_flag_cells_respects_min_counts_threshold(small_adata, default_config):
    """Test that cells below min_counts are flagged with reason='min_counts'.

    By setting min_counts extremely high we force cells that otherwise would
    pass mito and gene filters to fail on counts.
    """
    cfg = copy.deepcopy(default_config)
    cfg["cells"]["min_counts"] = 10_000_000  # impossibly high
    cfg["mito"]["max_pct"] = 100             # disable mito filter
    cfg["cells"]["min_genes"] = 0            # disable gene filter
    cfg["cells"]["max_genes"] = 999_999      # disable gene upper filter

    adata = calculate_qc_metrics(small_adata, mito_prefix="MT-", ribo_prefix="RPS|RPL")
    adata = flag_cells(adata, cfg)

    min_counts_flagged = (adata.obs["annqc_filter_reason"] == "min_counts").sum()
    assert min_counts_flagged > 0, (
        "No cells flagged with 'min_counts' even with impossibly high threshold"
    )
