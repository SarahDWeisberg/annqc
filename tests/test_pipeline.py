"""Integration and snapshot tests for the full annqc pipeline."""

import copy
import json
import os

import pytest

import annqc
from annqc.config import DEFAULT_CONFIG


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

REQUIRED_OBS_COLS = [
    "annqc_pass",
    "annqc_doublet_score",
    "annqc_is_doublet",
    "annqc_filter_reason",
]

REQUIRED_UNS_KEYS = [
    "version",
    "date",
    "seed",
    "config",
    "input_file",
    "thresholds",
    "cell_counts",
    "warnings",
    "status",
]

REQUIRED_THRESHOLD_KEYS = [
    "mito_max_pct",
    "min_genes",
    "max_genes",
    "min_counts",
    "doublet_threshold",
]

REQUIRED_CELL_COUNT_KEYS = [
    "input",
    "after_mito_filter",
    "after_gene_filter",
    "after_count_filter",
    "after_doublet_filter",
    "output",
]

_SNAPSHOT_PATH = os.path.join(
    os.path.dirname(__file__), "snapshots", "pbmc3k_summary.json"
)


# ---------------------------------------------------------------------------
# End-to-end pipeline tests (use pbmc3k session fixture)
# ---------------------------------------------------------------------------


def test_full_pipeline_runs(pbmc3k):
    """Test that the pipeline runs end-to-end on PBMC 3k without error.

    A successful run must return a non-None AnnData object with at least
    one cell.
    """
    adata = pbmc3k.copy()
    result = annqc.run(adata, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    assert result is not None
    assert result.n_obs > 0


def test_pipeline_adds_required_obs_columns(pbmc3k):
    """Test that all required obs columns are present after pipeline run.

    The per-cell contract defines four columns that must always be written
    regardless of which filters fired.
    """
    adata = pbmc3k.copy()
    result = annqc.run(adata, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    for col in REQUIRED_OBS_COLS:
        assert col in result.obs.columns, f"Missing required obs column: {col}"


def test_pipeline_uns_schema(pbmc3k):
    """Test that adata.uns['annqc'] contains all required fields.

    Every run must produce a complete reproducibility record.  A missing key
    means downstream tooling (report builder, snapshot tests) will break.
    """
    adata = pbmc3k.copy()
    result = annqc.run(adata, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    assert "annqc" in result.uns, "adata.uns['annqc'] was not written by the pipeline"

    record = result.uns["annqc"]
    for key in REQUIRED_UNS_KEYS:
        assert key in record, f"Missing key in adata.uns['annqc']: '{key}'"


def test_pipeline_thresholds_schema(pbmc3k):
    """Test that adata.uns['annqc']['thresholds'] contains all required threshold keys.

    The thresholds sub-dict must document every value that was actually applied
    so the run can be reproduced exactly.
    """
    adata = pbmc3k.copy()
    result = annqc.run(adata, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    thresholds = result.uns["annqc"]["thresholds"]
    for key in REQUIRED_THRESHOLD_KEYS:
        assert key in thresholds, f"Missing key in thresholds: '{key}'"


def test_pipeline_cell_counts_schema(pbmc3k):
    """Test that adata.uns['annqc']['cell_counts'] contains all pipeline step counts.

    Each filtering step must record its output count so the filtering decisions
    table in the HTML report can be built.
    """
    adata = pbmc3k.copy()
    result = annqc.run(adata, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    cell_counts = result.uns["annqc"]["cell_counts"]
    for key in REQUIRED_CELL_COUNT_KEYS:
        assert key in cell_counts, f"Missing key in cell_counts: '{key}'"


def test_pipeline_status_is_valid(pbmc3k):
    """Test that pipeline status is PASS, FAIL, or INCOMPLETE.

    The status field must be one of exactly three values.  Any other string
    indicates an implementation error.
    """
    adata = pbmc3k.copy()
    result = annqc.run(adata, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    assert result.uns["annqc"]["status"] in ("PASS", "FAIL", "INCOMPLETE"), (
        f"Unexpected status: '{result.uns['annqc']['status']}'"
    )


def test_pipeline_version_is_string(pbmc3k):
    """Test that the version recorded in uns is a non-empty string.

    The version field enables users and support teams to reproduce results
    with the exact same code.
    """
    adata = pbmc3k.copy()
    result = annqc.run(adata, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    version = result.uns["annqc"]["version"]
    assert isinstance(version, str) and len(version) > 0, (
        f"version should be a non-empty string, got {version!r}"
    )


def test_pipeline_seed_recorded(pbmc3k):
    """Test that the seed passed to annqc.run is recorded in uns.

    The recorded seed must match the value passed so the run can be reproduced.
    """
    adata = pbmc3k.copy()
    result = annqc.run(adata, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    assert result.uns["annqc"]["seed"] == 42, (
        f"Seed in uns should be 42, got {result.uns['annqc']['seed']}"
    )


def test_pipeline_config_recorded(pbmc3k):
    """Test that the config dict is stored verbatim in uns['annqc']['config'].

    The stored config must be a dict and must contain the top-level keys of
    DEFAULT_CONFIG so the run is fully reproducible from the record alone.
    """
    adata = pbmc3k.copy()
    result = annqc.run(adata, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    stored_config = result.uns["annqc"]["config"]
    assert isinstance(stored_config, dict), "uns['annqc']['config'] should be a dict"

    for key in DEFAULT_CONFIG:
        assert key in stored_config, (
            f"uns['annqc']['config'] is missing top-level key '{key}'"
        )


def test_pipeline_input_file_recorded(pbmc3k):
    """Test that the input_file value is written to uns.

    The recorded input_file must exactly match what was passed to annqc.run.
    """
    adata = pbmc3k.copy()
    result = annqc.run(
        adata, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k_test"
    )

    assert result.uns["annqc"]["input_file"] == "pbmc3k_test", (
        f"input_file mismatch: {result.uns['annqc']['input_file']!r}"
    )


def test_pipeline_warnings_is_list(pbmc3k):
    """Test that uns['annqc']['warnings'] is always a list.

    Even when there are no warnings the field must be an empty list, not None
    or a missing key.
    """
    adata = pbmc3k.copy()
    result = annqc.run(adata, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    warnings = result.uns["annqc"]["warnings"]
    assert isinstance(warnings, list), (
        f"uns['annqc']['warnings'] should be a list, got {type(warnings)}"
    )


def test_pipeline_doublet_threshold_is_float(pbmc3k):
    """Test that the resolved doublet_threshold in thresholds is a float.

    Even when threshold='auto' is configured, the record must store the
    resolved numeric value that was actually used.
    """
    adata = pbmc3k.copy()
    result = annqc.run(adata, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    threshold = result.uns["annqc"]["thresholds"]["doublet_threshold"]
    assert isinstance(threshold, float), (
        f"doublet_threshold should be a float, got {type(threshold).__name__}: {threshold!r}"
    )


# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------


def test_full_pipeline_deterministic(pbmc3k):
    """Test that seeded pipeline produces identical cell counts on repeated runs.

    Running with seed=42 twice on identical input must yield the same n_obs.
    This is the primary reproducibility guarantee of the package.
    """
    adata1 = pbmc3k.copy()
    adata2 = pbmc3k.copy()

    result1 = annqc.run(adata1, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")
    result2 = annqc.run(adata2, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    assert result1.n_obs == result2.n_obs, (
        f"Seeded pipeline is not deterministic: "
        f"run1 n_obs={result1.n_obs}, run2 n_obs={result2.n_obs}"
    )


def test_full_pipeline_deterministic_cell_counts(pbmc3k):
    """Test that all per-step cell counts are identical across seeded runs.

    Beyond the final cell count, every intermediate count checkpoint must also
    match, confirming that every step respects the seed.
    """
    adata1 = pbmc3k.copy()
    adata2 = pbmc3k.copy()

    result1 = annqc.run(adata1, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")
    result2 = annqc.run(adata2, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    cc1 = result1.uns["annqc"]["cell_counts"]
    cc2 = result2.uns["annqc"]["cell_counts"]

    for key in REQUIRED_CELL_COUNT_KEYS:
        assert cc1[key] == cc2[key], (
            f"cell_counts['{key}'] differs: run1={cc1[key]}, run2={cc2[key]}"
        )


# ---------------------------------------------------------------------------
# Monotonicity
# ---------------------------------------------------------------------------


def test_pipeline_cell_counts_decreasing(pbmc3k):
    """Test that cell counts decrease monotonically through the pipeline steps.

    Each step can only remove cells, never add them.  A count that increases
    between steps indicates a bug in the pipeline orchestrator.
    """
    adata = pbmc3k.copy()
    result = annqc.run(adata, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    cc = result.uns["annqc"]["cell_counts"]

    assert cc["input"] >= cc["after_mito_filter"], (
        f"after_mito_filter ({cc['after_mito_filter']}) > input ({cc['input']})"
    )
    assert cc["after_mito_filter"] >= cc["after_gene_filter"], (
        f"after_gene_filter ({cc['after_gene_filter']}) > "
        f"after_mito_filter ({cc['after_mito_filter']})"
    )
    assert cc["after_gene_filter"] >= cc["after_count_filter"], (
        f"after_count_filter ({cc['after_count_filter']}) > "
        f"after_gene_filter ({cc['after_gene_filter']})"
    )
    assert cc["after_count_filter"] >= cc["after_doublet_filter"], (
        f"after_doublet_filter ({cc['after_doublet_filter']}) > "
        f"after_count_filter ({cc['after_count_filter']})"
    )
    assert cc["after_doublet_filter"] >= cc["output"], (
        f"output ({cc['output']}) > after_doublet_filter ({cc['after_doublet_filter']})"
    )


def test_pipeline_output_count_equals_n_obs(pbmc3k):
    """Test that cell_counts['output'] equals the actual n_obs of the result.

    The recorded output count must reflect reality — if they differ, the
    snapshot test and the filtering table in the report will be wrong.
    """
    adata = pbmc3k.copy()
    result = annqc.run(adata, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    assert result.uns["annqc"]["cell_counts"]["output"] == result.n_obs, (
        f"cell_counts['output']={result.uns['annqc']['cell_counts']['output']} "
        f"!= result.n_obs={result.n_obs}"
    )


def test_pipeline_input_count_equals_original_n_obs(pbmc3k):
    """Test that cell_counts['input'] equals the n_obs of the input AnnData.

    The recorded input count must match what was passed in so the filtering
    table percentages are correct.
    """
    adata = pbmc3k.copy()
    original_n_obs = adata.n_obs

    result = annqc.run(adata, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    assert result.uns["annqc"]["cell_counts"]["input"] == original_n_obs, (
        f"cell_counts['input']={result.uns['annqc']['cell_counts']['input']} "
        f"!= original n_obs={original_n_obs}"
    )


# ---------------------------------------------------------------------------
# obs column contract
# ---------------------------------------------------------------------------


def test_pipeline_annqc_pass_is_bool(pbmc3k):
    """Test that the annqc_pass column is boolean dtype after a full pipeline run.

    Downstream code uses boolean indexing on this column; a non-bool dtype
    causes silent correctness bugs.
    """
    adata = pbmc3k.copy()
    result = annqc.run(adata, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    dtype = result.obs["annqc_pass"].dtype
    assert dtype == bool or str(dtype) in ("bool", "boolean"), (
        f"annqc_pass dtype should be bool, got {dtype}"
    )


def test_pipeline_annqc_is_doublet_is_bool(pbmc3k):
    """Test that annqc_is_doublet is boolean dtype after a full pipeline run."""
    adata = pbmc3k.copy()
    result = annqc.run(adata, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    dtype = result.obs["annqc_is_doublet"].dtype
    assert dtype == bool or str(dtype) in ("bool", "boolean"), (
        f"annqc_is_doublet dtype should be bool, got {dtype}"
    )


def test_pipeline_annqc_doublet_score_range(pbmc3k):
    """Test that annqc_doublet_score values are in [0, 1] for all output cells.

    After the pipeline completes, every remaining cell's doublet score must be
    a valid probability in [0, 1].
    """
    adata = pbmc3k.copy()
    result = annqc.run(adata, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    scores = result.obs["annqc_doublet_score"]
    valid = scores.notna()

    assert (scores[valid] >= 0.0).all(), "Some output doublet scores are < 0"
    assert (scores[valid] <= 1.0).all(), "Some output doublet scores are > 1"


def test_pipeline_all_output_cells_pass(pbmc3k):
    """Test that every cell in the output AnnData has annqc_pass == True.

    The pipeline applies filters and removes failing cells, so the output
    must contain only cells that passed all QC steps.
    """
    adata = pbmc3k.copy()
    result = annqc.run(adata, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    n_failing = (~result.obs["annqc_pass"]).sum()
    assert n_failing == 0, (
        f"{n_failing} cells in the output have annqc_pass=False; "
        "the pipeline should have removed them"
    )


# ---------------------------------------------------------------------------
# Snapshot test
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# New spec-field tests (v0.2.0)
# ---------------------------------------------------------------------------


def test_record_has_dry_run_field(pbmc3k):
    """Pipeline run with default args must record dry_run=False in uns['annqc']."""
    adata = pbmc3k.copy()
    result = annqc.run(adata, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    record = result.uns["annqc"]
    assert "dry_run" in record, "uns['annqc'] is missing 'dry_run' field"
    assert record["dry_run"] is False, (
        f"Expected dry_run=False for a default run, got {record['dry_run']!r}"
    )


def test_record_has_threshold_method(pbmc3k):
    """Pipeline run must record threshold_method as 'manual' or 'auto_mad'."""
    adata = pbmc3k.copy()
    result = annqc.run(adata, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    record = result.uns["annqc"]
    assert "threshold_method" in record, "uns['annqc'] is missing 'threshold_method' field"
    assert record["threshold_method"] in ("manual", "auto_mad"), (
        f"threshold_method must be 'manual' or 'auto_mad', got {record['threshold_method']!r}"
    )


def test_record_has_raw_obs_metrics(pbmc3k):
    """Pipeline run must store raw_obs_metrics with pct_counts_mt present."""
    adata = pbmc3k.copy()
    result = annqc.run(adata, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    record = result.uns["annqc"]
    assert "raw_obs_metrics" in record, "uns['annqc'] is missing 'raw_obs_metrics' field"
    assert "pct_counts_mt" in record["raw_obs_metrics"], (
        "uns['annqc']['raw_obs_metrics'] is missing 'pct_counts_mt'"
    )


def test_thresholds_has_method_field(pbmc3k):
    """Pipeline run must include a 'method' key inside thresholds."""
    adata = pbmc3k.copy()
    result = annqc.run(adata, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    thresholds = result.uns["annqc"]["thresholds"]
    assert "method" in thresholds, "uns['annqc']['thresholds'] is missing 'method' field"


def test_per_sample_thresholds_come_from_config(pbmc3k):
    """Per-sample FAIL threshold must come from config, not a hardcoded value."""
    adata = pbmc3k.copy()
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    cfg["thresholds"]["min_cells_pass_per_sample"] = 9999  # impossibly high

    # Add a fake sample key so per-sample mode is triggered
    adata.obs["sample"] = "s1"

    result = annqc.run(
        adata,
        config=cfg,
        seed=42,
        input_file="pbmc3k",
        sample_key="sample",
    )

    per_sample = result.uns["annqc"].get("per_sample", {})
    assert "s1" in per_sample, "per_sample key 's1' not found"
    assert per_sample["s1"]["status"] == "FAIL", (
        f"Expected FAIL with min_cells_pass_per_sample=9999, "
        f"got {per_sample['s1']['status']!r}"
    )


def test_auto_thresholds_sets_method(pbmc3k):
    """Running with auto_thresholds=True must set threshold_method and thresholds['method'] to 'auto_mad'."""
    adata = pbmc3k.copy()
    result = annqc.run(
        adata,
        config=copy.deepcopy(DEFAULT_CONFIG),
        seed=42,
        input_file="pbmc3k",
        auto_thresholds=True,
    )

    record = result.uns["annqc"]
    assert record["threshold_method"] == "auto_mad", (
        f"Expected threshold_method='auto_mad', got {record['threshold_method']!r}"
    )
    assert record["thresholds"]["method"] == "auto_mad", (
        f"Expected thresholds['method']='auto_mad', got {record['thresholds']['method']!r}"
    )


# ---------------------------------------------------------------------------
# Snapshot test (unchanged)
# ---------------------------------------------------------------------------


def _nan_to_none(obj):
    """Recursively replace float NaN with None for JSON-safe comparison."""
    if isinstance(obj, dict):
        return {k: _nan_to_none(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_nan_to_none(v) for v in obj]
    if isinstance(obj, float) and obj != obj:  # NaN != NaN
        return None
    return obj


def test_summary_snapshot(pbmc3k):
    """Test that pipeline output matches stored snapshot (or create it).

    On the first run the snapshot file is written.  On subsequent runs the
    output is compared to the stored values.  A mismatch means either the
    pipeline changed deterministically (update the snapshot intentionally) or
    a regression was introduced.
    """
    adata = pbmc3k.copy()
    result = annqc.run(adata, config=copy.deepcopy(DEFAULT_CONFIG), seed=42, input_file="pbmc3k")

    summary = {
        "n_cells_output": result.n_obs,
        "thresholds": result.uns["annqc"]["thresholds"],
        "cell_counts": result.uns["annqc"]["cell_counts"],
        "status": result.uns["annqc"]["status"],
    }
    # Normalize NaN → None so the dict is JSON-round-trippable and equality works
    summary_normalized = _nan_to_none(summary)

    if os.path.exists(_SNAPSHOT_PATH):
        with open(_SNAPSHOT_PATH) as f:
            content = f.read().strip()

        if content and content != "{}":
            expected = json.loads(content)
            assert summary_normalized == expected, (
                "Pipeline output does not match stored snapshot.\n"
                f"Expected: {json.dumps(expected, indent=2)}\n"
                f"Got:      {json.dumps(summary_normalized, indent=2)}\n"
                "If this change is intentional, delete tests/snapshots/pbmc3k_summary.json "
                "and re-run to regenerate it."
            )
        else:
            # Placeholder file — write real snapshot
            os.makedirs(os.path.dirname(_SNAPSHOT_PATH), exist_ok=True)
            with open(_SNAPSHOT_PATH, "w") as f:
                json.dump(summary_normalized, f, indent=2)
    else:
        os.makedirs(os.path.dirname(_SNAPSHOT_PATH), exist_ok=True)
        with open(_SNAPSHOT_PATH, "w") as f:
            json.dump(summary_normalized, f, indent=2)
