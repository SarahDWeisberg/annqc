"""Defines the adata.uns['annqc'] schema and validation helpers."""

import logging

logger = logging.getLogger(__name__)

REQUIRED_KEYS = [
    "version", "date", "seed", "config", "input_file",
    "thresholds", "cell_counts", "warnings", "status",
    "dry_run", "threshold_method", "raw_obs_metrics",
    "threshold_explanations",
]
REQUIRED_THRESHOLD_KEYS = [
    "mito_max_pct", "min_genes", "max_genes", "min_counts", "doublet_threshold",
    "method",
]
REQUIRED_CELL_COUNT_KEYS = [
    "input", "after_mito_filter", "after_gene_filter",
    "after_count_filter", "after_doublet_filter", "output",
]


def init_record(version: str, date: str, seed: int, config: dict, input_file: str) -> dict:
    """Return an initialized uns['annqc'] record with all required keys."""
    return {
        "version": version,
        "date": date,
        "seed": seed,
        "config": config,
        "input_file": input_file,
        "dry_run": False,
        "threshold_method": "manual",
        "raw_obs_metrics": {},
        "thresholds": {
            "method": "manual",
            "mito_max_pct": None,
            "min_genes": None,
            "max_genes": None,
            "min_counts": None,
            "doublet_threshold": None,
        },
        "cell_counts": {
            "input": 0,
            "after_mito_filter": 0,
            "after_gene_filter": 0,
            "after_count_filter": 0,
            "after_doublet_filter": 0,
            "output": 0,
        },
        "per_sample": {},
        "warnings": [],
        "status": "PASS",
        "doublet_status": None,
        "doublet_failure_reason": None,
        "threshold_explanations": {},
    }


def validate_record(record: dict) -> None:
    """Raise ValueError if record is missing required keys."""
    for key in REQUIRED_KEYS:
        if key not in record:
            raise ValueError(f"uns['annqc'] is missing required key: '{key}'")
    for key in REQUIRED_THRESHOLD_KEYS:
        if key not in record.get("thresholds", {}):
            raise ValueError(f"uns['annqc']['thresholds'] is missing required key: '{key}'")
    for key in REQUIRED_CELL_COUNT_KEYS:
        if key not in record.get("cell_counts", {}):
            raise ValueError(f"uns['annqc']['cell_counts'] is missing required key: '{key}'")
    logger.debug("uns['annqc'] record validation passed")
