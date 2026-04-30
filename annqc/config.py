"""YAML config loading, validation, and defaults for AnnQC."""

import copy
import logging
import os

import yaml

logger = logging.getLogger(__name__)

DEFAULT_CONFIG: dict = {
    "mito": {"prefix": "MT-", "max_pct": 20},
    "ribo": {"prefix": "RPS|RPL", "max_pct": 50},
    "cells": {
        "min_genes": 200,
        "max_genes": 6000,
        "min_counts": 500,
        "max_counts": None,
    },
    "genes": {"min_cells": 3},
    "doublets": {
        "method": "scrublet",
        "threshold": "auto",
        "simulate_doublet_ratio": 2.0,
    },
    "normalization": {"method": "log1p", "target_sum": 10000},
    # These are suggestions only. Appropriate values depend on tissue type,
    # sequencing depth, and experimental design.
    "thresholds": {
        "min_cells_pass": 100,
        "min_cells_warn": 500,
        "min_cells_pass_per_sample": 200,
        "min_cells_warn_per_sample": 500,
    },
    "report": {"title": "Single-Cell QC Report", "author": ""},
}


def _deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge override into base, returning a new dict."""
    result = copy.deepcopy(base)
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def load_config(path: str) -> dict:
    """Load YAML config from path and deep-merge with DEFAULT_CONFIG.

    Missing keys fall back to DEFAULT_CONFIG values.
    Raises FileNotFoundError if path does not exist.
    Raises ValueError with field name if a value has wrong type or range.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Config file not found: {path}")
    with open(path) as fh:
        user_cfg = yaml.safe_load(fh) or {}
    logger.debug("Loaded config from %s", path)
    merged = _deep_merge(DEFAULT_CONFIG, user_cfg)
    validate_config(merged)
    return merged


def validate_config(config: dict) -> None:
    """Validate config dict for required fields and value types.

    Raises ValueError with the exact missing or invalid field name.
    """
    required_sections = ["mito", "ribo", "cells", "genes", "doublets", "normalization", "thresholds", "report"]
    for section in required_sections:
        if section not in config:
            raise ValueError(f"Config is missing required section: '{section}'")

    mito_pct = config["mito"].get("max_pct")
    if mito_pct is not None and not (0 <= float(mito_pct) <= 100):
        raise ValueError("mito.max_pct must be between 0 and 100")

    ribo_pct = config["ribo"].get("max_pct")
    if ribo_pct is not None and not (0 <= float(ribo_pct) <= 100):
        raise ValueError("ribo.max_pct must be between 0 and 100")

    if config["cells"].get("min_genes", 0) < 0:
        raise ValueError("cells.min_genes must be >= 0")

    min_genes = config["cells"].get("min_genes")
    max_genes = config["cells"].get("max_genes")
    if min_genes is not None and max_genes is not None and min_genes > max_genes:
        raise ValueError(
            f"cells.min_genes ({min_genes}) must be <= cells.max_genes ({max_genes})"
        )

    norm_method = config["normalization"].get("method", "log1p")
    if norm_method not in ("log1p", "none"):
        raise ValueError(
            f"normalization.method must be 'log1p' or 'none', got '{norm_method}'"
        )

    doublet_method = config["doublets"].get("method", "scrublet")
    if doublet_method not in ("scrublet",):
        raise ValueError(
            f"doublets.method must be 'scrublet', got '{doublet_method}'"
        )

    threshold = config["doublets"].get("threshold", "auto")
    if threshold != "auto":
        try:
            val = float(threshold)
            if not (0.0 <= val <= 1.0):
                raise ValueError(
                    "doublets.threshold must be 'auto' or a float between 0 and 1"
                )
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "doublets.threshold must be 'auto' or a float between 0 and 1"
            ) from exc

    logger.debug("Config validation passed")


def config_to_yaml(config: dict) -> str:
    """Serialize config dict to a YAML string."""
    return yaml.dump(config, default_flow_style=False, sort_keys=False)


def get_default_config() -> dict:
    """Return a deep copy of DEFAULT_CONFIG."""
    return copy.deepcopy(DEFAULT_CONFIG)
