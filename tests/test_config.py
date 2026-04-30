"""Tests for annqc.config: loading, merging, validation, serialisation."""

import os
import tempfile

import pytest
import yaml

from annqc.config import (
    DEFAULT_CONFIG,
    config_to_yaml,
    get_default_config,
    load_config,
    validate_config,
)


# ---------------------------------------------------------------------------
# DEFAULT_CONFIG structure
# ---------------------------------------------------------------------------


def test_default_config_has_required_keys():
    """Test that DEFAULT_CONFIG has all required top-level keys.

    The spec defines six required sections.  Missing any of them would cause
    downstream code to crash with an unhelpful KeyError.
    """
    required_keys = {"mito", "ribo", "cells", "genes", "doublets", "normalization", "report"}
    missing = required_keys - set(DEFAULT_CONFIG.keys())
    assert not missing, f"DEFAULT_CONFIG is missing keys: {missing}"


def test_default_config_mito_section():
    """Test that the mito section contains expected sub-keys with correct defaults.

    The default mito prefix must be 'MT-' and max_pct must be 20.
    """
    mito = DEFAULT_CONFIG["mito"]
    assert "prefix" in mito, "DEFAULT_CONFIG['mito'] missing 'prefix'"
    assert "max_pct" in mito, "DEFAULT_CONFIG['mito'] missing 'max_pct'"
    assert mito["prefix"] == "MT-"
    assert mito["max_pct"] == 20


def test_default_config_cells_section():
    """Test that the cells section contains min_genes, max_genes, min_counts, max_counts.

    These four keys drive the core cell-level QC thresholds.
    """
    cells = DEFAULT_CONFIG["cells"]
    for key in ("min_genes", "max_genes", "min_counts", "max_counts"):
        assert key in cells, f"DEFAULT_CONFIG['cells'] missing '{key}'"


def test_default_config_doublets_method_is_scrublet():
    """Test that the default doublet detection method is 'scrublet'.

    v0.1 supports only scrublet; the default must reflect this.
    """
    assert DEFAULT_CONFIG["doublets"]["method"] == "scrublet"


# ---------------------------------------------------------------------------
# get_default_config
# ---------------------------------------------------------------------------


def test_get_default_config_returns_copy():
    """Test that get_default_config returns an independent copy of DEFAULT_CONFIG.

    Mutating the returned dict must not affect DEFAULT_CONFIG or a second call.
    """
    cfg1 = get_default_config()
    cfg2 = get_default_config()

    cfg1["mito"]["max_pct"] = 99
    assert DEFAULT_CONFIG["mito"]["max_pct"] != 99, (
        "Mutating get_default_config() result affected DEFAULT_CONFIG"
    )
    assert cfg2["mito"]["max_pct"] != 99, (
        "Mutating cfg1 affected cfg2 — they share the same object"
    )


# ---------------------------------------------------------------------------
# load_config
# ---------------------------------------------------------------------------


def test_load_config_file_not_found():
    """Test that load_config raises FileNotFoundError for missing files.

    Passing a path that doesn't exist must raise FileNotFoundError, not a
    generic KeyError or AttributeError.
    """
    with pytest.raises(FileNotFoundError):
        load_config("/nonexistent/path/to/config.yaml")


def test_load_config_returns_dict(tmp_path):
    """Test that load_config returns a dict when given a valid YAML file.

    A minimal config file with a single key must still produce a dict.
    """
    config_file = tmp_path / "minimal.yaml"
    config_file.write_text("mito:\n  max_pct: 15\n")

    result = load_config(str(config_file))
    assert isinstance(result, dict), f"load_config should return dict, got {type(result)}"


def test_load_config_merges_with_defaults(tmp_path):
    """Test that a partial YAML config merges correctly with defaults.

    When the user provides only a subset of keys, unspecified keys must
    fall back to their default values.  Specified keys must override defaults.
    """
    config_file = tmp_path / "partial.yaml"
    config_file.write_text("mito:\n  max_pct: 15\n")

    cfg = load_config(str(config_file))

    # Overridden value
    assert cfg["mito"]["max_pct"] == 15, (
        "load_config did not apply the user-specified mito.max_pct"
    )
    # Default value preserved for untouched key
    assert cfg["cells"]["min_genes"] == DEFAULT_CONFIG["cells"]["min_genes"], (
        "load_config dropped the default cells.min_genes when not specified"
    )


def test_load_config_full_override(tmp_path):
    """Test that a complete YAML config fully overrides all defaults.

    Writing every section to the YAML must result in those values being used
    rather than the defaults.
    """
    full_config = {
        "mito": {"prefix": "mt-", "max_pct": 30},
        "ribo": {"prefix": "Rps|Rpl", "max_pct": 40},
        "cells": {"min_genes": 100, "max_genes": 5000, "min_counts": 300, "max_counts": None},
        "genes": {"min_cells": 5},
        "doublets": {"method": "scrublet", "threshold": 0.25, "simulate_doublet_ratio": 3.0},
        "normalization": {"method": "log1p", "target_sum": 5000},
        "report": {"title": "My Report", "author": "Tester"},
    }
    config_file = tmp_path / "full.yaml"
    config_file.write_text(yaml.dump(full_config))

    cfg = load_config(str(config_file))

    assert cfg["mito"]["prefix"] == "mt-"
    assert cfg["mito"]["max_pct"] == 30
    assert cfg["cells"]["min_genes"] == 100
    assert cfg["normalization"]["target_sum"] == 5000


# ---------------------------------------------------------------------------
# validate_config
# ---------------------------------------------------------------------------


def test_validate_config_passes_default():
    """Test that DEFAULT_CONFIG passes validation without error.

    The built-in defaults must themselves be valid according to the validator.
    """
    import copy

    cfg = copy.deepcopy(DEFAULT_CONFIG)
    # Should not raise
    validate_config(cfg)


def test_validate_config_fails_bad_mito_pct():
    """Test that validate_config raises ValueError for mito.max_pct > 100.

    A percentage above 100 is physically impossible and must be rejected.
    """
    import copy

    cfg = copy.deepcopy(DEFAULT_CONFIG)
    cfg["mito"]["max_pct"] = 150

    with pytest.raises(ValueError, match=r"(?i)(mito|max_pct|percent|100)"):
        validate_config(cfg)


def test_validate_config_fails_negative_mito_pct():
    """Test that validate_config raises ValueError for mito.max_pct < 0.

    Negative percentages are meaningless and must be caught by validation.
    """
    import copy

    cfg = copy.deepcopy(DEFAULT_CONFIG)
    cfg["mito"]["max_pct"] = -5

    with pytest.raises(ValueError):
        validate_config(cfg)


def test_validate_config_fails_bad_method():
    """Test that validate_config raises ValueError for unknown normalization method.

    Only 'log1p' and 'none' are valid normalization methods in v0.1.
    Any other string must raise a ValueError.
    """
    import copy

    cfg = copy.deepcopy(DEFAULT_CONFIG)
    cfg["normalization"]["method"] = "scran"  # not supported in v0.1

    with pytest.raises(ValueError, match=r"(?i)(method|normalization|scran|unsupported|invalid)"):
        validate_config(cfg)


def test_validate_config_fails_min_genes_negative():
    """Test that validate_config raises ValueError when min_genes < 0.

    Negative cell thresholds are nonsensical and must be rejected.
    """
    import copy

    cfg = copy.deepcopy(DEFAULT_CONFIG)
    cfg["cells"]["min_genes"] = -10

    with pytest.raises(ValueError):
        validate_config(cfg)


def test_validate_config_fails_min_greater_than_max():
    """Test that validate_config raises ValueError when min_genes > max_genes.

    A window where the lower bound exceeds the upper bound can never contain
    any cell and must be rejected during validation.
    """
    import copy

    cfg = copy.deepcopy(DEFAULT_CONFIG)
    cfg["cells"]["min_genes"] = 8000
    cfg["cells"]["max_genes"] = 200

    with pytest.raises(ValueError):
        validate_config(cfg)


def test_validate_config_fails_unknown_doublet_method():
    """Test that validate_config raises ValueError for unsupported doublet methods.

    v0.1 supports only 'scrublet'.  Requesting another method must raise.
    """
    import copy

    cfg = copy.deepcopy(DEFAULT_CONFIG)
    cfg["doublets"]["method"] = "scDblFinder"  # R-based, not supported

    with pytest.raises(ValueError, match=r"(?i)(method|doublet|scrublet|unsupported|invalid)"):
        validate_config(cfg)


# ---------------------------------------------------------------------------
# config_to_yaml / round-trip
# ---------------------------------------------------------------------------


def test_config_to_yaml_returns_string():
    """Test that config_to_yaml returns a non-empty string.

    The serialised form must be a plain string that could be written to a file.
    """
    import copy

    cfg = copy.deepcopy(DEFAULT_CONFIG)
    result = config_to_yaml(cfg)

    assert isinstance(result, str), f"config_to_yaml should return str, got {type(result)}"
    assert len(result.strip()) > 0, "config_to_yaml returned an empty string"


def test_config_to_yaml_roundtrip():
    """Test that config serializes to YAML and back without data loss.

    Serialising DEFAULT_CONFIG to YAML and parsing the result must produce a
    dict equal to the original.  This ensures no data is dropped or corrupted.
    """
    import copy

    original = copy.deepcopy(DEFAULT_CONFIG)
    yaml_str = config_to_yaml(original)

    recovered = yaml.safe_load(yaml_str)

    # Top-level keys must be identical
    assert set(recovered.keys()) == set(original.keys()), (
        f"Round-trip lost keys. Original: {set(original)}, recovered: {set(recovered)}"
    )

    # Check nested values
    assert recovered["mito"]["max_pct"] == original["mito"]["max_pct"]
    assert recovered["cells"]["min_genes"] == original["cells"]["min_genes"]
    assert recovered["doublets"]["method"] == original["doublets"]["method"]
    assert recovered["normalization"]["method"] == original["normalization"]["method"]


def test_config_to_yaml_is_valid_yaml():
    """Test that config_to_yaml produces syntactically valid YAML.

    yaml.safe_load must be able to parse the output without raising any
    scanner or parser exception.
    """
    import copy

    cfg = copy.deepcopy(DEFAULT_CONFIG)
    yaml_str = config_to_yaml(cfg)

    # This must not raise
    parsed = yaml.safe_load(yaml_str)
    assert parsed is not None
