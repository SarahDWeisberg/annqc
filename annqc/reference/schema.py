"""Reference profile schema: load, validate, and index bundled profiles."""

import json
import logging
import os
from pathlib import Path

logger = logging.getLogger(__name__)

SCHEMA_VERSION = "1.0"
_PROFILES_DIR = Path(__file__).parent / "profiles"

CONFIDENCE_THRESHOLDS = {
    "high": 20,
    "medium": 10,
    "low": 3,
}

METRICS = ["pct_counts_mt", "n_genes_by_counts", "total_counts", "pct_counts_ribo"]

# Canonical assay name mapping from CELLxGENE raw strings
ASSAY_ALIASES = {
    "10x 3' v2": "10x_v2",
    "10x 3' v3": "10x_v3",
    "10x 3' v3.1": "10x_v3",
    "10x 5' v1": "10x_5p",
    "10x 5' v2": "10x_5p",
    "10x Multiome": "10x_multiome",
    "Parse Biosciences": "parse",
    "Seq-Well": "seq_well",
    "Drop-seq": "dropseq",
    "Smart-seq2": "smart_seq2",
    "Smart-seq v4": "smart_seq2",
}


def normalize_assay(raw_assay: str) -> str:
    """Return a canonical short assay key from a CELLxGENE assay string."""
    return ASSAY_ALIASES.get(raw_assay, raw_assay.lower().replace(" ", "_").replace("'", ""))


def confidence_level(n: int) -> str:
    """Return confidence level string based on number of reference datasets."""
    if n >= CONFIDENCE_THRESHOLDS["high"]:
        return "high"
    if n >= CONFIDENCE_THRESHOLDS["medium"]:
        return "medium"
    if n >= CONFIDENCE_THRESHOLDS["low"]:
        return "low"
    return "insufficient"


def profile_filename(tissue: str, assay: str, suspension_type: str) -> str:
    """Return the canonical filename for a reference profile."""
    return f"{tissue}__{assay}__{suspension_type}.json"


def load_profile(tissue: str, assay: str, suspension_type: str = "cell") -> dict | None:
    """Load a bundled reference profile, or None if not available."""
    fname = profile_filename(tissue, assay, suspension_type)
    path = _PROFILES_DIR / fname
    if not path.exists():
        return None
    with open(path) as fh:
        return json.load(fh)


def list_available_profiles() -> list[dict]:
    """Return a list of profile metadata dicts for all bundled profiles."""
    profiles = []
    for p in sorted(_PROFILES_DIR.glob("*.json")):
        try:
            with open(p) as fh:
                data = json.load(fh)
            profiles.append({
                "tissue": data.get("tissue"),
                "assay": data.get("assay"),
                "suspension_type": data.get("suspension_type"),
                "n_datasets": data.get("n_datasets"),
                "n_cells_total": data.get("n_cells_total"),
                "confidence": data.get("confidence"),
                "generated_date": data.get("generated_date"),
                "filename": p.name,
            })
        except Exception as exc:
            logger.warning("Could not read profile %s: %s", p.name, exc)
    return profiles


def validate_profile(profile: dict) -> None:
    """Raise ValueError if a profile dict is missing required fields."""
    required_top = [
        "schema_version", "tissue", "assay", "suspension_type", "species",
        "generated_date", "census_version", "n_datasets", "n_cells_total",
        "confidence", "metrics",
    ]
    for key in required_top:
        if key not in profile:
            raise ValueError(f"Profile missing required field: {key!r}")
    for metric in METRICS:
        if metric not in profile["metrics"]:
            continue
        m = profile["metrics"][metric]
        if "dataset_medians" not in m or "summary" not in m:
            raise ValueError(
                f"Profile metric {metric!r} missing 'dataset_medians' or 'summary'"
            )
    required_summary = ["median", "q25", "q75", "mad", "p5", "p95", "n"]
    for metric, m in profile["metrics"].items():
        s = m.get("summary", {})
        for key in required_summary:
            if key not in s:
                raise ValueError(
                    f"Profile metric {metric!r} summary missing field: {key!r}"
                )


def save_profile(profile: dict, output_dir: str | Path | None = None) -> Path:
    """Write a profile dict to the canonical filename. Returns path written."""
    validate_profile(profile)
    out_dir = Path(output_dir) if output_dir else _PROFILES_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    fname = profile_filename(
        profile["tissue"], profile["assay"], profile["suspension_type"]
    )
    out_path = out_dir / fname
    with open(out_path, "w") as fh:
        json.dump(profile, fh, indent=2)
    return out_path
