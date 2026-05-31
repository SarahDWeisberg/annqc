"""Compare an incoming dataset's QC metrics against reference atlas profiles."""

import logging
import numpy as np

from annqc.reference.schema import (
    METRICS, load_profile, normalize_assay, confidence_level,
)

logger = logging.getLogger(__name__)

# Percentile thresholds for flagging
_FLAG_HIGH_PCT = 85   # dataset median above this percentile → HIGH
_FLAG_LOW_PCT  = 15   # dataset median below this percentile → LOW


def find_profile(tissue: str, assay: str, suspension_type: str = "cell") -> dict | None:
    """Find the best matching bundled reference profile.

    Tries exact assay match first, then falls back to broader aliases.
    Returns None when no profile is available.
    """
    assay_key = normalize_assay(assay)
    profile = load_profile(tissue, assay_key, suspension_type)
    if profile is None:
        logger.debug(
            "No reference profile for tissue=%r assay=%r suspension=%r",
            tissue, assay_key, suspension_type,
        )
    return profile


def _percentile_rank(value: float, reference_values: list[float]) -> float | None:
    """Return the percentile rank (0-100) of value within reference_values.

    Returns None if reference_values has fewer than 3 entries.
    """
    if len(reference_values) < 3:
        return None
    arr = np.array(reference_values, dtype=float)
    return float(100.0 * (arr < value).sum() / len(arr))


def compare_to_reference(
    adata,
    tissue: str,
    assay: str,
    suspension_type: str = "cell",
    profile: dict | None = None,
) -> dict | None:
    """Compare adata's dataset-level QC medians to a reference atlas profile.

    Parameters
    ----------
    adata : AnnData
        Must have QC metrics in obs (pct_counts_mt, n_genes_by_counts, etc.)
        already computed by calculate_qc_metrics.
    tissue : str
        Tissue name matching CELLxGENE tissue_general vocabulary.
    assay : str
        Assay name (raw CELLxGENE string or canonical short key).
    suspension_type : str
        "cell" or "nucleus".
    profile : dict or None
        Pre-loaded profile dict. If None, looked up from bundled profiles.

    Returns
    -------
    dict or None
        Comparison report, or None if no matching profile or insufficient data.
        Keys: profile_used, tissue, assay, suspension_type, n_references,
              confidence, comparison (per-metric results).
    """
    if profile is None:
        profile = find_profile(tissue, assay, suspension_type)
    if profile is None:
        return None

    n_refs = profile.get("n_datasets", 0)
    conf = profile.get("confidence", confidence_level(n_refs))

    if conf == "insufficient":
        logger.info(
            "Reference profile for %s/%s/%s has insufficient data (n=%d) — skipping comparison",
            tissue, assay, suspension_type, n_refs,
        )
        return {
            "profile_used": f"{tissue}__{normalize_assay(assay)}__{suspension_type}",
            "tissue": tissue,
            "assay": normalize_assay(assay),
            "suspension_type": suspension_type,
            "n_references": n_refs,
            "confidence": "insufficient",
            "comparison": {},
            "note": f"Insufficient reference data (n={n_refs}, need ≥5).",
        }

    comparison = {}
    metrics_in_profile = profile.get("metrics", {})

    for metric in METRICS:
        if metric not in metrics_in_profile:
            continue
        if metric not in adata.obs.columns:
            continue

        metric_profile = metrics_in_profile[metric]
        reference_medians = metric_profile.get("dataset_medians", [])
        summary = metric_profile.get("summary", {})

        if len(reference_medians) < 3:
            continue

        # Dataset-level statistic: median of per-cell values
        dataset_median = float(np.median(adata.obs[metric].dropna().values))

        pct_rank = _percentile_rank(dataset_median, reference_medians)

        if pct_rank is None:
            flag = "UNKNOWN"
        elif pct_rank >= _FLAG_HIGH_PCT:
            flag = "HIGH"
        elif pct_rank <= _FLAG_LOW_PCT:
            flag = "LOW"
        else:
            flag = "NORMAL"

        comparison[metric] = {
            "dataset_median": round(dataset_median, 3),
            "reference_median": summary.get("median"),
            "reference_q25": summary.get("q25"),
            "reference_q75": summary.get("q75"),
            "reference_p5": summary.get("p5"),
            "reference_p95": summary.get("p95"),
            "percentile": round(pct_rank, 1) if pct_rank is not None else None,
            "flag": flag,
            "n_references": len(reference_medians),
        }

    if not comparison:
        return None

    return {
        "profile_used": profile.get("_filename",
            f"{tissue}__{normalize_assay(assay)}__{suspension_type}"),
        "tissue": tissue,
        "assay": normalize_assay(assay),
        "suspension_type": suspension_type,
        "n_references": n_refs,
        "confidence": conf,
        "census_version": profile.get("census_version"),
        "generated_date": profile.get("generated_date"),
        "comparison": comparison,
    }


def reference_warnings(comparison_result: dict) -> list[str]:
    """Generate human-readable warnings from a comparison result dict."""
    if not comparison_result or not comparison_result.get("comparison"):
        return []

    warnings = []
    conf = comparison_result.get("confidence", "high")
    n_refs = comparison_result.get("n_references", 0)
    tissue = comparison_result.get("tissue", "")
    assay = comparison_result.get("assay", "")

    conf_note = ""
    if conf == "medium":
        conf_note = f" (moderate confidence, n={n_refs} references)"
    elif conf == "low":
        conf_note = f" (low confidence, n={n_refs} references — treat with caution)"

    metric_labels = {
        "pct_counts_mt": "Mitochondrial %",
        "n_genes_by_counts": "Genes per cell",
        "total_counts": "UMI counts",
        "pct_counts_ribo": "Ribosomal %",
    }

    for metric, result in comparison_result["comparison"].items():
        label = metric_labels.get(metric, metric)
        flag = result.get("flag")
        pct = result.get("percentile")
        dataset_val = result.get("dataset_median")
        ref_median = result.get("reference_median")

        if flag == "HIGH":
            warnings.append(
                f"Reference: {label} dataset median ({dataset_val:.3g}) is at the "
                f"{pct:.0f}th percentile of {n_refs} {tissue}/{assay} reference datasets "
                f"(reference median: {ref_median:.3g}){conf_note}."
            )
        elif flag == "LOW":
            warnings.append(
                f"Reference: {label} dataset median ({dataset_val:.3g}) is at the "
                f"{pct:.0f}th percentile of {n_refs} {tissue}/{assay} reference datasets "
                f"(reference median: {ref_median:.3g}){conf_note}."
            )

    return warnings
