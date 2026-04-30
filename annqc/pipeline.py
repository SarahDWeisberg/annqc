"""Main AnnQC pipeline orchestrator."""

import copy
import logging
import os
from datetime import datetime, timezone

import numpy as np
import scanpy as sc

logger = logging.getLogger(__name__)


def run(
    adata_or_path,
    config=None,
    sample_key=None,
    seed: int = 0,
    output=None,
    report_path=None,
    input_file=None,
    dry_run: bool = False,
    auto_thresholds: bool = False,
    no_doublet_detection: bool = False,
):
    """Run the full AnnQC pipeline.

    Parameters
    ----------
    adata_or_path : AnnData or str
        Input data: an AnnData, a path to .h5ad, a 10x directory, or a .h5 file.
    config : dict, str, or None
        Config dict, path to YAML, or None to use built-in defaults.
    sample_key : str or None
        obs column for per-sample statistics.
    seed : int
        Random seed for reproducibility.
    output : str or None
        Path to save cleaned .h5ad (directory created if needed). Skipped if None.
    report_path : str or None
        Path to save HTML report. Skipped if None.
    input_file : str or None
        Label for the input used in the report and provenance record.
    dry_run : bool
        If True, skip actual cell removal, normalization, and h5ad writing.
    auto_thresholds : bool
        If True, override config thresholds with MAD-based suggestions.

    Returns
    -------
    AnnData
        Cleaned AnnData with adata.uns['annqc'] fully populated.
    """
    from annqc import __version__
    from annqc.config import DEFAULT_CONFIG, _deep_merge, get_default_config, load_config
    from annqc.doublets import detect_doublets
    from annqc.filter import apply_filters, filter_genes, flag_cells
    from annqc.qc import calculate_qc_metrics
    from annqc.spec import init_record
    from annqc.utils import ensure_dir

    # --- Step 1: Load input ---
    if isinstance(adata_or_path, (str, os.PathLike)):
        path = str(adata_or_path)
        if not os.path.exists(path):
            raise FileNotFoundError(f"Input file not found: {path}")
        if input_file is None:
            input_file = path
        if os.path.isdir(path):
            logger.info("Loading 10x directory: %s", path)
            adata = sc.read_10x_mtx(path, var_names="gene_symbols", cache=True)
        elif path.endswith(".h5"):
            logger.info("Loading 10x HDF5: %s", path)
            adata = sc.read_10x_h5(path)
        else:
            logger.info("Loading h5ad: %s", path)
            adata = sc.read_h5ad(path)
    else:
        adata = adata_or_path.copy()

    if input_file is None:
        input_file = "in-memory AnnData"

    logger.info("Input: %d cells x %d genes", adata.n_obs, adata.n_vars)

    # --- Step 2: Resolve config ---
    if config is None:
        cfg = get_default_config()
    elif isinstance(config, (str, os.PathLike)):
        cfg = load_config(str(config))
    elif isinstance(config, dict):
        cfg = _deep_merge(DEFAULT_CONFIG, config)
    else:
        raise TypeError(f"config must be dict, str/Path, or None; got {type(config)}")

    # --- Step 3: Initialize provenance record ---
    now = datetime.now(timezone.utc).replace(tzinfo=None).isoformat(timespec="seconds")
    record = init_record(
        version=__version__,
        date=now,
        seed=seed,
        config=copy.deepcopy(cfg),
        input_file=input_file,
    )
    adata.uns["annqc"] = record
    record["cell_counts"]["input"] = adata.n_obs
    record["dry_run"] = dry_run

    # --- Step 4: Calculate QC metrics ---
    adata = calculate_qc_metrics(
        adata,
        mito_prefix=cfg["mito"]["prefix"],
        ribo_prefix=cfg["ribo"]["prefix"],
    )

    # Snapshot of pre-filtering QC metrics — used by report for before/after plots
    _raw = {}
    for col in ("pct_counts_mt", "n_genes_by_counts", "total_counts", "pct_counts_ribo"):
        if col in adata.obs.columns:
            _raw[col] = adata.obs[col].tolist()
    record["raw_obs_metrics"] = _raw

    # --- Auto-thresholds (MAD-based) ---
    if auto_thresholds:
        from annqc.thresholds import suggest_thresholds
        suggested = suggest_thresholds(adata)["standard"]
        # Override config thresholds with MAD-based suggestions
        cfg["mito"]["max_pct"] = suggested.get("mito_max_pct") or cfg["mito"]["max_pct"]
        cfg["cells"]["min_genes"] = suggested.get("min_genes") or cfg["cells"]["min_genes"]
        cfg["cells"]["max_genes"] = suggested.get("max_genes") or cfg["cells"]["max_genes"]
        cfg["cells"]["min_counts"] = suggested.get("min_counts") or cfg["cells"]["min_counts"]
        cfg["cells"]["max_counts"] = suggested.get("max_counts")
        record["threshold_method"] = "auto_mad"
        record["thresholds"]["method"] = "auto_mad"
        logger.info("Auto-thresholds (MAD-based): %s", suggested)
    else:
        record["thresholds"]["method"] = "manual"

    # --- Step 5: First-pass cell flagging (no doublets yet) ---
    adata = flag_cells(adata, cfg)
    reasons = adata.obs["annqc_filter_reason"].values

    def _after_removing(bad_reasons):
        """Count cells not flagged by any of bad_reasons."""
        return int(sum(1 for r in reasons if r not in bad_reasons))

    record["cell_counts"]["after_mito_filter"] = _after_removing({"mito"})
    record["cell_counts"]["after_gene_filter"] = _after_removing(
        {"mito", "min_genes", "max_genes"}
    )
    record["cell_counts"]["after_count_filter"] = _after_removing(
        {"mito", "min_genes", "max_genes", "min_counts", "max_counts"}
    )

    # --- Step 6: Gene-level filtering (before cell removal) ---
    adata = filter_genes(adata, min_cells=cfg["genes"]["min_cells"])
    # sc.pp.filter_genes replaces adata's internal arrays, breaking the record reference
    record = adata.uns["annqc"]

    # --- Step 7: Doublet detection ---
    if no_doublet_detection:
        adata.obs["annqc_doublet_score"] = float("nan")
        adata.obs["annqc_is_doublet"] = False
        record["doublet_status"] = "SKIPPED"
        record["thresholds"]["doublet_threshold"] = float("nan")
        logger.info("Doublet detection skipped (--no-doublet-detection)")
    else:
        adata = detect_doublets(
            adata,
            threshold=cfg["doublets"]["threshold"],
            simulate_doublet_ratio=cfg["doublets"]["simulate_doublet_ratio"],
            seed=seed,
        )

    # --- Step 8: Re-flag cells including doublets ---
    # (runs regardless of whether doublet detection was skipped)
    adata = flag_cells(adata, cfg)
    reasons = adata.obs["annqc_filter_reason"].values

    record["cell_counts"]["after_doublet_filter"] = _after_removing(
        {"mito", "min_genes", "max_genes", "min_counts", "max_counts", "doublet"}
    )

    # --- Step 9: Remove flagged cells (or simulate in dry run) ---
    if dry_run:
        # In dry run, do not actually filter. Populate output counts from flags.
        n_would_keep = int(adata.obs["annqc_pass"].sum())
        record["cell_counts"]["output"] = n_would_keep
        logger.info(
            "DRY RUN — would keep %d / %d cells (no cells actually removed)",
            n_would_keep,
            adata.n_obs,
        )
    else:
        adata = apply_filters(adata)
        # apply_filters returns adata[mask].copy(), which deep-copies uns — re-sync record
        record = adata.uns["annqc"]
        record["cell_counts"]["output"] = adata.n_obs
        logger.info(
            "Filtering complete: %d / %d cells retained",
            adata.n_obs,
            record["cell_counts"]["input"],
        )

    # --- Step 10: Normalization ---
    if not dry_run:
        # Preserve raw counts before any normalization
        adata.layers["counts"] = adata.X.copy()
        record["raw_counts_layer"] = "counts"
        logger.info("Raw counts saved to adata.layers['counts']")

        if cfg["normalization"]["method"] == "log1p":
            target_sum = cfg["normalization"]["target_sum"]
            sc.pp.normalize_total(adata, target_sum=target_sum)
            sc.pp.log1p(adata)
            logger.info("Normalized to %d counts/cell + log1p", target_sum)
            # sc.pp.log1p may replace internal arrays, breaking the record reference
            record = adata.uns["annqc"]

    # --- Step 11: Finalize thresholds ---
    record["thresholds"]["mito_max_pct"] = cfg["mito"]["max_pct"]
    record["thresholds"]["min_genes"] = cfg["cells"].get("min_genes")
    record["thresholds"]["max_genes"] = cfg["cells"].get("max_genes")
    record["thresholds"]["min_counts"] = cfg["cells"].get("min_counts")
    # doublet_threshold was set by detect_doublets; guard against missing key
    if record["thresholds"].get("doublet_threshold") is None:
        record["thresholds"]["doublet_threshold"] = float("nan")

    # --- Step 11.5: Distribution analysis and threshold explanations ---
    try:
        from annqc.decisions import analyze_distributions
        analyze_distributions(adata)
        # re-sync record in case adata.uns was modified
        record = adata.uns["annqc"]
    except Exception as exc:
        logger.warning("Distribution analysis failed: %s", exc)

    # --- Step 12: Per-sample stats ---
    if sample_key is not None:
        if sample_key in adata.obs.columns:
            record["per_sample"] = _compute_per_sample(adata, sample_key, cfg)
        else:
            logger.warning(
                "sample_key '%s' not found in adata.obs — skipping per-sample mode",
                sample_key,
            )

    # --- Step 13: Warnings and status ---
    _generate_warnings(record, cfg)
    _set_status(record, cfg)
    logger.info("Pipeline complete. Status: %s", record["status"])

    # --- Step 14: Save cleaned h5ad ---
    if output is not None and not dry_run:
        out_dir = os.path.dirname(os.path.abspath(output))
        ensure_dir(out_dir)
        adata.write_h5ad(output)
        logger.info("Saved cleaned AnnData → %s", output)
    if dry_run:
        logger.info("DRY RUN — no output .h5ad written")

    # --- Step 15: Build HTML report ---
    if report_path is not None:
        try:
            from annqc.report.builder import build_report
            build_report(adata, output_path=report_path)
            logger.info("Saved HTML report → %s", report_path)
        except Exception as exc:
            logger.warning("Report generation failed: %s", exc)

    return adata


def _compute_per_sample(adata, sample_key: str, cfg: dict) -> dict:
    """Compute per-sample QC statistics from the cleaned AnnData."""
    per_sample = {}
    mito_threshold = cfg["mito"].get("max_pct", 20)
    thr_cfg = cfg.get("thresholds", {})
    min_cells_pass_per_sample = thr_cfg.get("min_cells_pass_per_sample", 200)
    min_cells_warn_per_sample = thr_cfg.get("min_cells_warn_per_sample", 500)
    for sample in adata.obs[sample_key].unique():
        mask = adata.obs[sample_key] == sample
        sub = adata[mask]
        n_out = int(mask.sum())

        median_genes = (
            float(np.median(sub.obs["n_genes_by_counts"].values))
            if "n_genes_by_counts" in sub.obs.columns
            else None
        )
        median_mito = (
            float(np.median(sub.obs["pct_counts_mt"].values))
            if "pct_counts_mt" in sub.obs.columns
            else None
        )
        doublet_pct = (
            100.0 * float(sub.obs["annqc_is_doublet"].sum()) / max(n_out, 1)
            if "annqc_is_doublet" in sub.obs.columns
            else None
        )

        if n_out < min_cells_pass_per_sample:
            sample_status = "FAIL"
        elif median_mito is not None and median_mito > mito_threshold * 0.85:
            sample_status = "WARN"
        elif n_out < min_cells_warn_per_sample:
            sample_status = "WARN"
        else:
            sample_status = "PASS"

        per_sample[str(sample)] = {
            "input": n_out,
            "output": n_out,
            "median_genes": median_genes,
            "median_mito_pct": median_mito,
            "doublet_pct": doublet_pct,
            "status": sample_status,
        }
    return per_sample


def _generate_warnings(record: dict, cfg: dict) -> None:
    """Append auto-generated warnings to record['warnings']."""
    warnings = record.setdefault("warnings", [])
    cc = record["cell_counts"]
    n_input = cc.get("input", 0)
    n_output = cc.get("output", 0)
    mito_threshold = cfg["mito"].get("max_pct", 20)

    # Very few cells remain
    min_cells_warn = cfg.get("thresholds", {}).get("min_cells_warn", 500)
    if n_output < min_cells_warn:
        warnings.append(
            f"⚠️ Only {n_output} cells remain after filtering. "
            f"Consider relaxing thresholds or checking input data quality."
        )

    # Mito filter too aggressive
    after_mito = cc.get("after_mito_filter", n_input)
    if n_input > 0:
        pct_removed_mito = (n_input - after_mito) / n_input
        if pct_removed_mito > 0.15:
            warnings.append(
                f"⚠️ Mito filter removed {pct_removed_mito:.1%} of cells — this is unusually high. "
                f"Consider raising mito.max_pct from {mito_threshold} or checking library quality."
            )

    # High doublet rate
    cells_before_doublet = cc.get("after_count_filter", 0)
    cells_after_doublet = cc.get("after_doublet_filter", 0)
    if cells_before_doublet > 0:
        doublet_rate = (cells_before_doublet - cells_after_doublet) / cells_before_doublet
        if doublet_rate > 0.08:
            warnings.append(
                f"⚠️ Doublet rate is {doublet_rate:.1%} — expected range is 1–5% for 10x. "
                f"This may indicate overloading of the microfluidics chip."
            )
        # Suspiciously low doublet rate
        if doublet_rate < 0.001 and n_input > 1000:
            warnings.append(
                f"⚠️ Doublet rate is unusually low ({doublet_rate:.1%}). "
                f"Verify that Scrublet ran correctly."
            )

    # Extreme gene count filtering
    n_genes_removed_pct = record.get("_genes_removed_pct")  # set by pipeline if available
    # (approximate from cell counts isn't reliable for gene counts — skip unless pipeline sets it)

    # Sample-level warnings
    for sample, stats in record.get("per_sample", {}).items():
        mito_pct = stats.get("median_mito_pct")
        if mito_pct is not None and mito_pct > mito_threshold * 0.85:
            warnings.append(
                f"⚠️ Sample '{sample}': median mito% is {mito_pct:.1f} "
                f"— approaching threshold of {mito_threshold}. Monitor this sample."
            )
        if stats.get("status") == "FAIL":
            warnings.append(
                f"⚠️ Sample '{sample}' has status FAIL — check filtering aggressiveness for this sample."
            )


def _set_status(record: dict, cfg: dict) -> None:
    """Set record['status'] to PASS, WARN, or FAIL."""
    min_cells_pass = cfg.get("thresholds", {}).get("min_cells_pass", 100)
    cc = record["cell_counts"]
    if cc["output"] < min_cells_pass:
        record["status"] = "FAIL"
        return
    for stats in record.get("per_sample", {}).values():
        if stats.get("status") == "FAIL":
            record["status"] = "FAIL"
            return
    # Never mark PASS if doublet detection failed
    if record.get("doublet_status") == "FAILED":
        record["status"] = "INCOMPLETE"
        return
    record["status"] = "PASS"
