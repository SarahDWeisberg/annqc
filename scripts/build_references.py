"""Build AnnQC reference profiles from CELLxGENE census.

Run once to generate reference profile JSONs in annqc/reference/profiles/.

Usage:
    python scripts/build_references.py
    python scripts/build_references.py --tissue pbmc --assay 10x_v3
    python scripts/build_references.py --dry-run  # show what would be built

Requirements:
    pip install cellxgene-census

What it does:
    For each (tissue, assay, suspension_type) target defined in TARGETS:
      1. Opens CELLxGENE census (streaming, no full download)
      2. Lists matching dataset_ids filtered by tissue/assay/disease/species
      3. For each dataset, streams only the obs metadata + MT- gene expression
      4. Computes per-cell QC metrics, then stores per-dataset medians
      5. Writes a versioned reference profile JSON

Profile files are small (~5-50 KB each) and should be committed to the repo.
"""

import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

# Add repo root to path so we can import annqc
sys.path.insert(0, str(Path(__file__).parent.parent))

from annqc.reference.schema import (
    METRICS, confidence_level, normalize_assay, save_profile, SCHEMA_VERSION,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Target definitions
# Each entry: (tissue_general, raw_assay_string, suspension_type)
# tissue_general must match CELLxGENE taxonomy; assay must match obs.assay
# ---------------------------------------------------------------------------
TARGETS = [
    # PBMC — standard benchmark tissue, whole-cell, common assays
    ("blood", "10x 3' v3",  "cell"),
    ("blood", "10x 3' v2",  "cell"),
    # Kidney
    ("kidney", "10x 3' v3", "cell"),
    ("kidney", "10x 3' v2", "cell"),
    # Brain — single nucleus
    ("brain",  "10x 3' v3", "nucleus"),
]

SPECIES = "Homo sapiens"
DISEASE_FILTER = "normal"            # restrict to healthy tissue
MAX_CELLS_PER_DATASET = 50_000       # subsample large datasets for speed
MIN_CELLS_PER_DATASET = 200          # skip tiny datasets


def _compute_mito_pct(obs_df, X_mt, total_counts_col="raw_sum"):
    """Compute per-cell mitochondrial % from obs table and MT gene expression matrix."""
    mt_sum = np.asarray(X_mt.sum(axis=1)).ravel()
    total = obs_df[total_counts_col].values.astype(float)
    total = np.where(total > 0, total, 1.0)  # avoid div-by-zero
    return 100.0 * mt_sum / total


def _dataset_medians(cells_df, mt_pct_arr):
    """Return per-dataset-level median statistics from a cells obs DataFrame."""
    medians = {}

    # pct_counts_mt
    if mt_pct_arr is not None and len(mt_pct_arr) > 0:
        medians["pct_counts_mt"] = float(np.median(mt_pct_arr))

    # n_genes_by_counts — census stores as 'nnz' or 'n_measured_vars'
    for col in ("nnz", "n_measured_vars"):
        if col in cells_df.columns:
            medians["n_genes_by_counts"] = float(cells_df[col].median())
            break

    # total_counts — census stores as 'raw_sum'
    if "raw_sum" in cells_df.columns:
        medians["total_counts"] = float(cells_df["raw_sum"].median())

    return medians


def _distribution_summary(values: list[float]) -> dict:
    """Compute summary stats for a list of per-dataset medians."""
    if not values:
        return {}
    arr = np.array(values, dtype=float)
    mad = float(np.median(np.abs(arr - np.median(arr))))
    return {
        "median": round(float(np.median(arr)), 4),
        "q25":    round(float(np.percentile(arr, 25)), 4),
        "q75":    round(float(np.percentile(arr, 75)), 4),
        "mad":    round(mad, 4),
        "p5":     round(float(np.percentile(arr, 5)), 4),
        "p95":    round(float(np.percentile(arr, 95)), 4),
        "n":      len(arr),
    }


def build_profile(
    tissue: str,
    raw_assay: str,
    suspension_type: str,
    census,
    census_version: str,
    dry_run: bool = False,
) -> dict | None:
    """Build a reference profile for one (tissue, assay, suspension_type) target."""
    import cellxgene_census
    import pandas as pd

    assay_key = normalize_assay(raw_assay)
    logger.info(
        "Building profile: tissue=%s  assay=%s (%s)  suspension=%s",
        tissue, assay_key, raw_assay, suspension_type,
    )

    # --- Query obs table ---
    value_filter = (
        f"tissue_general == '{tissue}' "
        f"and assay == '{raw_assay}' "
        f"and suspension_type == '{suspension_type}' "
        f"and disease == '{DISEASE_FILTER}' "
        f"and is_primary_data == True"
    )
    logger.debug("CELLxGENE filter: %s", value_filter)

    try:
        obs_df = (
            census["census_data"]["homo_sapiens"]
            .obs
            .read(
                value_filter=value_filter,
                column_names=[
                    "soma_joinid", "dataset_id", "raw_sum", "nnz",
                    "n_measured_vars", "disease", "development_stage",
                    "donor_id",
                ],
            )
            .concat()
            .to_pandas()
        )
    except Exception as exc:
        logger.error("Census obs query failed: %s", exc)
        return None

    if obs_df.empty:
        logger.warning("No cells found for %s/%s/%s — skipping", tissue, assay_key, suspension_type)
        return None

    logger.info("Found %d total cells across datasets", len(obs_df))

    # Get MT- gene var_ids for mito% calculation
    try:
        var_df = (
            census["census_data"]["homo_sapiens"]
            .ms["RNA"]
            .var
            .read(value_filter="feature_name like 'MT-%'")
            .concat()
            .to_pandas()
        )
        mt_soma_ids = var_df["soma_joinid"].tolist()
        logger.info("Found %d MT- genes", len(mt_soma_ids))
    except Exception as exc:
        logger.warning("Could not retrieve MT- gene IDs: %s", exc)
        mt_soma_ids = []

    # --- Per-dataset processing ---
    dataset_ids = obs_df["dataset_id"].unique().tolist()
    logger.info("Processing %d datasets", len(dataset_ids))

    per_dataset_medians = {m: [] for m in METRICS}
    per_dataset_medians["pct_counts_ribo"] = []  # will be empty (ribo not in census obs)
    n_processed = 0
    n_cells_total = 0
    disease_counts = {}
    dev_stage_counts = {}

    for dataset_id in dataset_ids:
        cells = obs_df[obs_df["dataset_id"] == dataset_id].copy()

        if len(cells) < MIN_CELLS_PER_DATASET:
            logger.debug("Skipping %s: only %d cells", dataset_id, len(cells))
            continue

        # Subsample large datasets
        if len(cells) > MAX_CELLS_PER_DATASET:
            cells = cells.sample(MAX_CELLS_PER_DATASET, random_state=42)
            logger.debug("Subsampled %s to %d cells", dataset_id, MAX_CELLS_PER_DATASET)

        # Compute mito% if MT genes available
        mt_pct_arr = None
        if mt_soma_ids and not dry_run:
            try:
                cell_soma_ids = cells["soma_joinid"].tolist()
                X_mt = (
                    census["census_data"]["homo_sapiens"]
                    .ms["RNA"]
                    .X["raw"]
                    .read(
                        coords=(cell_soma_ids, mt_soma_ids),
                    )
                    .tables()
                    .concat()
                    .to_pandas()
                )
                # Pivot to dense (cells × MT genes)
                import scipy.sparse as sp
                if not X_mt.empty:
                    # Build sparse sum of MT counts per cell
                    cell_id_to_idx = {cid: i for i, cid in enumerate(cell_soma_ids)}
                    mt_sums = np.zeros(len(cell_soma_ids))
                    for _, row in X_mt.iterrows():
                        idx = cell_id_to_idx.get(row["soma_dim_0"])
                        if idx is not None:
                            mt_sums[idx] += row["soma_data"]
                    total = cells["raw_sum"].values.astype(float)
                    total = np.where(total > 0, total, 1.0)
                    mt_pct_arr = 100.0 * mt_sums / total
                else:
                    mt_pct_arr = np.zeros(len(cell_soma_ids))
            except Exception as exc:
                logger.warning("MT% computation failed for %s: %s", dataset_id, exc)
                mt_pct_arr = None

        if dry_run:
            logger.info("  [dry-run] Would process dataset %s (%d cells)", dataset_id, len(cells))
            n_processed += 1
            continue

        medians = _dataset_medians(cells, mt_pct_arr)

        for metric, val in medians.items():
            if metric in per_dataset_medians:
                per_dataset_medians[metric].append(val)

        n_processed += 1
        n_cells_total += len(cells)

        # Track metadata composition
        if "disease" in cells.columns:
            for d in cells["disease"].unique():
                disease_counts[str(d)] = disease_counts.get(str(d), 0) + 1
        if "development_stage" in cells.columns:
            for s in cells["development_stage"].unique():
                dev_stage_counts[str(s)] = dev_stage_counts.get(str(s), 0) + 1

    if dry_run:
        logger.info("[dry-run] Would have processed %d datasets", n_processed)
        return None

    if n_processed == 0:
        logger.warning("No datasets passed minimum cell threshold — aborting")
        return None

    logger.info("Processed %d datasets, %d total cells", n_processed, n_cells_total)

    # Build metrics section (skip metrics with no data)
    metrics_section = {}
    for metric in METRICS:
        vals = per_dataset_medians.get(metric, [])
        if len(vals) >= 5:
            metrics_section[metric] = {
                "dataset_medians": [round(v, 4) for v in vals],
                "summary": _distribution_summary(vals),
            }

    profile = {
        "schema_version": SCHEMA_VERSION,
        "tissue": tissue,
        "assay": assay_key,
        "assay_raw": raw_assay,
        "suspension_type": suspension_type,
        "species": SPECIES,
        "generated_date": datetime.utcnow().date().isoformat(),
        "census_version": census_version,
        "inclusion_criteria": (
            f"tissue_general=='{tissue}', assay=='{raw_assay}', "
            f"suspension_type=='{suspension_type}', disease=='{DISEASE_FILTER}', "
            f"is_primary_data==True, min_cells_per_dataset={MIN_CELLS_PER_DATASET}"
        ),
        "n_datasets": n_processed,
        "n_cells_total": n_cells_total,
        "confidence": confidence_level(n_processed),
        "disease_composition": disease_counts,
        "development_stage_composition": dev_stage_counts,
        "metrics": metrics_section,
    }

    return profile


def main():
    parser = argparse.ArgumentParser(description="Build AnnQC reference profiles from CELLxGENE census")
    parser.add_argument("--tissue", default=None, help="Build only this tissue (optional)")
    parser.add_argument("--assay", default=None, help="Build only this assay canonical key (optional)")
    parser.add_argument("--suspension-type", default=None, help="Build only this suspension type")
    parser.add_argument("--output-dir", default=None, help="Output directory (default: annqc/reference/profiles/)")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be built without querying data")
    parser.add_argument("--census-version", default="stable", help="CELLxGENE census version (default: stable)")
    args = parser.parse_args()

    import cellxgene_census

    targets = TARGETS
    if args.tissue:
        targets = [t for t in targets if t[0] == args.tissue]
    if args.assay:
        targets = [t for t in targets if normalize_assay(t[1]) == args.assay]
    if args.suspension_type:
        targets = [t for t in targets if t[2] == args.suspension_type]

    if not targets:
        logger.error("No targets matched the specified filters")
        sys.exit(1)

    logger.info("Building %d reference profile(s)", len(targets))
    if args.dry_run:
        logger.info("DRY RUN — no data will be downloaded or written")

    written = []
    with cellxgene_census.open_soma(census_version=args.census_version) as census:
        census_version = census.get("census_info", {}).get("census_version", args.census_version)
        if hasattr(census_version, "read"):
            try:
                census_version = str(census["census_info"]["summary"]
                    .read().concat().to_pandas()
                    .set_index("label")["value"]["census_version"])
            except Exception:
                census_version = args.census_version

        for tissue, raw_assay, suspension_type in targets:
            profile = build_profile(
                tissue=tissue,
                raw_assay=raw_assay,
                suspension_type=suspension_type,
                census=census,
                census_version=str(census_version),
                dry_run=args.dry_run,
            )
            if profile is not None:
                path = save_profile(profile, output_dir=args.output_dir)
                logger.info("Wrote: %s (%d datasets, confidence=%s)",
                            path.name, profile["n_datasets"], profile["confidence"])
                written.append(path.name)

    if not args.dry_run:
        logger.info("Done. Wrote %d profile(s): %s", len(written), ", ".join(written))
    else:
        logger.info("Dry run complete. Would build %d profile(s).", len(targets))


if __name__ == "__main__":
    main()
