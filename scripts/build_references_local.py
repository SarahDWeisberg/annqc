"""Build AnnQC reference profiles from local h5ad files.

Architecture-independent alternative to build_references.py (which requires
tiledbsoma / CELLxGENE census). Use this to seed profiles from your own
processed datasets or downloaded h5ad files.

Usage:
    # Build a kidney profile from a directory of h5ad files
    python scripts/build_references_local.py \\
        --tissue kidney --assay 10x_v3 --suspension-type cell \\
        --files data/kidney_ref_*.h5ad

    # Or pass a JSON manifest
    python scripts/build_references_local.py --manifest manifests/kidney_refs.json

Manifest format (JSON list of objects):
    [
      {"path": "data/wu_kidney.h5ad", "dataset_id": "wu2019", "disease": "normal"},
      {"path": "data/lake_kidney.h5ad", "dataset_id": "lake2023", "disease": "normal"},
      ...
    ]
"""

import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from annqc.reference.schema import (
    METRICS, confidence_level, normalize_assay, save_profile, SCHEMA_VERSION,
)
from annqc.qc import calculate_qc_metrics

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

MIN_CELLS = 200


def _process_one(path: str, dataset_id: str, mito_prefix: str = "MT-") -> dict | None:
    """Load one h5ad or 10x HDF5 file, compute QC metrics, return per-dataset medians."""
    import scanpy as sc

    p = Path(path)
    try:
        if p.suffix == ".h5":
            adata = sc.read_10x_h5(path)
        else:
            adata = sc.read_h5ad(path)
        adata.var_names_make_unique()
    except Exception as exc:
        logger.warning("Could not load %s: %s", path, exc)
        return None

    if adata.n_obs < MIN_CELLS:
        logger.warning("Skipping %s: only %d cells", dataset_id, adata.n_obs)
        return None

    try:
        adata = calculate_qc_metrics(adata, mito_prefix=mito_prefix)
    except Exception as exc:
        logger.warning("QC metrics failed for %s: %s", dataset_id, exc)
        return None

    medians = {}
    for col in METRICS:
        if col in adata.obs.columns:
            vals = adata.obs[col].dropna().values
            if len(vals) > 0:
                medians[col] = float(np.median(vals))

    logger.info(
        "  %s: %d cells — %s",
        dataset_id,
        adata.n_obs,
        {k: f"{v:.2g}" for k, v in medians.items()},
    )
    return medians


def _build_summary(vals: list[float]) -> dict:
    a = np.array(vals, dtype=float)
    mad = float(np.median(np.abs(a - np.median(a))))
    return {
        "median": round(float(np.median(a)), 4),
        "q25":    round(float(np.percentile(a, 25)), 4),
        "q75":    round(float(np.percentile(a, 75)), 4),
        "mad":    round(mad, 4),
        "p5":     round(float(np.percentile(a, 5)), 4),
        "p95":    round(float(np.percentile(a, 95)), 4),
        "n":      len(vals),
    }


def main():
    parser = argparse.ArgumentParser(description="Build reference profiles from local h5ad files")
    parser.add_argument("--tissue", required=True, help="Tissue name (e.g. kidney, blood, brain)")
    parser.add_argument("--assay", required=True, help="Assay canonical key (e.g. 10x_v3, 10x_v2)")
    parser.add_argument("--suspension-type", default="cell", choices=["cell", "nucleus"])
    parser.add_argument("--files", nargs="+", metavar="PATH", help="h5ad files to include")
    parser.add_argument("--manifest", default=None, help="JSON manifest file (alternative to --files)")
    parser.add_argument("--mito-prefix", default="MT-", help="Mitochondrial gene prefix")
    parser.add_argument("--output-dir", default=None, help="Output directory (default: annqc/reference/profiles/)")
    parser.add_argument("--species", default="Homo sapiens")
    args = parser.parse_args()

    if not args.files and not args.manifest:
        parser.error("Provide either --files or --manifest")

    # Load file list
    if args.manifest:
        with open(args.manifest) as fh:
            entries = json.load(fh)
    else:
        entries = [{"path": f, "dataset_id": Path(f).stem, "disease": "normal"} for f in args.files]

    logger.info("Processing %d datasets for %s/%s/%s",
                len(entries), args.tissue, args.assay, args.suspension_type)

    per_metric = {m: [] for m in METRICS}
    n_processed = 0
    n_cells_total = 0

    for entry in entries:
        path = entry["path"]
        dataset_id = entry.get("dataset_id", Path(path).stem)
        logger.info("Processing: %s", dataset_id)
        medians = _process_one(path, dataset_id, mito_prefix=args.mito_prefix)
        if medians is None:
            continue
        for metric, val in medians.items():
            if metric in per_metric:
                per_metric[metric].append(val)
        n_processed += 1

    if n_processed == 0:
        logger.error("No datasets successfully processed")
        sys.exit(1)

    logger.info("Processed %d / %d datasets", n_processed, len(entries))

    metrics_section = {}
    for metric in METRICS:
        vals = per_metric.get(metric, [])
        if len(vals) >= 2:
            metrics_section[metric] = {
                "dataset_medians": [round(v, 4) for v in vals],
                "summary": _build_summary(vals),
            }

    profile = {
        "schema_version": SCHEMA_VERSION,
        "tissue": args.tissue,
        "assay": normalize_assay(args.assay),
        "assay_raw": args.assay,
        "suspension_type": args.suspension_type,
        "species": args.species,
        "generated_date": datetime.utcnow().date().isoformat(),
        "census_version": "local",
        "inclusion_criteria": f"local h5ad files, tissue={args.tissue}, assay={args.assay}",
        "n_datasets": n_processed,
        "n_cells_total": n_cells_total,
        "confidence": confidence_level(n_processed),
        "disease_composition": {},
        "development_stage_composition": {},
        "metrics": metrics_section,
    }

    path = save_profile(profile, output_dir=args.output_dir)
    logger.info("Wrote: %s (%d datasets, confidence=%s)",
                path.name, n_processed, profile["confidence"])


if __name__ == "__main__":
    main()
