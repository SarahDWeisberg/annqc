# PBMC 1k v3 — Dataset Source & Published QC

## Data Source

**Dataset:** 1,000 PBMCs from a Healthy Donor, Chromium 3' v3 Chemistry
**Provider:** 10x Genomics
**Download URL:** https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.h5
**Dataset page:** https://www.10xgenomics.com/datasets/1-k-pb-mc-s-from-a-healthy-donor-v-3-chemistry-3-1-standard-3-0-0

## Published QC Metrics

**Source:** Seurat v5 tutorial — "Guided Clustering Tutorial (PBMC)"
**URL:** https://satijalab.org/seurat/articles/pbmc3k_tutorial

**Also referenced in:** 10x Genomics Cell Ranger documentation
**URL:** https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-overview

| QC Metric          | Published Cutoff |
|--------------------|-----------------|
| min_genes          | 200             |
| max_genes          | 2,500           |
| max_mito_pct       | 5%              |
| min_cells per gene | 3               |
| min_counts         | not specified   |

**Cells before filtering:** ~1,222 (Cell Ranger filtered output)
**Cells after filtering:** not explicitly reported for 1k; Seurat tutorial on 3k retains ~2,638/2,700

**Note:** The Seurat tutorial was developed on PBMC 3k data, but identical thresholds are applied
across the v3 chemistry PBMC series (1k, 3k, 10k).
