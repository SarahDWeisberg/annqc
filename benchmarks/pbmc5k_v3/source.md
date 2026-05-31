# 5k PBMC v3 — Source & Published QC

## Data Source

**Dataset:** 5,000 Peripheral Blood Mononuclear Cells from a Healthy Donor, Chromium 3' v3
**Organism:** Homo sapiens
**Provider:** 10x Genomics (reference dataset)
**Download URL:** https://cf.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_v3/5k_pbmc_v3_filtered_feature_bc_matrix.h5
**Dataset page:** https://www.10xgenomics.com/datasets/5-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-1-standard-3-0-2

## Published QC Metrics

**Source:** NBIS Scanpy QC Workshop (2020) — applied to this exact dataset
**URL:** https://nbisweden.github.io/workshop-archive/workshop-scRNAseq/2020-01-27/labs/compiled/scanpy/scanpy_01_qc.html

| QC Metric          | Published Cutoff      |
|--------------------|-----------------------|
| min_genes          | 1,000                 |
| max_genes          | 4,100                 |
| max_mito_pct       | 25%                   |
| min_ribo_pct       | 5%                    |
| min_cells per gene | 3                     |
| min_counts         | not specified         |

**Cells before filtering:** ~2,931 (after Cell Ranger cell calling)
**Cells after filtering:** ~2,527 (~404 removed)

**Note:** The NBIS workshop explicitly targets v3 chemistry, hence the higher
max_genes (4,100) and max_mito_pct (25%) compared to v2 datasets.
