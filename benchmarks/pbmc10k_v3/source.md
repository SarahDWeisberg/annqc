# PBMC 10k v3 — Dataset Source & Published QC

## Data Source

**Dataset:** 10,000 PBMCs from a Healthy Donor, Chromium 3' v3 Chemistry
**Provider:** 10x Genomics
**Download URL:** https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5
**Dataset page:** https://www.10xgenomics.com/datasets/10-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-1-standard-3-0-0

## Published QC Metrics

**Source:** VSN-Pipelines documentation — "PBMC 10k tutorial"
**URL:** https://vsn-pipelines-examples.readthedocs.io/en/latest/PBMC10k.html

**Additional source (NBIS Scanpy QC workshop, 2020):**
**URL:** https://nbisweden.github.io/workshop-archive/workshop-scRNAseq/2020-01-27/labs/compiled/scanpy/scanpy_01_qc.html

| QC Metric          | Published Cutoff |
|--------------------|-----------------|
| min_genes          | 200             |
| max_genes          | 4,000           |
| max_mito_pct       | 15%             |
| min_cells per gene | 3               |
| min_counts         | not specified   |

**Cells before filtering:** ~10,194 (Cell Ranger filtered output)
**Cells after filtering:** not explicitly reported in source
