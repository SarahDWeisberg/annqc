# Mouse Thymic T-cells — Dataset Source & Published QC

## Data Source

**Dataset:** Mouse thymic T-cell development scRNA-seq (Bacon et al. 2018)
**Provider:** Zenodo (Galaxy Training Network)
**Download URL:** https://zenodo.org/records/7053673/files/Mito-counted_AnnData
**Zenodo record:** https://zenodo.org/records/7053673

## Published QC Metrics

**Source:** Galaxy Training Network — "Filter, plot and explore single-cell RNA-seq data (Scanpy)"
**URL:** https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-case-jupyter_basic-pipeline/tutorial.html

| QC Metric          | Published Cutoff        | Notes                              |
|--------------------|-------------------------|------------------------------------|
| min_genes (log1p)  | ≥ 5.7 (~300 genes)      | Applied in log1p space             |
| max_genes (log1p)  | ≤ 20.0                  |                                    |
| min_counts (log1p) | ≥ 6.3 (~500 UMIs)       | Applied in log1p space             |
| max_counts (log1p) | ≤ 20.0                  |                                    |
| max_mito_pct       | ≤ 4.5%                  |                                    |
| min_cells per gene | 3                       |                                    |

**Cells before filtering:** 31,178
**Cells after filtering:** 8,605
**% retained:** ~27.6%

**Original paper:** Bacon et al. (2018), murine thymic single-cell RNA-seq
**Species:** Mus musculus (mouse — use "mt-" prefix for mito genes)
