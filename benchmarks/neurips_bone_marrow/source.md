# NeurIPS 2021 Bone Marrow Multiome — Source & Published QC

## Data Source

**Dataset:** NeurIPS 2021 Multimodal Single Cell Data Integration Challenge — bone marrow (RNA only)
**Organism:** Homo sapiens
**Tissue:** Bone marrow
**Samples:** site1-donor1 (s1d1) + site1-donor3 (s1d3)
**Figshare:** https://figshare.com/articles/dataset/NeurIPS_2021_Multimodal_Single_Cell_Data_Integration_Challenge/22716739
**Download URLs:**
  s1d1: https://ndownloader.figshare.com/files/40347877
  s1d3: https://ndownloader.figshare.com/files/40347880

## Published QC Metrics

**Source:** sc-best-practices book, Quality Control chapter (Theis lab)
**URL:** https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html

| QC Metric             | Published Value                              |
|-----------------------|----------------------------------------------|
| max_mito_pct          | 8% (hard cap)                                |
| log1p_total_counts    | 5 MADs from median (adaptive)                |
| log1p_n_genes         | 5 MADs from median (adaptive)                |
| pct_counts_top20      | 5 MADs from median (adaptive)                |
| min_cells per gene    | 20                                           |

**Cells before filtering:** 16,934 (sc-best-practices; our combined h5ad has 17,125 — minor version difference)
**Cells after filtering:** 14,814

**Note:** The sc-best-practices approach combines a hard mito cap (8%) with adaptive MAD-based thresholds
on log1p-transformed metrics — different from both annqc manual (hard thresholds) and annqc auto (linear
MAD). The min_cells=20 per gene is more aggressive than the typical 3. annqc manual config uses max_mito=8%
matching the hard cap; gene/count thresholds are standard fallbacks since the published method uses adaptive values.
