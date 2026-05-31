# Mouse Small Intestinal Epithelium (Haber et al. 2017) — Source & Published QC

## Data Source

**Dataset:** Mouse small intestinal epithelium single-cell RNA-seq
**Organism:** Mus musculus
**Paper:** Haber et al. (2017) "A single-cell survey of the small intestinal epithelium"
**Journal:** Nature, DOI: 10.1038/nature24489
**GEO accession:** GSE92332
**Download URL (Zenodo repackage as h5ad):** https://zenodo.org/records/4447233/files/haber_raw.h5ad?download=1
**Zenodo record:** https://zenodo.org/records/4447233

## Published QC Metrics

**Source:** Haber et al. 2017 Methods section (PMC6022292) + Luecken & Theis 2019 tutorial
**Paper URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC6022292/
**Tutorial URL:** https://github.com/theislab/single-cell-tutorial

| QC Metric          | Paper (Haber 2017)                 | Tutorial (Luecken & Theis 2019)  |
|--------------------|------------------------------------|----------------------------------|
| min_genes          | 800                                | 500                              |
| max_genes          | top 1% by genes removed (doublets) | not stated                       |
| max_mito_pct       | not stated numerically             | 20%                              |
| min_counts         | not stated                         | not stated                       |
| min_cells per gene | not stated                         | 3                                |

**Cells before filtering:** 8,882 (droplet-based subset)
**Cells after filtering:** 7,216 (after removing 1,402 low-quality + 264 contaminating cells)

**Note:** The Luecken & Theis 2019 "Current best practices" paper (Mol. Syst. Biol. 15:e8746)
uses this dataset in their tutorial notebook with min_genes=500, max_mito=20%.
annqc config uses those tutorial thresholds since the original paper did not state mito cutoff.
