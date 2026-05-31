# Wu et al. 2019 Human Kidney — Source & Published QC

## Data Source

**Dataset:** Single-cell RNA sequencing of human kidney (3 donors)
**Organism:** Homo sapiens
**Tissue:** Kidney
**Paper:** Wu et al. (2019) "Single-cell RNA sequencing of human kidney"
**Journal:** Scientific Data, DOI: 10.1038/s41597-019-0351-8
**GEO accession:** GSE131685
**Figshare:** https://figshare.com/articles/Single-cell_RNA_sequencing_of_human_kidney/8131328
**GEO FTP:** https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131685/suppl/

## Published QC Metrics

**Source:** Wu et al. 2019 Methods section
**Paper URL:** https://www.nature.com/articles/s41597-019-0351-8

| QC Metric          | Published Value    |
|--------------------|--------------------|
| min_genes          | 200                |
| max_genes          | 2,500              |
| max_mito_pct       | 30%                |
| min_counts         | not stated         |
| min_cells per gene | not stated         |

**Cells before filtering:** 25,404 (kidney1=8164 + kidney2=6499 + kidney3=10741)
**Cells after filtering:** 23,366

**Note:** 30% mito threshold is notably high — kidney has elevated mitochondrial content vs blood.
This dataset is a good test case for the mito floor logic (annqc auto-MAD should not drop below 10%).
The paper removed cells with > 40,000 UMIs as likely doublets; annqc handles doublets via Scrublet.
