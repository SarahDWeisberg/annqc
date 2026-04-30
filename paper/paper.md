---
title: 'AnnQC: Reproducible, Auditable Quality Control for Single-Cell RNA-seq Data'
tags:
  - Python
  - single-cell RNA-seq
  - quality control
  - bioinformatics
  - reproducibility
authors:
  - name: Sarah D. Weisberg
    affiliation: 1
affiliations:
  - name: Independent Researcher
    index: 1
date: 2026-04-29
bibliography: paper.bib
---

# Summary

AnnQC is a Python package for reproducible, auditable quality control of single-cell RNA-seq (scRNA-seq) data. It wraps Scanpy's [@Wolf2018] QC primitives and adds configurable filtering thresholds, automatic MAD-based threshold selection, doublet detection via Scrublet [@Wolock2019], and a self-contained HTML report with per-metric visualizations, waterfall charts, and before/after distributions. It accepts `.h5ad`, 10x HDF5, and 10x Genomics MTX directory inputs.

A distinguishing feature of AnnQC is its provenance model: every threshold applied, software version used, cell count at each filtering step, and auto-generated methods text is recorded in `adata.uns['annqc']` and travels with the output file. This allows QC runs to be audited, compared, and reproduced exactly from the saved AnnData object alone.

AnnQC also provides a dedicated sensitivity analysis module (`annqc sensitivity`) that quantifies how cell retention and cluster composition change across a range of threshold values, enabling researchers to justify their QC choices with data rather than convention.

# Statement of Need

Python scRNA-seq workflows commonly lack a unified framework for documenting, explaining, and auditing QC decisions in a reproducible and shareable format. Tools such as Scanpy [@Wolf2018] provide excellent filtering primitives but do not record which thresholds were applied, why they were chosen, or how sensitive downstream results are to those choices. In practice, QC parameters are often copied from tutorials or set by convention without dataset-specific justification, and the thresholds used are frequently not reported in sufficient detail to reproduce a published analysis. AnnQC complements rather than replaces existing tools; it wraps Scanpy's QC primitives and adds provenance tracking, MAD-based threshold selection, and sensitivity analysis.

Separate tools exist for upstream steps that AnnQC intentionally does not cover. Ambient RNA correction (removal of background mRNA contamination) is addressed by tools such as CellBender [@Fleming2023] and SoupX [@Young2020]; these should be applied before AnnQC if ambient RNA removal is required.

A critical but underappreciated challenge in scRNA-seq QC is threshold sensitivity — the degree to which downstream biological conclusions depend on specific QC parameter choices. AnnQC addresses this with a dedicated sensitivity analysis module (`annqc sensitivity`) that quantifies the impact of threshold variation on cell retention and cluster composition, enabling researchers to justify their QC choices with data rather than convention.

# Implementation

**QC metric calculation.** AnnQC wraps `scanpy.pp.calculate_qc_metrics` to compute per-cell mitochondrial content (`pct_counts_mt`), ribosomal content (`pct_counts_ribo`), gene count (`n_genes_by_counts`), and total UMI count (`total_counts`). Gene annotations (mitochondrial, ribosomal) are derived from configurable prefix patterns.

**Threshold configuration.** Thresholds can be set manually via a YAML config file or computed automatically using a MAD-based (Median Absolute Deviation) method. Automatic thresholds are derived as `median ± N×MAD` for each metric, with three preset stringency levels (strict: 3-MAD, standard: 5-MAD, permissive: 7-MAD). The `annqc suggest` command allows threshold exploration before committing to a run.

**Cell filtering.** `flag_cells()` applies thresholds in a fixed priority order (mitochondrial content → gene count → UMI count → doublets), recording the first failing reason for each cell. This ensures waterfall cell counts are interpretable and non-overlapping. The `annqc run --dry-run` flag previews filtering without modifying data.

**Doublet detection.** Scrublet [@Wolock2019] is run after gene-level filtering so that low-variance genes do not inflate PCA-based doublet scores. The doublet threshold is determined automatically by Scrublet's bimodal score distribution fitting, or can be overridden. Scrublet failures are caught and logged without aborting the pipeline.

**Normalization.** After cell filtering, counts are normalized to 10,000 total counts per cell and log1p-transformed (standard Scanpy workflow). Normalization can be disabled for downstream integration pipelines.

**Sensitivity analysis.** The `annqc sensitivity` command sweeps a range of threshold values for each QC metric independently, computing cells removed and percentage removed at each value. If a full pipeline run succeeds, basic clustering (PCA → neighbors → Leiden) is performed on the cleaned data, and cluster impact is computed: for each threshold tested, what fraction of each cluster's cells would have been removed. This identifies clusters that are disproportionately sensitive to threshold choices — often indicating biologically meaningful subpopulations rather than low-quality cells.

**Report generation.** AnnQC renders a self-contained Jinja2/Plotly HTML report containing per-metric violin plots, before/after distributions, filtering waterfall, doublet score histogram, plain-English threshold explanations, and a ready-to-paste academic methods paragraph. A comparison report (`annqc compare`) supports side-by-side evaluation of two AnnQC runs.

**Provenance.** All pipeline inputs, thresholds, cell counts at each step, software versions, warnings, and PASS/FAIL status are stored in `adata.uns['annqc']`. This record is preserved in the output `.h5ad` file and used to auto-generate methods text via `annqc methods`.

# Validation

AnnQC includes a comprehensive test suite covering all pipeline components, including QC metric calculation, threshold configuration and validation, cell and gene filtering, doublet detection, normalization, report generation, and sensitivity analysis. Tests are run against both synthetic AnnData objects and the PBMC 3k dataset [@Zheng2017]. Sensitivity analysis tests verify that threshold sweep results are structurally correct (correct number of values, integer cell counts, percentages in range) and that the full sensitivity pipeline runs without error on standard datasets.

# Limitations

AnnQC is a cell QC and filtering tool and does not address all upstream or downstream processing steps.

1. **Doublet detection.** Only Scrublet [@Wolock2019] is currently supported. Alternative methods such as DoubletFinder or scDblFinder are not integrated.

2. **Ambient RNA correction.** AnnQC does not perform ambient RNA correction. Contamination from ambient mRNA (e.g., from lysed cells) should be addressed with CellBender [@Fleming2023] or SoupX [@Young2020] before running AnnQC.

3. **Memory scaling.** Scrublet converts sparse count matrices to dense for PCA-based doublet simulation. For very large datasets (>50,000 cells × >10,000 genes), this may require >4 GB of RAM. The `--no-doublet-detection` flag allows users to skip this step when resources are constrained.

4. **Cluster impact analysis.** The cluster impact analysis in `annqc sensitivity` uses an exploratory clustering (PCA → neighbors → Leiden) computed from the current QC configuration. These clusters are diagnostic only and should not be interpreted as ground-truth cell types. Users with established cell type labels should supply them via `--cluster-labels`.

# Acknowledgements

The author thanks the Scanpy and AnnData development teams for the foundational data structures and QC utilities on which AnnQC is built.

# References

[@Wolf2018]: Wolf, F.A., Angerer, P., & Theis, F.J. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*, 19, 15.

[@Zheng2017]: Zheng, G.X.Y., et al. (2017). Massively parallel digital transcriptional profiling of single cells. *Nature Communications*, 8, 14049.

[@Wolock2019]: Wolock, S.L., Lopez, R., & Klein, A.M. (2019). Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data. *Cell Systems*, 8(4), 281–291.

[@Fleming2023]: Fleming, S.J., et al. (2023). Unsupervised removal of systematic background noise from droplet-based single-cell experiments using CellBender. *Nature Methods*, 20, 1323–1335.

[@Young2020]: Young, M.D., & Behjati, S. (2020). SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data. *GigaScience*, 9(12), giaa151.
