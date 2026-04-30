# AnnQC

**Reproducible QC reports for AnnData single-cell RNA-seq workflows.**

Run standardized quality control on your scRNA-seq data and get a shareable HTML report that documents every decision made.

---

## Install

```bash
pip install annqc
```

## Quickstart

```python
import annqc
import scanpy as sc

adata = sc.read_h5ad("raw.h5ad")
cleaned = annqc.run(adata, output="cleaned.h5ad", report_path="qc_report.html")
```

Or from the command line:

```bash
annqc run raw.h5ad --output cleaned.h5ad --report qc_report.html --seed 42
```

---

## What it does

AnnQC applies best-practice QC steps to your AnnData object — mito/gene/count filtering, gene-level filtering, and doublet detection — then writes a self-contained HTML report showing exactly what was filtered and why. Every run is reproducible via a fixed seed and a full provenance record stored in `adata.uns["annqc"]`.

---

## CLI Reference

```
annqc run INPUT [OPTIONS]

  Run the full QC pipeline on a .h5ad file, 10x .h5 file, or 10x directory.

Arguments:
  INPUT             Path to input file or 10x directory [required]

Options:
  --output PATH     Output .h5ad path  [default: annqc_cleaned.h5ad]
  --report PATH     Output HTML report path  [default: annqc_report.html]
  --config PATH     YAML config file (omit to use built-in defaults)
  --sample-key COL  obs column name for per-sample statistics
  --seed INT        Random seed for reproducibility  [default: 0]
  --verbose         Show detailed log output
  --version         Show version and exit
  --help            Show this message and exit
```

```bash
# Validate a config without running
annqc validate-config my_config.yaml

# Print the default config to stdout
annqc init-config > my_config.yaml
```

---

## YAML Config Reference

```yaml
mito:
  prefix: "MT-"    # Gene prefix for mito genes (use "mt-" for mouse)
  max_pct: 20      # Max % mito reads per cell

ribo:
  prefix: "RPS|RPL"
  max_pct: 50

cells:
  min_genes: 200   # Min genes per cell
  max_genes: 6000  # Max genes per cell
  min_counts: 500  # Min total UMI counts
  max_counts: null # Max UMI counts (null = disabled)

genes:
  min_cells: 3     # Min cells a gene must appear in

doublets:
  method: "scrublet"
  threshold: "auto"
  simulate_doublet_ratio: 2.0

normalization:
  method: "log1p"   # "log1p" or "none"
  target_sum: 10000

report:
  title: "Single-Cell QC Report"
  author: ""
```

---

## Python API

```python
import annqc

# Full pipeline
cleaned = annqc.run(
    "raw.h5ad",
    config="my_config.yaml",   # or a dict, or None for defaults
    sample_key="sample",
    seed=42,
    output="cleaned.h5ad",
    report_path="qc_report.html",
)

# Inspect the provenance record
record = cleaned.uns["annqc"]
print(record["status"])          # "PASS" or "FAIL"
print(record["cell_counts"])     # input → output counts at each step
print(record["thresholds"])      # thresholds actually applied
```

---

## Reproducibility

Every run stores a full provenance record in `adata.uns["annqc"]`:

```python
{
    "version": "0.1.0",
    "date": "2025-01-15T10:30:00",
    "seed": 42,
    "config": { ... },           # full config used
    "input_file": "raw.h5ad",
    "thresholds": { ... },       # resolved threshold values
    "cell_counts": { ... },      # cells at each filtering step
    "warnings": [ ... ],
    "status": "PASS"
}
```

To reproduce a run: save the YAML config from the report's reproducibility section and rerun with the same `--seed`.

---

## Citation

If you use AnnQC in your work, please cite:

> *Citation pending JOSS submission.*

---

## License

MIT
