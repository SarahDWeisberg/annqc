# AnnQC — Full Project Specification for Claude Code

## What We Are Building

AnnQC is a Python package that takes raw single-cell RNA sequencing data, runs
standardized quality control steps, and produces:
1. A cleaned AnnData (.h5ad) file ready for downstream analysis
2. A beautiful, shareable HTML report showing exactly what was done and why

The philosophy: **AnnQC standardizes QC decisions, records them, and makes them shareable.**
It is not inventing new algorithms. It is packaging best-practice QC into a reproducible,
CLI-friendly, AnnData-first workflow.

The report is the product. Everything else serves the report.

---

## Package Name
`annqc`

## Install command (target)
```bash
pip install annqc
```

## Primary CLI command
```bash
annqc run raw.h5ad \
  --output cleaned.h5ad \
  --report annqc_report.html \
  --config qc_config.yaml \
  --sample-key sample \
  --seed 42
```

---

## Tech Stack

- **Python >= 3.9**
- **scanpy** — QC metrics, filtering, normalization, plotting
- **scrublet** — doublet detection
- **anndata** — data format throughout
- **jinja2** — HTML report templating
- **plotly** — interactive plots in report
- **pyyaml** — config parsing
- **click** — CLI
- **pytest** — testing

No R. No rpy2. No GPU dependencies. Pure Python only for v0.1.

---

## Repository Structure

```
annqc/
├── annqc/
│   ├── __init__.py
│   ├── cli.py              # click CLI entrypoint
│   ├── pipeline.py         # main orchestrator - runs all steps in order
│   ├── qc.py               # QC metric calculations (wraps scanpy)
│   ├── filter.py           # filtering logic with threshold decisions
│   ├── doublets.py         # scrublet wrapper
│   ├── config.py           # YAML config loading and validation
│   ├── spec.py             # defines adata.uns["annqc"] schema
│   ├── report/
│   │   ├── __init__.py
│   │   ├── builder.py      # assembles report data
│   │   ├── plots.py        # all plotly figures
│   │   └── templates/
│   │       └── report.html # jinja2 HTML template
│   └── utils.py            # logging, helpers
├── tests/
│   ├── conftest.py         # shared fixtures, PBMC 3k download
│   ├── test_qc.py
│   ├── test_filter.py
│   ├── test_doublets.py
│   ├── test_config.py
│   ├── test_pipeline.py    # seeded integration test
│   └── snapshots/          # snapshot test outputs
│       └── pbmc3k_summary.json
├── examples/
│   └── qc_config.yaml      # example config with all options documented
├── pyproject.toml
├── README.md               # README with screenshot on line 1
└── .github/
    └── workflows/
        └── ci.yml          # GitHub Actions: test on push
```

---

## YAML Config Schema

This is the public interface. Define it carefully — every field must be documented.

```yaml
# qc_config.yaml — AnnQC configuration

# Mitochondrial gene filtering
mito:
  prefix: "MT-"          # gene name prefix for mito genes (use "mt-" for mouse)
  max_pct: 20            # maximum % mitochondrial reads (cells above this are removed)

# Ribosomal gene filtering (optional)
ribo:
  prefix: "RPS|RPL"      # regex prefix for ribo genes
  max_pct: 50            # maximum % ribosomal reads (null = disabled)

# Count/gene thresholds
cells:
  min_genes: 200         # minimum number of genes per cell
  max_genes: 6000        # maximum number of genes per cell (catches doublets)
  min_counts: 500        # minimum total UMI counts per cell
  max_counts: null       # maximum total UMI counts (null = disabled)

# Gene-level filtering
genes:
  min_cells: 3           # minimum number of cells a gene must appear in

# Doublet detection
doublets:
  method: "scrublet"     # currently only "scrublet" supported
  threshold: "auto"      # auto = scrublet's default threshold, or set a float e.g. 0.25
  simulate_doublet_ratio: 2.0

# Normalization (applied after filtering)
normalization:
  method: "log1p"        # "log1p" or "none"
  target_sum: 10000      # normalize each cell to this total count

# Report options
report:
  title: "Single-Cell QC Report"
  author: ""             # optional author name shown in report
```

---

## AnnData Spec (`adata.uns["annqc"]`)

This is written to the output AnnData. It is the reproducibility record.
Every field must be present in every run.

```python
adata.uns["annqc"] = {
    "version": "0.1.0",           # annqc version used
    "date": "2025-01-15T10:30:00",# ISO timestamp of run
    "seed": 42,                   # random seed used
    "config": { ... },            # full config dict that was used (copy of YAML)
    "input_file": "raw.h5ad",     # path to input file
    "thresholds": {               # actual thresholds applied
        "mito_max_pct": 20,
        "min_genes": 200,
        "max_genes": 6000,
        "min_counts": 500,
        "doublet_threshold": 0.25  # resolved value even if "auto" was set
    },
    "cell_counts": {
        "input": 10000,
        "after_mito_filter": 9800,
        "after_gene_filter": 9500,
        "after_count_filter": 9400,
        "after_doublet_filter": 9100,
        "output": 9100
    },
    "per_sample": {               # only present if --sample-key used
        "sample_A": { "input": 5000, "output": 4500, "status": "PASS" },
        "sample_B": { "input": 5000, "output": 4600, "status": "PASS" }
    },
    "warnings": [                 # list of any auto-generated warnings
        "Sample B: high median mito% (18.2) — approaching threshold"
    ],
    "status": "PASS"             # overall PASS or FAIL
}
```

Per-cell columns added to `adata.obs`:
```python
adata.obs["annqc_pass"]           # bool: True if cell passed all filters
adata.obs["annqc_doublet_score"]  # float: scrublet doublet score
adata.obs["annqc_is_doublet"]     # bool: True if flagged as doublet
adata.obs["annqc_filter_reason"]  # str: reason for removal, "" if kept
                                  # values: "mito", "min_genes", "max_genes",
                                  #         "min_counts", "max_counts", "doublet"
```

---

## Pipeline Logic (`pipeline.py`)

The orchestrator runs these steps in order. Each step reads from and writes to
the same AnnData object. Each step logs to `adata.uns["annqc"]`.

```
Step 1: Load input (.h5ad or 10x .h5 directory)
Step 2: Calculate QC metrics (scanpy.pp.calculate_qc_metrics)
         - annotate mito genes, ribo genes
         - compute pct_counts_mt, pct_counts_ribo, n_genes_by_counts, total_counts
Step 3: Mito filtering — flag cells exceeding max_pct mito
Step 4: Gene/count filtering — flag cells outside min/max gene and count thresholds
Step 5: Gene-level filtering — remove genes appearing in too few cells
Step 6: Doublet detection — run scrublet, add scores and calls to obs
Step 7: Apply all flags — actually remove flagged cells from AnnData
Step 8: Normalization — log1p normalize if requested
Step 9: Write adata.uns["annqc"] record
Step 10: Save cleaned .h5ad
Step 11: Build and save HTML report
```

Each step must:
- Be independently testable
- Update `adata.uns["annqc"]["cell_counts"]` with cells remaining
- Log to Python's standard `logging` module
- Never silently fail — raise clear errors with helpful messages

---

## HTML Report Specification

This is the most important part of the project. The report must be good enough
that a scientist would screenshot it and put it in a lab meeting presentation.

### Technology
- Single self-contained HTML file (no external dependencies)
- Plotly figures embedded as JSON (interactive)
- Jinja2 template rendering
- Inline CSS — clean, modern, readable. Think clinical dashboard, not academic paper.

### Report Sections (in order)

#### 1. Header
- AnnQC logo/name
- Report title (from config)
- Run date, annqc version, input file name

#### 2. Executive Summary (most important section — top of page)
Large, scannable cards showing:
- Total input cells → Total output cells → % retained
- Overall status: big green PASS or red FAIL badge
- If per-sample mode: number of samples PASS / FAIL

#### 3. Per-Sample QC Table (if sample_key provided)
A table with one row per sample:
```
| Sample | Input Cells | Output Cells | Median Genes | Mito % | Doublet % | Status |
```
- Status column: green PASS / red FAIL / yellow WARN badge
- Rows sortable by clicking header
- Flagged samples highlighted in yellow

#### 4. QC Metric Distributions
For each metric (n_genes_by_counts, total_counts, pct_counts_mt, pct_counts_ribo):
- Violin plot or histogram
- Threshold line(s) drawn on the plot with a label showing the value
- "Before filtering" shown in gray, "After filtering" shown in color
- Per-sample traces if sample_key is provided

#### 5. Filtering Decisions
A step-by-step table:
```
| Step               | Cells Before | Cells Removed | Cells After | % Removed |
| Mito filtering     | 10,000       | 200           | 9,800       | 2.0%      |
| Gene filtering     | 9,800        | 300           | 9,500       | 3.1%      |
| Count filtering    | 9,500        | 100           | 9,400       | 1.1%      |
| Doublet removal    | 9,400        | 300           | 9,100       | 3.2%      |
```

#### 6. Doublet Detection Summary
- Distribution plot of doublet scores with threshold line
- Total cells flagged as doublets
- Doublet rate per sample (if applicable)

#### 7. ⚠️ Warnings (only shown if any exist)
Auto-generated warnings in yellow boxes, e.g.:
- "Sample B: median mito% is 18.2 — approaching threshold of 20"
- "High doublet rate detected (8.2%) — expected rate is 3-5% for 10x"
- "Very few cells remain after filtering (n=450) — check thresholds"

#### 8. Reproducibility Footer
Collapsible section containing:
- Full YAML config used for this run (copy-pasteable)
- Command that was run
- Software versions (annqc, scanpy, scrublet, python)
- A note: "To reproduce this run exactly, save this config and use --seed 42"

---

## CLI Specification (`cli.py`)

```bash
# Main command
annqc run INPUT \
  --output OUTPUT.h5ad \
  --report report.html \
  --config config.yaml \
  --sample-key COLUMN_NAME \
  --seed 42 \
  --verbose

# Validate a config file without running
annqc validate-config config.yaml

# Print default config to stdout (so user can customize it)
annqc init-config > my_config.yaml

# Show version
annqc --version
```

Arguments:
- `INPUT` — positional, required. Path to .h5ad file OR path to 10x directory (auto-detected)
- `--output` — path for cleaned output .h5ad. Default: `annqc_cleaned.h5ad`
- `--report` — path for HTML report. Default: `annqc_report.html`
- `--config` — path to YAML config. If not provided, use built-in defaults.
- `--sample-key` — column name in adata.obs for per-sample mode. Optional.
- `--seed` — integer random seed for reproducibility. Default: 0
- `--verbose` — show detailed logging output

---

## Testing Specification

### Test data
Use the PBMC 3k dataset from scanpy's built-in datasets:
```python
import scanpy as sc
adata = sc.datasets.pbmc3k()
```
This is ~2.7k cells, small enough to run fast in CI.

### Three layers of tests required:

**Layer 1: Unit tests**
- `test_qc.py` — test that QC metrics are calculated correctly
- `test_filter.py` — test that filtering removes exactly the right cells
- `test_doublets.py` — test that scrublet runs and adds correct obs columns
- `test_config.py` — test YAML loading, validation, defaults, error cases

**Layer 2: Seeded integration test**
```python
def test_full_pipeline_deterministic():
    adata = sc.datasets.pbmc3k()
    result = annqc.run(adata, config=default_config, seed=42)
    assert result.n_obs == EXPECTED_CELL_COUNT  # hardcode this after first run
    assert "annqc" in result.uns
    assert result.uns["annqc"]["status"] in ["PASS", "FAIL"]
    assert all(col in result.obs for col in [
        "annqc_pass", "annqc_doublet_score", "annqc_is_doublet", "annqc_filter_reason"
    ])
```

**Layer 3: Snapshot test**
```python
def test_summary_snapshot():
    # Run pipeline, extract summary dict, compare to stored JSON
    result = annqc.run(adata, config=default_config, seed=42)
    summary = {
        "n_cells_output": result.n_obs,
        "thresholds": result.uns["annqc"]["thresholds"],
        "cell_counts": result.uns["annqc"]["cell_counts"],
        "status": result.uns["annqc"]["status"]
    }
    with open("tests/snapshots/pbmc3k_summary.json") as f:
        expected = json.load(f)
    assert summary == expected
```

---

## Python API (in addition to CLI)

The package must also be importable and usable as a library:

```python
import annqc

# Run full pipeline
cleaned_adata = annqc.run(
    adata,
    config="config.yaml",    # or a dict
    sample_key="sample",
    seed=42
)

# Run individual steps
adata = annqc.calculate_qc_metrics(adata, mito_prefix="MT-")
adata = annqc.filter_cells(adata, min_genes=200, max_mito_pct=20)
adata = annqc.detect_doublets(adata, threshold="auto", seed=42)

# Generate report from an already-processed AnnData
annqc.generate_report(adata, output="report.html")
```

---

## Error Handling Rules

- If input file doesn't exist: clear error message with the path
- If a required config field is missing: tell the user exactly which field
- If scrublet fails: catch, warn, continue without doublet scores (don't crash)
- If output directory doesn't exist: create it automatically
- If adata has no cells after filtering: raise an error with a helpful message
  ("All cells were removed — check your thresholds")
- Never print raw Python tracebacks to the user unless --verbose is set

---

## README Requirements

The README must have:
1. A screenshot of the HTML report — FIRST THING after the title
2. One-line install: `pip install annqc`
3. Quickstart: 3 lines of code to run it
4. What it does (2 sentences max)
5. Full CLI reference
6. Full YAML config reference
7. How to use as a Python library
8. How to cite (placeholder for JOSS)

---

## pyproject.toml Requirements

```toml
[project]
name = "annqc"
version = "0.1.0"
description = "Reproducible QC reports for AnnData single-cell RNA-seq workflows"
requires-python = ">=3.9"
dependencies = [
    "scanpy>=1.9",
    "anndata>=0.9",
    "scrublet>=0.2",
    "click>=8.0",
    "pyyaml>=6.0",
    "jinja2>=3.0",
    "plotly>=5.0",
]

[project.scripts]
annqc = "annqc.cli:main"

[project.optional-dependencies]
dev = ["pytest", "pytest-cov", "black", "ruff"]
```

---

## GitHub Actions CI

On every push and pull request:
1. Run `pytest` with coverage
2. Check that coverage is above 80%
3. Run `ruff` linting
4. Test on Python 3.9, 3.10, 3.11

---

## What NOT to Build in v0.1

Do not include any of these — they are for later versions:
- SoupX (R dependency)
- CellBender (GPU dependency)
- scran normalization (R dependency)
- scDblFinder (R dependency)
- rpy2 of any kind
- Nextflow integration
- Multi-modal data support
- Cell type annotation

If someone asks for these during development, add a placeholder error:
`"CellBender support is planned for v0.2. See GitHub issues."`

---

## Definition of Done for v0.1

The package is complete when:
- [ ] `pip install annqc` works
- [ ] `annqc run pbmc3k.h5ad --output clean.h5ad --report report.html` works end to end
- [ ] The HTML report has all 8 sections listed above
- [ ] All unit tests pass
- [ ] Seeded integration test is deterministic (same output every run)
- [ ] Snapshot test is committed
- [ ] README has a report screenshot
- [ ] GitHub Actions CI passes
- [ ] `annqc init-config` prints a valid default config
- [ ] `adata.uns["annqc"]` contains all required fields
