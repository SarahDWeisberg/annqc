"""Command-line interface for AnnQC."""

import os
import sys
import logging

import click

import annqc
from annqc import __version__
from annqc.utils import setup_logging
from annqc.config import load_config, config_to_yaml, get_default_config

logger = logging.getLogger(__name__)


@click.group()
@click.version_option(version=__version__, prog_name="annqc")
def main():
    """AnnQC: reproducible QC reports for single-cell RNA-seq data."""


@main.command("run")
@click.argument("input", metavar="INPUT")
@click.option(
    "--output",
    default="annqc_cleaned.h5ad",
    show_default=True,
    help="Path for the cleaned output .h5ad file.",
)
@click.option(
    "--report",
    default="annqc_report.html",
    show_default=True,
    help="Path for the HTML QC report.",
)
@click.option(
    "--config",
    default=None,
    help="Path to a YAML config file. Omit to use built-in defaults.",
)
@click.option(
    "--sample-key",
    default=None,
    help="obs column name for per-sample statistics mode.",
)
@click.option(
    "--seed",
    default=0,
    show_default=True,
    type=int,
    help="Random seed for reproducibility.",
)
@click.option(
    "--verbose",
    is_flag=True,
    default=False,
    help="Enable debug logging and show full tracebacks on error.",
)
@click.option(
    "--auto-thresholds",
    is_flag=True,
    default=False,
    help="Auto-calculate QC thresholds using MAD-based statistics from the data.",
)
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help="Preview what would be filtered without writing output. No .h5ad is written.",
)
@click.option(
    "--no-doublet-detection",
    is_flag=True,
    default=False,
    help="Skip Scrublet doublet detection. All cells are kept as non-doublets.",
)
def run_cmd(input, output, report, config, sample_key, seed, verbose, auto_thresholds, dry_run, no_doublet_detection):
    """Run the AnnQC pipeline on INPUT.

    INPUT may be a path to a .h5ad file, a 10x HDF5 (.h5) file, or a
    10x Genomics MTX directory.
    """
    setup_logging(verbose)

    try:
        adata = annqc.run(
            input,
            config=config,
            sample_key=sample_key,
            seed=seed,
            output=output,
            report_path=report,
            input_file=input,
            dry_run=dry_run,
            auto_thresholds=auto_thresholds,
            no_doublet_detection=no_doublet_detection,
        )
    except Exception as exc:
        if verbose:
            raise
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    cc = adata.uns["annqc"]["cell_counts"]
    n_in = cc["input"]
    n_out = cc["output"]
    pct_retained = 100.0 * n_out / n_in if n_in > 0 else 0.0

    if dry_run:
        click.echo("AnnQC DRY RUN complete ✓ — no cells removed, no output file written")
        click.echo(f"  Input:    {n_in:,} cells")
        click.echo(f"  Would remove: {n_in - n_out:,} cells ({(1-n_out/n_in)*100:.2f}%)")
        click.echo(f"  Would retain: {n_out:,} cells ({pct_retained:.2f}%)")
        click.echo(f"  Preview report: {report}")
    else:
        click.echo("AnnQC complete ✓")
        click.echo(f"  Input:    {n_in:,} cells")
        click.echo(f"  Output:   {n_out:,} cells ({pct_retained:.2f}% retained)")
        click.echo(f"  Report:   {report}")

    if no_doublet_detection:
        click.echo("\n⚠️  WARNING: Doublet detection skipped.")
        click.echo("    QC is incomplete with respect to multiplets.")
        click.echo("    Do not use --no-doublet-detection for publication-quality analysis.")

    if auto_thresholds:
        try:
            thr = adata.uns["annqc"]["thresholds"]
            click.echo("\nAuto-thresholds used (MAD-based):")
            try:
                click.echo(f"  mito max pct:  {thr.get('mito_max_pct'):.2f}%")
            except (TypeError, ValueError):
                click.echo(f"  mito max pct:  {thr.get('mito_max_pct')}")
            try:
                click.echo(f"  min genes:     {thr.get('min_genes'):.0f}")
            except (TypeError, ValueError):
                click.echo(f"  min genes:     {thr.get('min_genes')}")
            try:
                click.echo(f"  max genes:     {thr.get('max_genes'):.0f}")
            except (TypeError, ValueError):
                click.echo(f"  max genes:     {thr.get('max_genes')}")
            try:
                click.echo(f"  min counts:    {thr.get('min_counts'):.0f}")
            except (TypeError, ValueError):
                click.echo(f"  min counts:    {thr.get('min_counts')}")
        except (KeyError, TypeError):
            pass


@main.command("validate-config")
@click.argument("config_path", metavar="CONFIG_PATH")
def validate_config_cmd(config_path):
    """Validate a YAML config file and report any errors.

    CONFIG_PATH is the path to the YAML config file to validate.
    """
    try:
        load_config(config_path)
        click.echo("✓ Config is valid.")
    except Exception as exc:
        click.echo(f"✗ Config error: {exc}", err=True)
        sys.exit(1)


@main.command("init-config")
def init_config_cmd():
    """Print the default AnnQC config YAML to stdout with inline documentation.

    Redirect to a file to create a starter config:

        annqc init-config > my_config.yaml
    """
    default_yaml = """\
# AnnQC configuration file
# All fields are shown. Unspecified fields fall back to built-in defaults.
# Pass this file to `annqc run --config <path>`.

# ---------------------------------------------------------------------------
# Mitochondrial gene filtering
# ---------------------------------------------------------------------------
mito:
  prefix: "MT-"     # Gene name prefix used to identify mitochondrial genes.
                    # Use "MT-" for human datasets, "mt-" for mouse datasets.
  max_pct: 20       # Maximum % mitochondrial reads allowed per cell (0-100).
                    # Cells above this threshold are removed.
                    # Typical range: 10-25 for human 10x data.

# ---------------------------------------------------------------------------
# Ribosomal gene filtering
# ---------------------------------------------------------------------------
ribo:
  prefix: "RPS|RPL" # Regex prefix pattern to identify ribosomal protein genes.
                    # Matches RPS* (small subunit) and RPL* (large subunit).
  max_pct: 50       # Maximum % ribosomal reads allowed per cell (0-100).
                    # High ribo% can indicate stressed or dying cells.
                    # Set to 100 to effectively disable this filter.

# ---------------------------------------------------------------------------
# Cell-level count/gene filters
# ---------------------------------------------------------------------------
cells:
  min_genes: 200    # Minimum number of genes detected per cell.
                    # Cells with fewer unique genes are removed (likely empty
                    # droplets or very low-quality cells). Typical: 200-500.
  max_genes: 6000   # Maximum number of genes detected per cell.
                    # Cells above this may be multiplets. Typical: 4000-8000.
                    # Set to null to disable the upper gene-count filter.
  min_counts: 500   # Minimum total UMI counts per cell.
                    # Cells below this are considered low-quality. Typical: 500-1000.
  max_counts: null  # Maximum total UMI counts per cell.
                    # Set to an integer to cap UMI counts (removes potential doublets).
                    # null = no upper limit applied.

# ---------------------------------------------------------------------------
# Gene-level filters
# ---------------------------------------------------------------------------
genes:
  min_cells: 3      # Minimum number of cells in which a gene must be detected.
                    # Genes detected in fewer cells are dropped from the matrix.
                    # Typical: 3-10.

# ---------------------------------------------------------------------------
# Doublet detection
# ---------------------------------------------------------------------------
doublets:
  method: "scrublet" # Doublet detection algorithm. Currently only "scrublet"
                     # is supported.
  threshold: "auto"  # Doublet score threshold. "auto" lets Scrublet choose
                     # automatically (recommended). Alternatively, set to a
                     # float between 0.0 and 1.0 to override (e.g. 0.25).
  simulate_doublet_ratio: 2.0
                     # Ratio of simulated doublets to real cells used by
                     # Scrublet. Higher values improve threshold estimation
                     # but increase runtime. Typical: 1.5-3.0.

# ---------------------------------------------------------------------------
# Normalization
# ---------------------------------------------------------------------------
normalization:
  method: "log1p"   # Normalization strategy applied after filtering.
                    # "log1p": normalize each cell to target_sum counts, then
                    #           apply log(x + 1) — the standard Scanpy workflow.
                    # "none":  skip normalization entirely (useful if you plan
                    #           to normalize downstream).
  target_sum: 10000 # Total counts per cell after normalization (used only
                    # when method = "log1p"). Standard value: 10000 (10k).

# ---------------------------------------------------------------------------
# Report metadata
# ---------------------------------------------------------------------------
report:
  title: "Single-Cell QC Report"
                    # Title displayed at the top of the HTML report.
  author: ""        # Optional author name shown in the report header.
                    # Leave blank to omit.
"""
    click.echo(default_yaml, nl=False)


@main.command("suggest")
@click.argument("input", metavar="INPUT")
@click.option(
    "--verbose",
    is_flag=True,
    default=False,
    help="Show debug output.",
)
def suggest_cmd(input, verbose):
    """Calculate and print suggested QC thresholds at three confidence levels.

    Does NOT run the full pipeline. Use this to explore thresholds before committing.

    INPUT may be a path to a .h5ad file, a 10x HDF5 (.h5) file, or a
    10x Genomics MTX directory.
    """
    setup_logging(verbose)

    try:
        import scanpy as sc
        from annqc.qc import calculate_qc_metrics
        from annqc.thresholds import suggest_thresholds
        from annqc.config import get_default_config

        # Load data
        if not os.path.exists(input):
            raise FileNotFoundError(f"Input not found: {input}")
        if os.path.isdir(input):
            adata = sc.read_10x_mtx(input, var_names="gene_symbols", cache=True)
        elif input.endswith(".h5"):
            adata = sc.read_10x_h5(input)
        else:
            adata = sc.read_h5ad(input)

        cfg = get_default_config()
        adata = calculate_qc_metrics(
            adata,
            mito_prefix=cfg["mito"]["prefix"],
            ribo_prefix=cfg["ribo"]["prefix"],
        )
        suggestions = suggest_thresholds(adata)

    except Exception as exc:
        if verbose:
            raise
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    click.echo(f"\nSuggested QC thresholds for: {input}")
    click.echo(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes\n")

    for level in ("strict", "standard", "permissive"):
        thr = suggestions[level]
        click.echo(f"  [{level.upper()}]")
        for k, v in thr.items():
            if v is not None:
                click.echo(f"    {k}: {v:.2f}")
        click.echo()

    click.echo("Use with: annqc run INPUT --auto-thresholds")
    click.echo("Or set manually in a config file: annqc init-config > my_config.yaml")


@main.command("compare")
@click.argument("file_a", metavar="FILE_A")
@click.argument("file_b", metavar="FILE_B")
@click.option(
    "--report",
    default="annqc_comparison.html",
    show_default=True,
    help="Output path for the comparison HTML report.",
)
@click.option("--verbose", is_flag=True, default=False)
def compare_cmd(file_a, file_b, report, verbose):
    """Compare two AnnQC-processed datasets and generate a comparison report.

    Both FILE_A and FILE_B must be .h5ad files processed by AnnQC
    (i.e., must have adata.uns['annqc'] populated).
    """
    setup_logging(verbose)

    try:
        import scanpy as sc
        from annqc.report.builder import build_comparison_report

        adata_a = sc.read_h5ad(file_a)
        adata_b = sc.read_h5ad(file_b)

        if "annqc" not in adata_a.uns:
            raise ValueError(f"{file_a} does not have adata.uns['annqc'] — was it processed by AnnQC?")
        if "annqc" not in adata_b.uns:
            raise ValueError(f"{file_b} does not have adata.uns['annqc'] — was it processed by AnnQC?")

        build_comparison_report(
            adata_a, adata_b,
            label_a=os.path.basename(file_a),
            label_b=os.path.basename(file_b),
            output_path=report,
        )

    except Exception as exc:
        if verbose:
            raise
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    click.echo(f"Comparison report written: {report}")


@main.command("methods")
@click.argument("input", metavar="INPUT")
@click.option(
    "--short",
    is_flag=True,
    default=False,
    help="Print the short (one-sentence) version instead of the full paragraph.",
)
def methods_cmd(input, short):
    """Print a ready-to-paste academic methods paragraph for INPUT.

    INPUT must be a .h5ad file previously processed by AnnQC
    (i.e., must have adata.uns['annqc'] populated).

    Example:
        annqc methods cleaned.h5ad
        annqc methods cleaned.h5ad --short
    """
    try:
        import scanpy as sc
        from annqc.methods_text import generate_methods_text, generate_methods_short

        if not os.path.exists(input):
            raise FileNotFoundError(f"File not found: {input}")

        adata = sc.read_h5ad(input)

        if "annqc" not in adata.uns:
            raise ValueError(
                f"{input} does not have adata.uns['annqc']. "
                "Was this file processed by AnnQC?"
            )

        if short:
            text = generate_methods_short(adata)
            click.echo("\n--- Short version (for abstracts) ---")
        else:
            text = generate_methods_text(adata)
            click.echo("\n--- Full methods paragraph (paste into your paper) ---")

        click.echo(text)
        click.echo()

    except Exception as exc:
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)


@main.command("sensitivity")
@click.argument("input", metavar="INPUT")
@click.option("--output", default=None, help="Path to save sensitivity results as JSON.")
@click.option(
    "--report",
    default="annqc_sensitivity.html",
    show_default=True,
    help="Path for the HTML sensitivity report.",
)
@click.option("--sample-key", default=None, help="obs column name for per-sample breakdown.")
@click.option("--seed", default=0, show_default=True, type=int, help="Random seed.")
@click.option(
    "--profiles",
    is_flag=True,
    default=False,
    help="Run strict/standard/permissive MAD profile comparison in addition to per-metric sweeps.",
)
@click.option(
    "--cluster-labels",
    default=None,
    metavar="LABELS_CSV",
    help="CSV file with columns cell_barcode,cluster. If provided, use these labels instead of running Leiden clustering.",
)
@click.option("--verbose", is_flag=True, default=False)
def sensitivity_cmd(input, output, report, sample_key, seed, profiles, cluster_labels, verbose):
    """Run threshold sensitivity analysis on INPUT.

    Tests a range of threshold values for each QC metric and shows
    how cell retention and cluster composition change.

    INPUT may be a path to a .h5ad file, a 10x HDF5 (.h5) file, or a
    10x Genomics MTX directory.
    """
    setup_logging(verbose)
    try:
        from annqc.sensitivity import run_sensitivity_analysis
        results = run_sensitivity_analysis(
            input,
            sample_key=sample_key,
            seed=seed,
            output_path=output,
            report_path=report,
            profiles=profiles,
            cluster_labels_path=cluster_labels,
        )
    except Exception as exc:
        if verbose:
            raise
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    n_input = results.get("n_input", 0)
    click.echo(f"\nSensitivity analysis complete — {n_input:,} input cells")
    for metric in ("mito", "min_genes", "max_genes", "min_counts"):
        m = results.get(metric, {})
        rec = m.get("mad_suggested")
        if rec is not None:
            click.echo(f"  {metric}: MAD-suggested = {rec:.3g}  ({m.get('mad_suggested_reason', '')})")

    if profiles and results.get("profiles"):
        click.echo("\nProfile comparison:")
        click.echo(f"  {'Profile':<12} {'Mito':>8} {'Min genes':>10} {'Max genes':>10} {'Retained':>10} {'%':>7}")
        click.echo(f"  {'-'*12} {'-'*8} {'-'*10} {'-'*10} {'-'*10} {'-'*7}")
        for row in results["profiles"]:
            if "error" not in row:
                click.echo(
                    f"  {row['profile']:<12} "
                    f"{row.get('mito_max_pct', '?'):>7.2f}% "
                    f"{row.get('min_genes', '?'):>10} "
                    f"{row.get('max_genes', '?'):>10} "
                    f"{row.get('n_output', '?'):>10,} "
                    f"{row.get('pct_retained', 0):>6.1f}%"
                )

    if profiles and results.get("profiles"):
        standard_row = next(
            (r for r in results["profiles"] if r.get("profile") == "standard" and "error" not in r),
            None,
        )
        click.echo(
            "\nNote: These profiles are exploratory. Standard 10x PBMC datasets typically "
            "retain 90-99% of cells with well-tuned thresholds."
        )
        if standard_row and standard_row.get("pct_retained", 100) < 85:
            click.echo(
                f"  Standard profile retained {standard_row['pct_retained']:.1f}% — "
                "review your data quality."
            )
        else:
            click.echo("  If standard profile retains <85%, review your data quality.")

    if report:
        click.echo(f"\nReport: {report}")
    if output:
        click.echo(f"JSON:   {output}")
