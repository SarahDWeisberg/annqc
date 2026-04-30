"""
annqc.methods_text
~~~~~~~~~~~~~~~~~~
Auto-generate academic methods paragraphs from a completed AnnQC run.
"""

import math


def generate_methods_text(adata) -> str:
    """Generate a full academic methods paragraph from adata.uns['annqc'].

    The returned string is ready to paste into a paper's Methods section.
    Works on both freshly-run and saved-then-loaded AnnData objects.
    """
    record = adata.uns.get("annqc", {})
    version = record.get("version", "unknown")
    cfg = record.get("config", {})
    thr = record.get("thresholds", {})
    cc = record.get("cell_counts", {})
    threshold_method = record.get("threshold_method", "manual")

    n_input = cc.get("input", 0)
    n_output = cc.get("output", 0)
    pct_retained = 100.0 * n_output / n_input if n_input > 0 else 0.0

    mito_prefix = cfg.get("mito", {}).get("prefix", "MT-")
    ribo_prefix = cfg.get("ribo", {}).get("prefix", "RPS|RPL")
    min_cells_gene = cfg.get("genes", {}).get("min_cells", 3)

    mito_max = thr.get("mito_max_pct")
    min_genes = thr.get("min_genes")
    max_genes = thr.get("max_genes")
    min_counts = thr.get("min_counts")

    doublet_status = record.get("doublet_status")
    dbl_thr = thr.get("doublet_threshold")
    dbl_thr_valid = (
        dbl_thr is not None
        and not (isinstance(dbl_thr, float) and math.isnan(float(dbl_thr)))
    )

    # Doublet stats
    n_after_count = cc.get("after_count_filter", n_input)
    n_after_doublet = cc.get("after_doublet_filter", n_after_count)
    n_doublets = max(0, n_after_count - n_after_doublet)
    doublet_rate = 100.0 * n_doublets / n_input if n_input > 0 else 0.0

    # Threshold method phrase
    if threshold_method == "auto_mad":
        thr_phrase = "MAD-based adaptive thresholds derived from the data distribution"
    else:
        thr_phrase = "manually specified thresholds"

    # Build threshold description
    thr_parts = []
    if mito_max is not None:
        thr_parts.append(f"mitochondrial content ≤ {mito_max:.3g}%")
    if min_genes is not None and max_genes is not None:
        thr_parts.append(
            f"gene count between {min_genes:.0f} and {max_genes:.0f} genes per cell"
        )
    elif min_genes is not None:
        thr_parts.append(f"minimum {min_genes:.0f} genes per cell")
    elif max_genes is not None:
        thr_parts.append(f"maximum {max_genes:.0f} genes per cell")
    if min_counts is not None and min_counts > 0:
        thr_parts.append(f"minimum UMI count of {min_counts:.0f}")

    thr_sentence = (
        f"Cells were filtered using {thr_phrase}"
        + (f": {', '.join(thr_parts)}" if thr_parts else "")
        + "."
    )

    # Gene filter sentence
    gene_filter_sentence = (
        f"Gene-level filtering removed genes expressed in fewer than "
        f"{min_cells_gene} cells."
    )

    # Doublet sentence
    if doublet_status == "SKIPPED":
        doublet_sentence = (
            "Doublet detection was not performed. Results may contain cell multiplets."
        )
    elif dbl_thr_valid:
        doublet_sentence = (
            f"Doublet detection was performed using Scrublet (Wolock et al., 2019) "
            f"with an automatically determined threshold of {float(dbl_thr):.3f}, "
            f"identifying and removing {n_doublets} doublets ({doublet_rate:.1f}% doublet rate)."
        )
    else:
        doublet_sentence = (
            "Doublet detection was not performed or was skipped."
        )

    # Summary sentence
    summary_sentence = (
        f"After quality control, following established best practices, "
        f"{n_output:,} of {n_input:,} cells were retained ({pct_retained:.1f}%). "
        "Sensitivity analysis was performed using annqc sensitivity to assess the impact "
        "of threshold choices on cell retention and cluster composition, using independent "
        "threshold sweeps and combined strict/standard/permissive MAD-based profiles. "
        "AnnQC does not perform ambient RNA correction; tools such as CellBender or SoupX "
        "should be applied upstream if ambient RNA removal is required. "
        "All QC parameters and software versions are recorded in "
        "adata.uns['annqc'] for full reproducibility."
    )

    paragraph = (
        f"Quality control was performed using AnnQC v{version} (Weisberg, 2026). "
        f"Mitochondrial genes were identified by the '{mito_prefix}' prefix and "
        f"ribosomal genes by the '{ribo_prefix}' prefix. "
        f"{thr_sentence} "
        f"{gene_filter_sentence} "
        f"{doublet_sentence} "
        f"{summary_sentence}"
    )

    return paragraph


def generate_methods_short(adata) -> str:
    """Generate a one-sentence abstract-length methods description."""
    record = adata.uns.get("annqc", {})
    version = record.get("version", "unknown")
    cc = record.get("cell_counts", {})
    threshold_method = record.get("threshold_method", "manual")

    n_input = cc.get("input", 0)
    n_output = cc.get("output", 0)
    pct_retained = 100.0 * n_output / n_input if n_input > 0 else 0.0

    thr_phrase = "MAD-based adaptive thresholds" if threshold_method == "auto_mad" else "standard thresholds"

    return (
        f"scRNA-seq data were quality controlled using AnnQC v{version} "
        f"with {thr_phrase}, retaining {n_output:,} of {n_input:,} cells "
        f"({pct_retained:.1f}%) after removal of low-quality cells and doublets "
        f"(ambient RNA correction not performed). "
        f"Threshold sensitivity was assessed with annqc sensitivity."
    )
