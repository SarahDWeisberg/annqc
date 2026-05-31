"""Tissue-aware QC presets for AnnQC."""

import copy

TISSUE_PRESETS = {
    "pbmc": {
        "description": "Peripheral blood mononuclear cells (standard 10x)",
        "mito": {"max_pct": 20},
        "cells": {"min_genes": 200, "max_genes": 6000, "min_counts": 500, "max_counts": None},
    },
    "tumor": {
        "description": "Tumor microenvironment (malignant + stromal + immune cells)",
        "mito": {"max_pct": 25},
        "cells": {"min_genes": 200, "max_genes": 8000, "min_counts": 500, "max_counts": None},
    },
    "brain_sn": {
        "description": "Single-nucleus RNA-seq from brain/CNS (nuclei, not whole cells)",
        "mito": {"max_pct": 5},
        "cells": {"min_genes": 200, "max_genes": 6000, "min_counts": 500, "max_counts": None},
    },
    "organoid": {
        "description": "Organoid cultures (elevated stress, mixed maturity)",
        "mito": {"max_pct": 25},
        "cells": {"min_genes": 200, "max_genes": 8000, "min_counts": 500, "max_counts": None},
    },
    "gut": {
        "description": "Intestinal / gut epithelium (high baseline mitochondrial content)",
        "mito": {"max_pct": 30},
        "cells": {"min_genes": 200, "max_genes": 6000, "min_counts": 500, "max_counts": None},
    },
    "heart": {
        "description": "Cardiac tissue (cardiomyocytes have very high mitochondrial content)",
        "mito": {"max_pct": 40},
        "cells": {"min_genes": 200, "max_genes": 8000, "min_counts": 500, "max_counts": None},
    },
    "kidney": {
        "description": "Kidney (elevated mito in proximal tubule cells)",
        "mito": {"max_pct": 30},
        "cells": {"min_genes": 200, "max_genes": 6000, "min_counts": 500, "max_counts": None},
    },
    "embryo": {
        "description": "Embryonic / developmental tissues (broad cell type heterogeneity)",
        "mito": {"max_pct": 20},
        "cells": {"min_genes": 200, "max_genes": 8000, "min_counts": 300, "max_counts": None},
    },
    "low_input": {
        "description": "Low-input / low-coverage experiments",
        "mito": {"max_pct": 20},
        "cells": {"min_genes": 100, "max_genes": None, "min_counts": 200, "max_counts": None},
    },
}


def get_preset(name: str) -> dict:
    """Return a config overlay dict for the named tissue preset.

    The returned dict has "mito" and "cells" keys only — suitable for
    deep-merging into a full config. Raises ValueError for unknown names.
    """
    if name not in TISSUE_PRESETS:
        valid = ", ".join(sorted(TISSUE_PRESETS))
        raise ValueError(f"Unknown tissue preset {name!r}. Valid presets: {valid}")
    preset = TISSUE_PRESETS[name]
    return {
        "mito": copy.deepcopy(preset["mito"]),
        "cells": copy.deepcopy(preset["cells"]),
    }


def list_presets() -> list:
    """Return a list of (name, description) tuples for all presets."""
    return [(name, info["description"]) for name, info in TISSUE_PRESETS.items()]
