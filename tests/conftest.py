import copy

import anndata as ad
import numpy as np
import pytest
import scanpy as sc

from annqc.config import DEFAULT_CONFIG


@pytest.fixture(scope="session")
def pbmc3k():
    """Download and return the PBMC 3k dataset from scanpy.

    Session-scoped so the dataset is downloaded only once per test run.
    """
    adata = sc.datasets.pbmc3k()
    return adata


@pytest.fixture
def small_adata():
    """Create a small synthetic AnnData (200 cells x 500 genes) for fast unit tests.

    The first 20 genes are given MT- prefixes to simulate mitochondrial genes.
    Cells 0-9 have their mitochondrial gene counts multiplied by 10, making them
    high-mito cells that should be flagged by the mito filter.
    """
    np.random.seed(42)
    n_cells, n_genes = 200, 500

    X = np.random.negative_binomial(5, 0.5, size=(n_cells, n_genes)).astype(float)

    var_names = [f"Gene{i}" for i in range(n_genes)]
    # Make first 20 genes mitochondrial
    var_names[:20] = [f"MT-Gene{i}" for i in range(20)]

    obs_names = [f"Cell{i}" for i in range(n_cells)]

    adata = ad.AnnData(
        X=X,
        obs={"cell_id": obs_names},
        var={"gene_id": var_names},
    )
    adata.obs_names = obs_names
    adata.var_names = var_names

    # Make cells 0-9 high-mito so they will be flagged
    adata.X[:10, :20] = adata.X[:10, :20] * 10

    return adata


@pytest.fixture
def small_adata_with_ribo():
    """Create a small synthetic AnnData with both MT- and RPS/RPL genes.

    200 cells x 500 genes.  Genes 0-19 are MT-, genes 20-39 are RPS/RPL.
    """
    np.random.seed(0)
    n_cells, n_genes = 200, 500

    X = np.random.negative_binomial(5, 0.5, size=(n_cells, n_genes)).astype(float)

    var_names = [f"Gene{i}" for i in range(n_genes)]
    var_names[:20] = [f"MT-Gene{i}" for i in range(20)]
    var_names[20:30] = [f"RPS{i}" for i in range(10)]
    var_names[30:40] = [f"RPL{i}" for i in range(10)]

    obs_names = [f"Cell{i}" for i in range(n_cells)]

    adata = ad.AnnData(
        X=X,
        obs={"cell_id": obs_names},
        var={"gene_id": var_names},
    )
    adata.obs_names = obs_names
    adata.var_names = var_names

    return adata


@pytest.fixture
def default_config():
    """Return a deep copy of DEFAULT_CONFIG so each test gets an independent instance."""
    return copy.deepcopy(DEFAULT_CONFIG)
