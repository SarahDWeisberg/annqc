"""Root conftest: stub numba + patch scanpy helpers that need JIT (build-free)."""
import sys
from unittest.mock import MagicMock

# --- 1. Stub numba before scanpy loads it ---
if "numba" not in sys.modules:
    _numba = MagicMock()
    _numba.prange = range

    def _njit(func=None, **kwargs):
        if func is not None:
            return func
        return lambda f: f

    _numba.njit = _njit
    sys.modules["numba"] = _numba
    sys.modules["numba.core"] = MagicMock()

# --- 2. Patch sc.pp.normalize_total to avoid numba-only code paths ---
import numpy as np
import scipy.sparse as _sp


def _normalize_total_plain(adata, target_sum=None, **kwargs):
    """Pure-Python normalize_total (no numba) for the test environment."""
    X = adata.X
    if _sp.issparse(X):
        counts = np.asarray(X.sum(axis=1)).flatten()
    else:
        counts = np.asarray(X).sum(axis=1).flatten()

    if target_sum is None:
        target_sum = np.median(counts[counts > 0])

    scale = np.where(counts > 0, target_sum / counts, 1.0)

    if _sp.issparse(X):
        from scipy.sparse import diags
        adata.X = diags(scale) @ X
    else:
        adata.X = X * scale[:, np.newaxis]


def pytest_configure(config):
    """Patch scanpy after it has been imported."""
    try:
        import scanpy as sc
        sc.pp.normalize_total = _normalize_total_plain
    except Exception:
        pass
