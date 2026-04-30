"""AnnQC: single-cell RNA-seq quality control pipeline."""

import sys

# Stub numba if it's not installed (llvmlite build failure on some platforms).
# scanpy imports numba at module load time; this must happen before scanpy loads.
if "numba" not in sys.modules:
    try:
        import numba  # noqa: F401 — real numba available, nothing to do
    except ImportError:
        from unittest.mock import MagicMock

        _stub = MagicMock()
        _stub.prange = range
        _stub.njit = lambda f=None, **kw: (f if f else lambda g: g)
        _stub.get_num_threads.return_value = 1
        sys.modules["numba"] = _stub
        sys.modules["numba.core"] = MagicMock()

        # scanpy's normalize_total uses a @njit function that only works under
        # real numba JIT; patch it with a plain-Python equivalent.
        import numpy as _np
        import scipy.sparse as _sp

        def _normalize_total_compat(adata, target_sum=None, **kwargs):
            X = adata.X
            if _sp.issparse(X):
                counts = _np.asarray(X.sum(axis=1)).flatten()
            else:
                counts = _np.asarray(X).sum(axis=1).flatten()
            if target_sum is None:
                target_sum = float(_np.median(counts[counts > 0]))
            scale = _np.where(counts > 0, target_sum / counts, 1.0)
            if _sp.issparse(X):
                adata.X = _sp.diags(scale) @ X
            else:
                adata.X = X * scale[:, _np.newaxis]

        import scanpy as _sc
        _sc.pp.normalize_total = _normalize_total_compat

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("annqc")
except PackageNotFoundError:
    __version__ = "0.0.0"

from annqc.pipeline import run

__all__ = ["run", "__version__"]
