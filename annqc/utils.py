"""Logging helpers and utility functions for AnnQC."""

import importlib.metadata
import logging
import os
import sys


def get_logger(name: str) -> logging.Logger:
    """Return a configured logger with the given name."""
    return logging.getLogger(name)


def setup_logging(verbose: bool = False) -> None:
    """Configure root logging. INFO level by default, DEBUG if verbose=True."""
    root = logging.getLogger()
    if not root.handlers:
        handler = logging.StreamHandler(sys.stderr)
        fmt = "%(asctime)s [%(name)s] %(levelname)s: %(message)s"
        handler.setFormatter(logging.Formatter(fmt))
        root.addHandler(handler)
    root.setLevel(logging.DEBUG if verbose else logging.INFO)


def ensure_dir(path: str) -> None:
    """Create directory and parents if they don't exist."""
    os.makedirs(path, exist_ok=True)


def format_number(n: int) -> str:
    """Format integer with thousands separator, e.g. 10000 -> '10,000'."""
    return f"{int(n):,}"


def get_software_versions() -> dict:
    """Return dict of package versions: annqc, scanpy, scrublet, python."""
    versions = {"python": sys.version.split()[0]}
    for pkg in ("annqc", "scanpy", "scrublet", "anndata", "plotly", "jinja2"):
        try:
            versions[pkg] = importlib.metadata.version(pkg)
        except importlib.metadata.PackageNotFoundError:
            versions[pkg] = "unknown"
    return versions
