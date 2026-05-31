#!/bin/bash
# Run any Python script using the native ARM64 Python with the correct libexpat.
# Needed because the x86_64 Homebrew Python crashes (SIGILL) when loading tiledbsoma
# due to AVX instructions in the x86 binary not being supported under Rosetta 2.
#
# Usage:
#   bash scripts/run_arm64.sh scripts/build_references.py [args...]
#   bash scripts/run_arm64.sh scripts/build_references_local.py --tissue kidney ...
export DYLD_LIBRARY_PATH=/opt/homebrew/Cellar/expat/2.8.1/lib
exec /opt/homebrew/bin/python3.12 "$@"
