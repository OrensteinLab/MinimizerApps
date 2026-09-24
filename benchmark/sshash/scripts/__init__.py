"""sshash benchmark Python package."""

import sys
from pathlib import Path

# Shared YAML-experiment engine + xlsx primitives (benchmark/common), also
# used by kraken2's and kmc's own scripts packages.
_shared = Path(__file__).resolve().parents[3] / "benchmark" / "common"
if str(_shared) not in sys.path:
    sys.path.insert(0, str(_shared))

from scripts.config import benchmark_dir, repo_root

__all__ = ["benchmark_dir", "repo_root"]
