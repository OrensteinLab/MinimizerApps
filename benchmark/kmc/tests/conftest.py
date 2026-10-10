"""Make `scripts` (benchmark/kmc) importable however pytest is invoked, e.g.
`python -m pytest benchmark/kmc/tests benchmark/common` from the repo root."""

import sys
from pathlib import Path

_KMC_BENCH = str(Path(__file__).resolve().parents[1])
if _KMC_BENCH not in sys.path:
    sys.path.insert(0, _KMC_BENCH)
