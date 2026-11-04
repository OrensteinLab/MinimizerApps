"""Minimizer density helpers of the sshash benchmark: the density factor and the
GreedyMini table checks of scripts/schemes.py. The expected density is counted
by sshash's own code (scripts/expected_density.py, benchmark/common/expected_density/).

density_factor = density * (w + 1), where density = n_selected / n_kmers and
w = k - m + 1 is the window size (see orders/GreedyMini/w11_k15_canonical.*:
the "w11" in that filename *is* this w, for k=31/m=21). A uniformly random
total order has expected density_factor -> 2 as w grows (Roberts et al.
2004); GreedyMini's whole point is to push it below that.
"""

from __future__ import annotations

import struct
from pathlib import Path

# benchmark/common (on sys.path via scripts/__init__): the extension of a table to an
# l-mer needs l = gm_k + 2i, the same number of extra bases on both sides
from scheme_names import gm_center_fits as gm_center_valid  # noqa: F401

INF_RANK = 0xFFFF


def density_factor(n_selected: int, n_kmers: int, w: int) -> float | None:
    if n_kmers <= 0 or n_selected <= 0:
        return None
    return round(n_selected / n_kmers * (w + 1), 10)


def read_gm_ranks(path: Path) -> list[int]:
    """Read a classic .gm file ([u64 len][len x u64] raw ranks; entries equal
    to `len` are the "not in UHS" sentinel) and return the dense uint16-style
    rank table (INF_RANK for non-UHS), matching GreedyMiniTable::load_from_raw."""
    data = path.read_bytes()
    (n,) = struct.unpack_from("<Q", data, 0)
    raw = struct.unpack_from(f"<{n}Q", data, 8)
    return [INF_RANK if r >= n else r for r in raw]


def _gm_k_from_table(gm_ranks: list[int]) -> int:
    """Binary-alphabet GreedyMini tables have exactly 2**gm_k entries."""
    n = len(gm_ranks)
    gm_k = n.bit_length() - 1
    assert (1 << gm_k) == n, f"GreedyMini table size {n} is not a power of 2"
    return gm_k


def gm_table_mirror_symmetric(gm_ranks: list[int]) -> bool:
    gm_k = _gm_k_from_table(gm_ranks)
    return all(
        gm_ranks[key] == gm_ranks[int(format(key, f"0{gm_k}b")[::-1], 2)] for key in range(len(gm_ranks))
    )
