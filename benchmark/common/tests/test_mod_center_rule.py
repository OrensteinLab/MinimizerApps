"""The shared mod-minimizer reference (density.mod_center_selections): the center
tie rule is canonical (a sequence and its reverse complement select mirrored
k-mers) and forward (selections never move backwards), with t == k (mod w)."""
from __future__ import annotations

import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from density import (  # noqa: E402
    canonical_mmer,
    is_forward,
    mod_center_selections,
    mod_sampling_density,
    rolling_values,
    xor_hash,
)

SEED = 0x2545F4914F6CDD1D


def selections(seq: list[int], w: int, k: int, t: int) -> list[int]:
    L, T = w + k - 1, w + k - t
    ranks = [xor_hash(canonical_mmer(v, t), SEED) for v in rolling_values(seq, t)]

    def le_rc(s: int) -> bool:
        win = seq[s:s + L]
        return win <= [3 - b for b in reversed(win)]

    return mod_center_selections(ranks, w, T, le_rc)


@pytest.mark.parametrize("w,k,t", [(10, 15, 5), (10, 19, 9), (5, 35, 5), (4, 12, 4), (3, 8, 2), (11, 21, 10)])
@pytest.mark.parametrize("alphabet", [4, 2])  # i.i.d. ACGT, or {A, C} only (many ties)
def test_canonical_and_forward(w, k, t, alphabet):
    assert t % w == k % w
    rng = random.Random(w * 1000 + k * 10 + t + alphabet)
    seq = [rng.randrange(alphabet) for _ in range(3000)]
    rc = [3 - b for b in reversed(seq)]
    fwd, rev = selections(seq, w, k, t), selections(rc, w, k, t)
    assert is_forward(fwd) and is_forward(rev)
    L, n = w + k - 1, len(seq)
    for s, pos in enumerate(fwd):
        win = seq[s:s + L]
        if win == [3 - b for b in reversed(win)]:
            continue  # its own reverse complement: both strands take the left offset
        mirror = (n - k) - rev[len(rev) - 1 - s]  # k-mer start in the forward sequence
        assert pos == mirror, (s, pos, mirror)


def test_unique_minimum_needs_no_comparison():
    calls = []
    ranks = list(range(100, 0, -1)) + list(range(100))  # one minimum per window
    out = mod_center_selections(ranks, 4, 8, lambda s: calls.append(s) or True)
    assert calls == [] and is_forward(out)


def test_density_counts_distinct_positions():
    rng = random.Random(7)
    seq = [rng.randrange(4) for _ in range(5000)]
    n_kmers, n_sel = mod_sampling_density(seq, 15, 10, 5, xor_hash, SEED)
    assert n_kmers == 5000 - 15 + 1
    assert n_sel == len(set(selections(seq, 10, 15, 5)))
