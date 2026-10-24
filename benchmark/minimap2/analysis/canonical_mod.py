#!/usr/bin/env python3
"""Canonical (strand-symmetric) mod-minimizer: brute-force check and density.

Schemes, per window W of L = w + k - 1 bases (w k-mers, T = w + k - t t-mers;
every rank is a strand-symmetric hash of the canonical t-mer / k-mer):

- minimizer:          leftmost minimal canonical k-mer (minimap2/kraken2's
                      plain minimizer, up to ties).
- canonical-minimizer: simd-minimizers' rule (Groot Koerkamp & Martayan, SEA
                      2025): leftmost minimum on a forward-preferred window,
                      rightmost on a reverse-preferred one.
- mod-minimizer:      the plain WABI 2024 scheme: k-mer at offset x_L mod w,
                      x_L = leftmost minimal t-mer (not strand-symmetric).
- canonical-mod:      the G/T-strand construction the minimap2 and kraken2
                      mod-minimizers used before the centre tie rule
                      (mod_minimizer_rule.h, 2026-10): forward-preferred
                      windows as mod-minimizer; reverse-preferred windows the
                      mirror rule, x' = (T-1) - x_R (x_R = rightmost minimal
                      t-mer), k-mer at offset (w-1) - (x' mod w).
Preferred strand: the one with more G/T bases (equivalently W has more G/T
than its reverse complement); a tie (even L) is broken by W <= rc(W)
lexicographically, so exactly one of W, rc(W) is forward-preferred unless W
is its own reverse complement.

Reported per (k, w, t, scheme): density factor (selected k-mer positions /
k-mers * (w+1)) on i.i.d. random DNA, and strand agreement: the fraction of
loci selected on a sequence whose mirror locus is also selected on its
reverse complement (1.0 = canonical).
"""
from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
from density import pcg_rxs_m_xs64  # noqa: E402

COMP = "TGCA"  # complement of ACGT by index


def revcomp(s: str) -> str:
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def canon_hash(s: str, seed: int) -> int:
    """Strand-symmetric rank of a substring: hash of its lexicographically
    smaller strand."""
    c = min(s, revcomp(s))
    v = 0
    for ch in c:
        v = (v << 2) | "ACGT".index(ch)
    return pcg_rxs_m_xs64(v ^ seed)


def forward_preferred(win: str) -> bool:
    gt = sum(1 for ch in win if ch in "GT")
    L = len(win)
    if 2 * gt != L:
        return 2 * gt > L
    return win <= revcomp(win)


def select(seq: str, k: int, w: int, t: int, scheme: str, seed: int) -> set[int]:
    """Start positions of the selected k-mers."""
    L = w + k - 1
    T = w + k - t
    kr = [canon_hash(seq[i:i + k], seed) for i in range(len(seq) - k + 1)]
    tr = [canon_hash(seq[i:i + t], seed) for i in range(len(seq) - t + 1)] if "mod" in scheme else None
    out = set()
    for s in range(len(seq) - L + 1):
        if scheme in ("minimizer", "canonical-minimizer"):
            vals = kr[s:s + w]
            m = min(vals)
            left = vals.index(m)
            right = w - 1 - vals[::-1].index(m)
            if scheme == "minimizer" or forward_preferred(seq[s:s + L]):
                p = left
            else:
                p = right
        else:
            vals = tr[s:s + T]
            m = min(vals)
            x_l = vals.index(m)
            x_r = T - 1 - vals[::-1].index(m)
            if scheme == "mod-minimizer" or forward_preferred(seq[s:s + L]):
                p = x_l % w
            else:
                p = (w - 1) - (((T - 1) - x_r) % w)
        out.add(s + p)
    return out


def evaluate(n: int, k: int, w: int, t: int, scheme: str, seed: int, rng: random.Random) -> tuple[float, float]:
    seq = "".join(rng.choice("ACGT") for _ in range(n))
    fwd = select(seq, k, w, t, scheme, seed)
    rc = select(revcomp(seq), k, w, t, scheme, seed)
    mirrored = {(n - k) - q for q in rc}
    density = len(fwd) / (n - k + 1) * (w + 1)
    agree = len(fwd & mirrored) / len(fwd)
    return density, agree


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=50_000, help="random sequence length")
    ap.add_argument("--reps", type=int, default=3)
    args = ap.parse_args(argv)
    configs = [(15, 10), (19, 10), (21, 11), (25, 10), (31, 11)]
    schemes = ["minimizer", "canonical-minimizer", "mod-minimizer", "canonical-mod"]
    print("k  w  t  scheme               density(mean)  strand_agreement(mean)")
    for k, w in configs:
        t = 4 + (k - 4) % w
        for scheme in schemes:
            ds, ag = [], []
            for rep in range(args.reps):
                rng = random.Random(1000 * k + 10 * w + rep)
                d, a = evaluate(args.n, k, w, t, scheme, 0x9E3779B97F4A7C15 + rep, rng)
                ds.append(d); ag.append(a)
            print(f"{k:<2} {w:<2} {t:<2} {scheme:20s} {sum(ds)/len(ds):.4f}         {sum(ag)/len(ag):.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
