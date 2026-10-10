#!/usr/bin/env python3
"""Generate the default "smoke" dataset: a small synthetic genome with an
engineered long, near-identical repeat region, plus simulated long reads with
known ground-truth origin.

Synthetic rather than a downloaded genome (per the approved plan): fast to
regenerate, no download dependency, and the repeat-vs-density comparison the
Winnowmap-style evaluation cares about only needs *some* long near-identical
repeat to exist, not biological realism. See README.md's "Default dataset"
section. A real genome (e.g. a human chromosome) can be swapped in later by
pointing MINIMIZER_ORDER_GENOME at it and skipping this generator.

Read simulation here is a lightweight, dependency-free stand-in for a real
long-read simulator (pbsim2 -- see environment.yml/setup.sh): it samples
contiguous windows from the genome and injects substitutions/indels at
ONT-like rates, recording true origin coordinates in truth.tsv. Swap
`simulate_reads_pbsim2()` in once the conda env's pbsim2 is available; this
keeps the harness runnable without it in the meantime.
"""
from __future__ import annotations

import argparse
import gzip
import random
import sys
from pathlib import Path

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def random_seq(rng: random.Random, n: int, gc: float = 0.5) -> str:
    # gc controls P(G or C) per base; independent bases (no k-mer structure
    # beyond what falls out of chance), which is exactly what "not trying to
    # be biologically realistic" means here.
    weights = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]  # A C G T
    return "".join(rng.choices(BASES, weights=weights, k=n))


def revcomp(s: str) -> str:
    return s.translate(COMPLEMENT)[::-1]


def make_genome(rng: random.Random, total_len: int, repeat_len: int) -> tuple[str, list[tuple[str, int, int]]]:
    """Returns (sequence, [(label, start, end), ...]) for the injected
    repeat's two copies (0-indexed, half-open), for the README/reports."""
    assert repeat_len * 2 < total_len
    unique_flanks = total_len - 2 * repeat_len
    left = unique_flanks // 3
    right = unique_flanks - left - (unique_flanks // 3)
    middle = unique_flanks // 3

    repeat_unit = random_seq(rng, repeat_len)
    a = random_seq(rng, left)
    b = random_seq(rng, middle)
    c = random_seq(rng, right)

    seq = a + repeat_unit + b + repeat_unit + c
    r1_start = len(a)
    r1_end = r1_start + repeat_len
    r2_start = r1_end + len(b)
    r2_end = r2_start + repeat_len
    assert len(seq) == total_len, (len(seq), total_len)
    return seq, [("repeat_copy_1", r1_start, r1_end), ("repeat_copy_2", r2_start, r2_end)]


def inject_errors(rng: random.Random, s: str, sub_rate: float, ins_rate: float, del_rate: float) -> str:
    out = []
    for c in s:
        r = rng.random()
        if r < del_rate:
            continue  # deletion: skip this base
        if rng.random() < sub_rate:
            c = rng.choice([b for b in BASES if b != c])
        out.append(c)
        if rng.random() < ins_rate:
            out.append(rng.choice(BASES))
    return "".join(out)


def simulate_reads(
    rng: random.Random,
    genome: str,
    n_reads: int,
    min_len: int,
    max_len: int,
    sub_rate: float,
    ins_rate: float,
    del_rate: float,
    repeat_spans: list[tuple[str, int, int]],
    repeat_frac: float,
) -> list[tuple[str, str, int, int, str]]:
    """Returns [(read_id, read_seq, true_start, true_end, true_strand), ...].
    A `repeat_frac` fraction of reads are forced to start inside one of the
    two repeat copies (uniformly split between them), matching the paper's
    own headline comparison (repeat-region mapping specifically)."""
    n = len(genome)
    reads = []
    n_repeat = int(round(n_reads * repeat_frac))
    for i in range(n_reads):
        read_len = rng.randint(min_len, max_len)
        if i < n_repeat and repeat_spans:
            label, rs, re_ = repeat_spans[i % len(repeat_spans)]
            lo = max(0, rs - read_len // 2)
            hi = min(n - read_len, re_ - read_len // 2)
            if hi <= lo:
                start = rs
            else:
                start = rng.randint(lo, hi)
        else:
            start = rng.randint(0, max(0, n - read_len))
        end = min(n, start + read_len)
        true_seq = genome[start:end]
        strand = rng.choice("+-")
        if strand == "-":
            true_seq = revcomp(true_seq)
        err_seq = inject_errors(rng, true_seq, sub_rate, ins_rate, del_rate)
        if len(err_seq) < 50:
            continue
        reads.append((f"sim_read_{i}", err_seq, start, end, strand))
    return reads


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for j in range(0, len(seq), 80):
                f.write(seq[j:j + 80] + "\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--genome-len", type=int, default=300_000)
    ap.add_argument("--repeat-len", type=int, default=20_000, help="length of the injected near-identical repeat")
    ap.add_argument("--n-reads", type=int, default=400)
    ap.add_argument("--min-read-len", type=int, default=2_000)
    ap.add_argument("--max-read-len", type=int, default=10_000)
    ap.add_argument("--sub-rate", type=float, default=0.06, help="ONT-ish per-base substitution rate")
    ap.add_argument("--ins-rate", type=float, default=0.03)
    ap.add_argument("--del-rate", type=float, default=0.03)
    ap.add_argument("--repeat-read-frac", type=float, default=0.3, help="fraction of reads forced to originate inside the repeat")
    args = ap.parse_args(argv)

    rng = random.Random(args.seed)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    genome, repeat_spans = make_genome(rng, args.genome_len, args.repeat_len)
    write_fasta(args.out_dir / "genome.fa", [("smoke_genome", genome)])

    reads = simulate_reads(
        rng, genome, args.n_reads, args.min_read_len, args.max_read_len,
        args.sub_rate, args.ins_rate, args.del_rate, repeat_spans, args.repeat_read_frac,
    )
    write_fasta(args.out_dir / "reads.fa", [(rid, seq) for rid, seq, *_ in reads])

    with (args.out_dir / "truth.tsv").open("w") as f:
        f.write("read_id\ttrue_contig\ttrue_start\ttrue_end\ttrue_strand\tin_repeat\n")
        for rid, _seq, start, end, strand in reads:
            # half-open overlap: also catches a read that fully spans a repeat
            in_repeat = any(start < re_ and rs < end for _l, rs, re_ in repeat_spans)
            f.write(f"{rid}\tsmoke_genome\t{start}\t{end}\t{strand}\t{1 if in_repeat else 0}\n")

    with (args.out_dir / "repeats.bed").open("w") as f:
        for label, rs, re_ in repeat_spans:
            f.write(f"smoke_genome\t{rs}\t{re_}\t{label}\n")

    print(
        f"wrote {args.out_dir}: genome={args.genome_len}bp (repeat={args.repeat_len}bp x2), "
        f"{len(reads)} reads ({args.repeat_read_frac:.0%} targeted at the repeat)",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
