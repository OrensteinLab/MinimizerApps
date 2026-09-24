#!/usr/bin/env python3
"""Simulate reads for an existing dataset (genome.fa + repeats.bed, from
make_dataset.py) using a real ONT-model long-read simulator (pbsim2, see
environment.yml/setup.sh), replacing make_dataset.py's own lightweight
built-in simulator's reads.fa/truth.tsv with pbsim2's.

pbsim2 doesn't target specific loci, so repeat-region coverage isn't forced
the way make_dataset.py's own simulator does; with --depth high enough over
a genome whose repeat is ~7% of its length (20kb/300kb here), a comparable
number of reads still land there by chance -- reported in the output.

pbsim2's own <prefix>_0001.maf carries each read's true (start, end, strand)
directly (the `s ref <start> <alnlen> <strand> ...` line right before each
read's own `s <read_id> ...` line): parsed here into the same truth.tsv
shape make_dataset.py produces, so score_accuracy.py/lib/run.sh need no
changes to consume either.
"""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path


def load_repeats(path: Path) -> list[tuple[int, int]]:
    spans = []
    if not path.is_file():
        return spans
    with path.open() as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 3:
                spans.append((int(cols[1]), int(cols[2])))
    return spans


def parse_maf(path: Path) -> list[tuple[str, int, int, str]]:
    """Returns [(read_id, true_start, true_end, true_strand), ...]."""
    out = []
    ref_start = ref_len = None
    with path.open() as f:
        for line in f:
            if not line.startswith("s "):
                continue
            cols = line.split()
            # s <name> <start> <alnlen> <strand> <src_size> <seq>
            name, start, alnlen, strand = cols[1], int(cols[2]), int(cols[3]), cols[4]
            if name == "ref":
                ref_start, ref_len = start, alnlen
            else:
                assert ref_start is not None, "read line before its ref line"
                out.append((name, ref_start, ref_start + ref_len, strand))
                ref_start = ref_len = None
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset-dir", required=True, type=Path, help="existing data/<name>/ dir (needs genome.fa)")
    ap.add_argument("--pbsim-bin", default="pbsim")
    ap.add_argument("--model", required=True, help="path to a pbsim2 .model file, e.g. .../data/R94.model")
    ap.add_argument("--depth", type=float, default=8.0)
    ap.add_argument("--length-min", type=int, default=2000)
    ap.add_argument("--length-max", type=int, default=10000)
    ap.add_argument("--difference-ratio", default="23:31:46", help="sub:ins:del, ONT-like by default")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args(argv)

    genome = args.dataset_dir / "genome.fa"
    if not genome.is_file():
        print(f"error: {genome} not found (run make_dataset.py first)", file=sys.stderr)
        return 1

    work = args.dataset_dir / "_pbsim_work"
    work.mkdir(exist_ok=True)
    prefix = work / "sd"
    subprocess.run(
        [args.pbsim_bin, f"--prefix={prefix}", f"--depth={args.depth}",
         f"--length-min={args.length_min}", f"--length-max={args.length_max}",
         f"--difference-ratio={args.difference_ratio}", f"--seed={args.seed}",
         f"--hmm_model={args.model}", str(genome)],
        check=True,
    )

    fastqs = sorted(work.glob("sd_*.fastq"))
    mafs = sorted(work.glob("sd_*.maf"))
    if not fastqs:
        print("error: pbsim produced no .fastq output", file=sys.stderr)
        return 1

    reads_fa = args.dataset_dir / "reads.fa"
    n_written = 0
    with reads_fa.open("w") as out:
        for fq in fastqs:
            lines = fq.read_text().splitlines()
            for i in range(0, len(lines), 4):
                name, seq = lines[i][1:], lines[i + 1]
                out.write(f">{name}\n{seq}\n")
                n_written += 1

    truth_rows: list[tuple[str, int, int, str]] = []
    for maf in mafs:
        truth_rows.extend(parse_maf(maf))

    repeat_spans = load_repeats(args.dataset_dir / "repeats.bed")
    truth_tsv = args.dataset_dir / "truth.tsv"
    n_in_repeat = 0
    with truth_tsv.open("w") as f:
        f.write("read_id\ttrue_contig\ttrue_start\ttrue_end\ttrue_strand\tin_repeat\n")
        for rid, start, end, strand in truth_rows:
            in_repeat = any(rs <= start < re_ or rs < end <= re_ for rs, re_ in repeat_spans)
            n_in_repeat += in_repeat
            f.write(f"{rid}\tsmoke_genome\t{start}\t{end}\t{strand}\t{1 if in_repeat else 0}\n")

    shutil.rmtree(work)
    print(
        f"wrote {reads_fa} ({n_written} reads) and {truth_tsv} "
        f"({len(truth_rows)} truth rows, {n_in_repeat} landed in the repeat by chance)",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
