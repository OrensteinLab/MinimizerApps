#!/usr/bin/env python3
"""Simulate reads for an existing dataset (genome.fa + repeats.bed, from
make_dataset.py) using a real ONT-model long-read simulator (pbsim2, see
environment.yml/setup.sh), replacing make_dataset.py's own lightweight
built-in simulator's reads.fa/truth.tsv with pbsim2's.

pbsim2 doesn't target specific loci, so repeat-region coverage isn't forced
the way make_dataset.py's own simulator does; with --depth high enough over
a genome whose repeat is ~7% of its length (20kb/300kb here), a comparable
number of reads still land there by chance -- reported in the output.

pbsim2 runs each reference sequence as its own batch (its own stderr
"Simulation stats (ref.N)" block, own sd_NNNN.fastq/.maf pair, numbered in
genome.fa's own sequence order) -- <prefix>_NNNN.maf's `s ref <start>
<alnlen> <strand> ...` line right before each read's own `s <read_id> ...`
line carries that read's true (start, end, strand) *within that one
contig*. Parsed here into the same truth.tsv shape make_dataset.py
produces, so score_accuracy.py/lib/run.sh need no changes to consume
either -- except unlike make_dataset.py's own single-contig synthetic
genome, a real multi-contig reference (e.g. grch38_wg's 195 sequences)
means true_contig must be each read's own real source contig, not a
placeholder, both for correct repeat-overlap checking below and for
score_accuracy.py's own correctness check (a read's primary alignment must
land on the *same contig* it truly came from, not just an overlapping
coordinate range on any contig).

Reads whose true span overlaps an assembly gap (a run of N bases in the
reference) are dropped: such a read is simulated partly or entirely from
unsequenced placeholder bases, so it has no real biological sequence to be
correctly mapped back to -- not a mapping-quality question at all, just
simulator noise. This matters most for GRCh38 (whose centromeres remain
large unresolved N-blocks, unlike T2T-CHM13); e.g. grch38_wg's chr1 has an
~17Mbp, 98.6%-N block around 125-142Mbp that otherwise produces reads that
are >90% N and unmappable by construction.
"""
from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path


def load_repeats(path: Path) -> dict[str, list[tuple[int, int]]]:
    """{contig: [(start, end), ...]}, from a 3+-column BED file."""
    spans: dict[str, list[tuple[int, int]]] = {}
    if not path.is_file():
        return spans
    with path.open() as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 3:
                spans.setdefault(cols[0], []).append((int(cols[1]), int(cols[2])))
    return spans


def contig_names(genome: Path) -> list[str]:
    """Sequence names in FASTA header order -- matches pbsim2's own per-
    reference batch numbering (ref.1 = first header, ref.2 = second, ...)."""
    names = []
    with genome.open() as f:
        for line in f:
            if line.startswith(">"):
                names.append(line[1:].split()[0].strip())
    return names


def parse_maf(path: Path) -> list[tuple[str, int, int, str]]:
    """Returns [(read_id, true_start, true_end, true_strand), ...], contig-
    local coordinates (this one maf file is for a single reference contig)."""
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


def n_runs(samtools: str, genome: Path, contig: str) -> list[tuple[int, int]]:
    """[(start, end), ...] 0-based half-open spans of N (any case) in one
    contig, read once via `samtools faidx` (no .fai index write needed if
    one doesn't already exist -- faidx creates it on the fly)."""
    out = subprocess.run(
        [samtools, "faidx", str(genome), contig], capture_output=True, text=True, check=True
    ).stdout
    seq = "".join(out.splitlines()[1:])
    return [(m.start(), m.end()) for m in re.finditer(r"[Nn]+", seq)]


def overlaps_any(start: int, end: int, spans: list[tuple[int, int]]) -> bool:
    return any(start < s_end and end > s_start for s_start, s_end in spans)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset-dir", required=True, type=Path, help="existing data/<name>/ dir (needs genome.fa)")
    ap.add_argument("--pbsim-bin", default="pbsim")
    ap.add_argument("--samtools", default="samtools", help="for N-gap detection (faidx); see module docstring")
    ap.add_argument("--model", required=True, help="path to a pbsim2 .model file, e.g. .../data/R94.model")
    ap.add_argument("--depth", type=float, default=8.0)
    ap.add_argument("--length-min", type=int, default=2000, help="ignored if --length-mean is given")
    ap.add_argument("--length-max", type=int, default=10000, help="ignored if --length-mean is given")
    ap.add_argument("--length-mean", type=float, default=None, help="pbsim2's own model-based length parameterization (mutually exclusive with --length-min/--length-max) -- use to match a paper's own stated mean read length exactly, e.g. Winnowmap's D2 (15kbp)")
    ap.add_argument("--length-sd", type=float, default=None, help="paired with --length-mean; pbsim2 default (7000.0) if omitted")
    ap.add_argument("--accuracy-mean", type=float, default=None, help="pbsim2's own model-based mean accuracy (e.g. 0.90 for a paper's stated '10% mean error-rate'); pbsim2 default (0.85) if omitted")
    ap.add_argument("--difference-ratio", default="23:31:46", help="sub:ins:del; ONT-like by default, pass 6:50:54 (pbsim2's own default) for PacBio-style reads")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--keep-n-gap-reads", action="store_true", help="skip N-gap filtering (debugging only)")
    args = ap.parse_args(argv)

    genome = args.dataset_dir / "genome.fa"
    if not genome.is_file():
        print(f"error: {genome} not found (run make_dataset.py first)", file=sys.stderr)
        return 1

    work = args.dataset_dir / "_pbsim_work"
    work.mkdir(exist_ok=True)
    prefix = work / "sd"
    cmd = [args.pbsim_bin, f"--prefix={prefix}", f"--depth={args.depth}",
           f"--difference-ratio={args.difference_ratio}", f"--seed={args.seed}",
           f"--hmm_model={args.model}"]
    if args.length_mean is not None:
        cmd.append(f"--length-mean={args.length_mean}")
        if args.length_sd is not None:
            cmd.append(f"--length-sd={args.length_sd}")
    else:
        cmd += [f"--length-min={args.length_min}", f"--length-max={args.length_max}"]
    if args.accuracy_mean is not None:
        cmd.append(f"--accuracy-mean={args.accuracy_mean}")
    cmd.append(str(genome))
    subprocess.run(cmd, check=True)

    fastqs = sorted(work.glob("sd_*.fastq"))
    mafs = sorted(work.glob("sd_*.maf"))
    if not fastqs:
        print("error: pbsim produced no .fastq output", file=sys.stderr)
        return 1

    names = contig_names(genome)
    if len(names) != len(mafs):
        print(
            f"error: genome.fa has {len(names)} sequences but pbsim produced {len(mafs)} "
            "ref batches -- contig<->maf pairing would be wrong, refusing to guess",
            file=sys.stderr,
        )
        return 1

    # (read_id, contig, start, end, strand) -- contig-local coordinates,
    # contig assigned positionally (maf i <-> genome.fa sequence i, both in
    # the same pbsim2-preserved order).
    truth_rows: list[tuple[str, str, int, int, str]] = []
    for contig, maf in zip(names, mafs):
        for rid, start, end, strand in parse_maf(maf):
            truth_rows.append((rid, contig, start, end, strand))

    n_gap_ids: set[str] = set()
    if not args.keep_n_gap_reads:
        by_contig: dict[str, list[tuple[str, int, int]]] = {}
        for rid, contig, start, end, _strand in truth_rows:
            by_contig.setdefault(contig, []).append((rid, start, end))
        for contig, reads in by_contig.items():
            gaps = n_runs(args.samtools, genome, contig)
            if not gaps:
                continue
            for rid, start, end in reads:
                if overlaps_any(start, end, gaps):
                    n_gap_ids.add(rid)

    repeat_spans = load_repeats(args.dataset_dir / "repeats.bed")

    reads_fa = args.dataset_dir / "reads.fa"
    n_written = 0
    with reads_fa.open("w") as out:
        for fq in fastqs:
            lines = fq.read_text().splitlines()
            for i in range(0, len(lines), 4):
                name, seq = lines[i][1:], lines[i + 1]
                if name in n_gap_ids:
                    continue
                out.write(f">{name}\n{seq}\n")
                n_written += 1

    truth_tsv = args.dataset_dir / "truth.tsv"
    n_in_repeat = 0
    n_kept = 0
    with truth_tsv.open("w") as f:
        f.write("read_id\ttrue_contig\ttrue_start\ttrue_end\ttrue_strand\tin_repeat\n")
        for rid, contig, start, end, strand in truth_rows:
            if rid in n_gap_ids:
                continue
            n_kept += 1
            in_repeat = overlaps_any(start, end, repeat_spans.get(contig, []))
            n_in_repeat += in_repeat
            f.write(f"{rid}\t{contig}\t{start}\t{end}\t{strand}\t{1 if in_repeat else 0}\n")

    shutil.rmtree(work)
    print(
        f"wrote {reads_fa} ({n_written} reads) and {truth_tsv} ({n_kept} truth rows, "
        f"{n_in_repeat} landed in a repeat by chance, {len(n_gap_ids)} N-gap-overlapping "
        f"reads dropped out of {len(truth_rows)} simulated)",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
