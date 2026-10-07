#!/usr/bin/env python3
"""Simulate reads for an existing dataset (genome.fa + repeats.bed, from
make_dataset.py) with a real long-read simulator, replacing make_dataset.py's
own lightweight built-in simulator's reads.fa/truth.tsv:

- --simulator pbsim1 (PBSIM 1.0.3, Ono et al. 2013; the simulator the
  Winnowmap paper used for its D1/D2): model-based CLR (--model = its
  model_qc_clr file). Its CLR defaults are what the paper ran with besides
  the stated mean length / error rate: log-normal lengths sd 2300, min 100,
  max 25,000 (hard cap), accuracy sd 0.02, sub:ins:del 10:60:30.
- --simulator pbsim2 (default; HMM models, e.g. P6C4.model / R94.model).
Both write the same sd_NNNN.{fastq,maf} per reference sequence.

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

Reads containing "NN" are dropped (both reads.fa and truth.tsv), exactly as
minimap2's `paftools.js pbsim2fq` does when it turns PBSIM output into reads
with truth names -- the tool the Winnowmap paper used to evaluate with
paftools. PBSIM copies reference N bases into the read, so these are reads
simulated from assembly gaps (e.g. GRCh38's chr1 has a ~17 Mbp N block at
125-142 Mbp): unsequenced placeholder, not mappable sequence. Without this the
paper's D2 would not show 0.00% unmapped. --keep-nn-reads keeps them.
"""
from __future__ import annotations

import argparse
import csv
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


def overlaps_any(start: int, end: int, spans: list[tuple[int, int]]) -> bool:
    return any(start < s_end and end > s_start for s_start, s_end in spans)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset-dir", required=True, type=Path, help="existing data/<name>/ dir (needs genome.fa)")
    ap.add_argument("--pbsim-bin", default="pbsim")
    ap.add_argument("--simulator", choices=("pbsim2", "pbsim1"), default="pbsim2")
    ap.add_argument("--model", required=True, help="pbsim2: a .model file (e.g. .../data/R94.model); pbsim1: model_qc_clr")
    ap.add_argument("--depth", type=float, default=8.0)
    # Length/accuracy options are passed only when given (all of them together,
    # e.g. mean+sd+min+max); otherwise the simulator's own defaults apply.
    ap.add_argument("--length-mean", type=float, default=None)
    ap.add_argument("--length-sd", type=float, default=None)
    ap.add_argument("--length-min", type=int, default=None)
    ap.add_argument("--length-max", type=int, default=None)
    ap.add_argument("--accuracy-mean", type=float, default=None, help="e.g. 0.90 for a paper's stated '10%% mean error-rate'")
    ap.add_argument("--accuracy-sd", type=float, default=None)
    ap.add_argument("--difference-ratio", default=None, help="sub:ins:del (simulator default if omitted: pbsim1 CLR 10:60:30, pbsim2 6:50:54); pbsim2's ONT-like mix is 23:31:46")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--keep-nn-reads", action="store_true", help="keep reads containing 'NN' (default: drop them, as paftools pbsim2fq does)")
    args = ap.parse_args(argv)

    genome = args.dataset_dir / "genome.fa"
    if not genome.is_file():
        print(f"error: {genome} not found (run make_dataset.py first)", file=sys.stderr)
        return 1

    work = args.dataset_dir / "_pbsim_work"
    work.mkdir(exist_ok=True)
    prefix = work / "sd"
    cmd = [args.pbsim_bin, f"--prefix={prefix}", f"--depth={args.depth}", f"--seed={args.seed}"]
    if args.simulator == "pbsim1":
        cmd += ["--data-type=CLR", f"--model_qc={args.model}"]
    else:
        cmd.append(f"--hmm_model={args.model}")
    for flag, val in (("length-mean", args.length_mean), ("length-sd", args.length_sd),
                      ("length-min", args.length_min), ("length-max", args.length_max),
                      ("accuracy-mean", args.accuracy_mean), ("accuracy-sd", args.accuracy_sd),
                      ("difference-ratio", args.difference_ratio)):
        if val is not None:
            cmd.append(f"--{flag}={val}")
    cmd.append(str(genome))
    print(" ".join(cmd), file=sys.stderr)
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
    # the same simulator-preserved order).
    truth_rows: list[tuple[str, str, int, int, str]] = []
    for contig, maf in zip(names, mafs):
        for rid, start, end, strand in parse_maf(maf):
            truth_rows.append((rid, contig, start, end, strand))

    # paftools.js pbsim2fq rule: drop any read whose sequence contains "NN".
    nn_ids: set[str] = set()
    if not args.keep_nn_reads:
        for fq in fastqs:
            lines = fq.read_text().splitlines()
            for i in range(0, len(lines), 4):
                if "NN" in lines[i + 1]:
                    nn_ids.add(lines[i][1:])

    repeat_spans = load_repeats(args.dataset_dir / "repeats.bed")

    reads_fa = args.dataset_dir / "reads.fa"
    n_written = 0
    with reads_fa.open("w") as out:
        for fq in fastqs:
            lines = fq.read_text().splitlines()
            for i in range(0, len(lines), 4):
                name, seq = lines[i][1:], lines[i + 1]
                if name in nn_ids:
                    continue
                out.write(f">{name}\n{seq}\n")
                n_written += 1

    truth_tsv = args.dataset_dir / "truth.tsv"
    n_in_repeat = 0
    n_kept = 0
    with truth_tsv.open("w") as f:
        f.write("read_id\ttrue_contig\ttrue_start\ttrue_end\ttrue_strand\tin_repeat\n")
        for rid, contig, start, end, strand in truth_rows:
            if rid in nn_ids:
                continue
            n_kept += 1
            in_repeat = overlaps_any(start, end, repeat_spans.get(contig, []))
            n_in_repeat += in_repeat
            f.write(f"{rid}\t{contig}\t{start}\t{end}\t{strand}\t{1 if in_repeat else 0}\n")

    shutil.rmtree(work)
    print(
        f"wrote {reads_fa} ({n_written} reads) and {truth_tsv} ({n_kept} truth rows, "
        f"{n_in_repeat} landed in a repeat by chance, {len(nn_ids)} reads containing 'NN' "
        f"dropped out of {len(truth_rows)} simulated)",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
