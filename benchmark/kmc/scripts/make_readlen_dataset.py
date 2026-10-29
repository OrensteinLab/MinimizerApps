#!/usr/bin/env python3
"""Build a derived dataset for the read-length confound check.

GreedyMini's wall-time win only reproduced reliably on f_vesca, which is both
the smallest dataset (~4.1B k-mers) *and* the one with by far the longest
reads (avg ~350bp vs ~64-101bp for everything else). With only 5 datasets,
size and read length are confounded and it's not possible to tell which one
actually matters from the existing runs alone (see README: "Read-length
confound check").

This script concatenates every `--join` consecutive FASTQ records from an
already-downloaded dataset into one pseudo-read, optionally capping the
output at `--max-output-reads` pseudo-reads. That gives two controlled
variants built from the same source data:

  --join 1 --max-output-reads <f_vesca's read count>
      Native read length, f_vesca's *read count* -> isolates size.

  --join 4 --max-output-reads <count matching f_vesca's total bases>
      ~4x native read length, f_vesca's *total bases* -> isolates length
      (this is the direct "does GreedyMini's win reappear" replication).

This is a mechanical transform for stressing KMC's splitter/signature code
path with a given read length, not simulated long-read sequencing: quality
scores are just concatenated and the sequence at each join point has no
biological meaning. Only read length and total size are controlled; treat
any k-mer-content difference from a "real" long-read dataset as expected.

The output is registered as a new dataset id under the same
$KMC3_DATA_DIR (= $KMC3_BENCH_ROOT/data) tree that `download`-ed datasets use, so
`kmc3_bench.sh run -d <out-id> ...` works unmodified. Derived datasets are
not in data/datasets.csv, so there is no mem_gb to look up -- pass --mem
explicitly when running kmc3_bench.sh against one.
"""

from __future__ import annotations

import argparse
import gzip
import math
import sys
from pathlib import Path

from scripts.common import die, load_dotenv
from paths import bench_paths


def _open_maybe_gz(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def iter_fastq_records(paths: list[Path]):
    """Yield (header, seq, qual) across all files, back to back."""
    for p in paths:
        with _open_maybe_gz(p) as fh:
            while True:
                header = fh.readline()
                if not header:
                    break
                seq = fh.readline()
                plus = fh.readline()
                qual = fh.readline()
                if not (seq and plus and qual):
                    break
                yield header.rstrip("\n"), seq.rstrip("\n"), qual.rstrip("\n")


def probe_avg_read_len(source_paths: list[Path], n_probe: int) -> float:
    """Average native record length over the first `n_probe` source records."""
    total = 0
    n = 0
    for _header, seq, _qual in iter_fastq_records(source_paths):
        total += len(seq)
        n += 1
        if n >= n_probe:
            break
    if n == 0:
        die("source has no records to probe")
    return total / n


def max_output_reads_for_target(
    source_paths: list[Path], join_n: int, k: int, target_kmers: int, n_probe: int
) -> int:
    """Pseudo-reads needed so join_n*probe_avg_len - k + 1 per read reaches target_kmers."""
    avg_native_len = probe_avg_read_len(source_paths, n_probe)
    kmers_per_output_read = join_n * avg_native_len - k + 1
    if kmers_per_output_read <= 0:
        die(
            f"join={join_n} x avg native read length {avg_native_len:.1f} is shorter than "
            f"k={k}; raise --join"
        )
    n_reads = math.ceil(target_kmers / kmers_per_output_read)
    print(
        f"probed avg native read length {avg_native_len:.1f} over {n_probe} records -> "
        f"{n_reads} pseudo-reads needed for ~{target_kmers} k-mers at join={join_n}, k={k}"
    )
    return n_reads


def output_paths_for(out_path: Path, num_files: int) -> list[Path]:
    """Single `out_path` for num_files=1, else `{stem}_part{i:02d}{suffixes}`."""
    if num_files == 1:
        return [out_path]
    stem = out_path.name.split(".", 1)[0]
    suffix = out_path.name[len(stem):]  # e.g. ".fastq.gz"
    width = len(str(num_files - 1))
    return [out_path.with_name(f"{stem}_part{i:0{width}d}{suffix}") for i in range(num_files)]


def build(
    source_paths: list[Path],
    out_path: Path,
    join_n: int,
    max_output_reads: int | None,
    num_files: int = 1,
) -> tuple[int, int, int, list[Path]]:
    """Write joined pseudo-reads, round-robin across `num_files` gzip files
    (matching multi-file native downloads so KMC's readers can decompress in
    parallel; see README "Read-length confound check" I/O-parallelism note).
    Return (n_input_records, n_output_reads, total_bases, out_paths_written)."""
    if num_files < 1:
        die("num_files must be >= 1")
    out_paths = output_paths_for(out_path, num_files)
    n_in = 0
    n_out = 0
    total_bases = 0
    group: list[tuple[str, str]] = []
    handles = [gzip.open(p, "wt", compresslevel=6) for p in out_paths]
    try:
        for _header, seq, qual in iter_fastq_records(source_paths):
            n_in += 1
            group.append((seq, qual))
            if len(group) < join_n:
                continue
            seq_j = "".join(s for s, _ in group)
            qual_j = "".join(q for _, q in group)
            group = []
            out = handles[n_out % num_files]
            n_out += 1
            out.write(f"@readlen_join{join_n}_{n_out}\n{seq_j}\n+\n{qual_j}\n")
            total_bases += len(seq_j)
            if max_output_reads and n_out >= max_output_reads:
                break
    finally:
        for h in handles:
            h.close()
    return n_in, n_out, total_bases, out_paths


def main(argv: list[str] | None = None) -> int:
    load_dotenv()
    ap = argparse.ArgumentParser(prog="scripts make-readlen-dataset", description=__doc__)
    ap.add_argument("--source", required=True, help="existing (already downloaded) dataset id")
    ap.add_argument("--out-id", required=True, help="new dataset id to create under data/")
    ap.add_argument(
        "--join", type=int, required=True,
        help="consecutive source reads concatenated into 1 pseudo-read (1 = passthrough/subsample only)",
    )
    ap.add_argument(
        "--max-output-reads", type=int, default=None,
        help="stop after this many pseudo-reads (controls total dataset size); "
        "alternative to --target-kmers",
    )
    ap.add_argument(
        "--target-kmers", type=int, default=None,
        help="aim for this many total k-mers instead of a literal --max-output-reads: "
        "probes the source's native avg read length and computes the pseudo-read count "
        "automatically from --join and --k (no manual arithmetic needed)",
    )
    ap.add_argument("--k", type=int, default=28, help="k, only used with --target-kmers (default: 28)")
    ap.add_argument(
        "--probe-reads", type=int, default=200_000,
        help="source records sampled to estimate avg native read length for --target-kmers (default: 200000)",
    )
    ap.add_argument("--input", choices=("gzip", "plain"), default="gzip", help="which source file list to read")
    ap.add_argument(
        "--num-files", type=int, default=1,
        help="split output round-robin across this many gzip files instead of one "
        "(default: 1); use to match a native download's file count so KMC's readers "
        "can decompress in parallel -- isolates dataset content from single-file I/O "
        "throttling, see README 'Read-length confound check'",
    )
    ap.add_argument("--data-dir", default=None, help="override KMC3_DATA_DIR (default $KMC3_BENCH_ROOT/data)")
    ap.add_argument("--force", action="store_true", help="overwrite an existing out-id")
    args = ap.parse_args(argv)

    if args.join < 1:
        die("--join must be >= 1")
    if args.num_files < 1:
        die("--num-files must be >= 1")
    if args.max_output_reads is not None and args.target_kmers is not None:
        die("pass only one of --max-output-reads / --target-kmers")

    data_root = Path(args.data_dir) if args.data_dir else bench_paths("kmc")["data"]

    list_name = "files_gz.lst" if args.input == "gzip" else "files_plain.lst"
    src_list = data_root / args.source / list_name
    if not src_list.is_file():
        die(f"missing {src_list}; run: benchmark/kmc3_bench.sh download -d {args.source}")
    src_paths = [Path(p) for p in src_list.read_text().splitlines() if p.strip()]
    if not src_paths:
        die(f"{src_list} is empty")
    for p in src_paths:
        if not p.is_file():
            die(f"missing source file {p}")

    max_output_reads = args.max_output_reads
    if args.target_kmers is not None:
        max_output_reads = max_output_reads_for_target(
            src_paths, args.join, args.k, args.target_kmers, args.probe_reads
        )

    out_dir = data_root / args.out_id / "fastq_gz"
    out_path = out_dir / f"{args.out_id}.fastq.gz"
    out_list = data_root / args.out_id / "files_gz.lst"
    planned_paths = output_paths_for(out_path, args.num_files)
    existing = [p for p in planned_paths if p.exists()]
    if existing and not args.force:
        die(f"{existing[0]} already exists; pass --force to overwrite")
    out_dir.mkdir(parents=True, exist_ok=True)

    n_in, n_out, total_bases, out_paths = build(
        src_paths, out_path, args.join, max_output_reads, args.num_files
    )
    if n_out == 0:
        die("no output reads written (source empty, or --join larger than available input records)")

    out_list.write_text("".join(f"{p.resolve()}\n" for p in out_paths))

    avg_len = total_bases / n_out
    print(f"wrote {len(out_paths)} file(s) under {out_dir}")
    print(f"  {n_in} source records -> {n_out} pseudo-reads (join={args.join})")
    print(f"  avg_read_len={avg_len:.1f}  total_bases={total_bases}")
    print(f"  files_gz.lst -> {out_list} ({len(out_paths)} file(s))")
    print()
    print("next (pass --mem explicitly; derived datasets aren't in data/datasets.csv):")
    print(
        f"  benchmark/kmc3_bench.sh run -d {args.out_id} -k 28 -S all -n 5 "
        "--mem 12 --bins 512 -t 12 --force"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
