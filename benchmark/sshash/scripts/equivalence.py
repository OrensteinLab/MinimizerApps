"""Fixed positive/negative k-mer query sets + the order/seed equivalence gate.

Generated once per dataset from a throwaway *native*-order reference index
(sampling real k-mers from the input FASTA as positives, and random
33-mers verified absent via `sshash query` as negatives), then re-checked
against every other order/seed's index for that dataset with the exact same
fixed sets -- so "does this order answer membership queries the same way"
is tested with the same yardstick everywhere, not re-derived per run.

See export.py:validate_orders_equivalence for the hard-fail gate that reads
back the `query_pos_match`/`query_neg_match`/`num_kmers` fields this module
(via lib/run.sh) writes into each run's meta.json.
"""

from __future__ import annotations

import argparse
import json
import random
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ALPHABET = "ACGT"
_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _read_fasta_sequences(path: Path) -> list[str]:
    """Read a (gzip-or-not) FASTA's sequence lines, uppercase, one per record line."""
    opener = _opener_for(path)
    seqs: list[str] = []
    with opener(path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seqs.append(line.upper())
    return seqs


def _opener_for(path: Path):
    import gzip

    return gzip.open if str(path).endswith(".gz") else open


def sample_positive_kmers(fasta_path: Path, k: int, n_pos: int, seed: int) -> list[str]:
    """Distinct valid (ACGT-only) k-mers actually present in the input, sampled."""
    rng = random.Random(seed)
    kmers: set[str] = set()
    for s in _read_fasta_sequences(fasta_path):
        for i in range(len(s) - k + 1):
            sub = s[i : i + k]
            if all(c in "ACGT" for c in sub):
                kmers.add(sub)
                if len(kmers) >= n_pos * 50:
                    # Enough distinct k-mers to subsample from without walking
                    # the entire input (matters for the ~2.5B-kmer human set).
                    break
        if len(kmers) >= n_pos * 50:
            break
    pool = list(kmers)
    if not pool:
        raise RuntimeError(f"no valid {k}-mers found in {fasta_path}")
    rng.shuffle(pool)
    return pool[:n_pos]


def _run_sshash_query(sshash_bin: str, index_path: Path, query_fasta: Path) -> dict:
    proc = subprocess.run(
        [sshash_bin, "query", "-i", str(index_path), "-q", str(query_fasta), "--verbose"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    combined = (proc.stdout or "") + (proc.stderr or "")
    m = re.search(r"\{.*\}", combined)
    if not m:
        raise RuntimeError(
            f"sshash query produced no JSON summary (exit {proc.returncode}): {combined[-2000:]}"
        )
    return json.loads(m.group(0))


def write_kmers_fasta(kmers: list[str], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        for i, kmer in enumerate(kmers):
            f.write(f">{i}\n{kmer}\n")


def sample_negative_kmers(
    sshash_bin: str,
    index_path: Path,
    positive_pool: set[str],
    k: int,
    n_neg: int,
    seed: int,
    tmp_dir: Path,
) -> list[str]:
    """Random k-mers verified absent from `index_path` via `sshash query`.

    Verification is via the index's own aggregate query counts (batched, not
    one-by-one) -- collisions between a random k-mer and a real one are
    astronomically unlikely at k=31 (4^31 space), so a batch is normally
    accepted whole; any batch that doesn't verify falls back to a slower
    per-kmer check to isolate exactly which candidates were actually present.
    """
    rng = random.Random(seed ^ 0x5EED)
    canon_positive = {p for kmer in positive_pool for p in (kmer, revcomp(kmer))}
    negatives: list[str] = []
    attempt = 0
    while len(negatives) < n_neg:
        attempt += 1
        need = n_neg - len(negatives)
        batch: list[str] = []
        seen_batch: set[str] = set()
        while len(batch) < need:
            cand = "".join(rng.choice(ALPHABET) for _ in range(k))
            if cand in canon_positive or cand in seen_batch:
                continue
            seen_batch.add(cand)
            batch.append(cand)
        qfa = tmp_dir / f"neg_candidates_{attempt}.fa"
        write_kmers_fasta(batch, qfa)
        result = _run_sshash_query(sshash_bin, index_path, qfa)
        qfa.unlink(missing_ok=True)
        if int(result["num_negative_kmers"]) == int(result["num_kmers"]):
            negatives.extend(batch)
            continue
        # Fallback: isolate which candidates in this batch were actually present.
        for cand in batch:
            qfa1 = tmp_dir / "neg_one.fa"
            write_kmers_fasta([cand], qfa1)
            r1 = _run_sshash_query(sshash_bin, index_path, qfa1)
            qfa1.unlink(missing_ok=True)
            if int(r1["num_negative_kmers"]) == 1:
                negatives.append(cand)
            if len(negatives) >= n_neg:
                break
    return negatives[:n_neg]


def build_reference_index(
    sshash_bin: str, fasta_path: Path, k: int, m: int, threads: int, out_path: Path
) -> None:
    cmd = [
        sshash_bin,
        "build",
        "-i",
        str(fasta_path),
        "-k",
        str(k),
        "-m",
        str(m),
        "-O",
        "native",
        "-t",
        str(threads),
        "-o",
        str(out_path),
    ]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"reference build failed (exit {proc.returncode}): {(proc.stdout or '') + (proc.stderr or '')}"
        )


def ensure_query_sets(
    *,
    sshash_bin: str,
    fasta_path: Path,
    k: int,
    m: int,
    threads: int,
    out_dir: Path,
    n_pos: int = 10000,
    n_neg: int = 10000,
    seed: int = 1,
    force: bool = False,
) -> tuple[Path, Path]:
    """Generate (or reuse) `{out_dir}/pos.txt` and `{out_dir}/neg.txt`."""
    pos_path = out_dir / "pos.txt"
    neg_path = out_dir / "neg.txt"
    if pos_path.is_file() and neg_path.is_file() and not force:
        return pos_path, neg_path

    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_root = Path(tempfile.mkdtemp(prefix="sshash_refidx_"))
    try:
        ref_index = tmp_root / "reference_native.sshash"
        build_reference_index(sshash_bin, fasta_path, k, m, threads, ref_index)
        positives = sample_positive_kmers(fasta_path, k, n_pos, seed)
        negatives = sample_negative_kmers(
            sshash_bin, ref_index, set(positives), k, n_neg, seed, tmp_root
        )
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)

    pos_path.write_text("\n".join(positives) + "\n")
    neg_path.write_text("\n".join(negatives) + "\n")
    return pos_path, neg_path


def kmers_txt_to_fasta(txt_path: Path, fasta_path: Path) -> Path:
    """Cache a k-mer-per-line txt file as a FASTA sshash can query directly."""
    if fasta_path.is_file() and fasta_path.stat().st_mtime >= txt_path.stat().st_mtime:
        return fasta_path
    kmers = [ln.strip() for ln in txt_path.read_text().splitlines() if ln.strip()]
    write_kmers_fasta(kmers, fasta_path)
    return fasta_path


def check_query_file(sshash_bin: str, index_path: Path, query_fasta: Path, expect: str) -> dict:
    """Query `query_fasta` against `index_path`; `expect` is "positive" or "negative"."""
    result = _run_sshash_query(sshash_bin, index_path, query_fasta)
    total = int(result["num_kmers"])
    n_pos = int(result["num_positive_kmers"])
    n_neg = int(result["num_negative_kmers"])
    n_invalid = int(result["num_invalid_kmers"])
    if expect == "positive":
        match = n_invalid == 0 and n_pos == total and n_neg == 0
    elif expect == "negative":
        match = n_invalid == 0 and n_neg == total and n_pos == 0
    else:
        raise ValueError(f"expect must be 'positive' or 'negative', got {expect!r}")
    return {
        "match": match,
        "num_kmers": total,
        "num_positive_kmers": n_pos,
        "num_negative_kmers": n_neg,
        "num_invalid_kmers": n_invalid,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="scripts equivalence", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_ensure = sub.add_parser("ensure-queries", help="generate (or reuse) fixed pos/neg query sets")
    p_ensure.add_argument("--sshash", required=True)
    p_ensure.add_argument("--fasta", required=True)
    p_ensure.add_argument("--k", type=int, required=True)
    p_ensure.add_argument("--m", type=int, required=True)
    p_ensure.add_argument("--threads", type=int, default=8)
    p_ensure.add_argument("--out-dir", required=True)
    p_ensure.add_argument("--n-pos", type=int, default=10000)
    p_ensure.add_argument("--n-neg", type=int, default=10000)
    p_ensure.add_argument("--seed", type=int, default=1)
    p_ensure.add_argument("--force", action="store_true")

    p_tofasta = sub.add_parser("to-fasta", help="cache a k-mer txt list as a query FASTA")
    p_tofasta.add_argument("--txt", required=True)
    p_tofasta.add_argument("--fasta", required=True)

    p_check = sub.add_parser("check", help="query a fixed set against one index; prints JSON")
    p_check.add_argument("--sshash", required=True)
    p_check.add_argument("--index", required=True)
    p_check.add_argument("--query-fasta", required=True)
    p_check.add_argument("--expect", choices=("positive", "negative"), required=True)

    args = ap.parse_args(argv)

    if args.cmd == "ensure-queries":
        pos_path, neg_path = ensure_query_sets(
            sshash_bin=args.sshash,
            fasta_path=Path(args.fasta),
            k=args.k,
            m=args.m,
            threads=args.threads,
            out_dir=Path(args.out_dir),
            n_pos=args.n_pos,
            n_neg=args.n_neg,
            force=args.force,
        )
        print(json.dumps({"pos": str(pos_path), "neg": str(neg_path)}))
        return 0

    if args.cmd == "to-fasta":
        out = kmers_txt_to_fasta(Path(args.txt), Path(args.fasta))
        print(str(out))
        return 0

    if args.cmd == "check":
        result = check_query_file(args.sshash, Path(args.index), Path(args.query_fasta), args.expect)
        print(json.dumps(result))
        return 0 if result["match"] else 1

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
