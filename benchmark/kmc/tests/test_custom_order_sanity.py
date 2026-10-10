#!/usr/bin/env python3
"""Sanity checks for custom signature order (greedymini table, xor, xor-pcg-hash)."""

from __future__ import annotations

import json
import os
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.config import repo_root
from scripts.expected_density import rxs_m_xs64
from scripts.gm import ensure_gm_tables, gm_binary_path

ROOT = repo_root()
GM = gm_binary_path(28, 9)
GM_DNA = ensure_gm_tables(28, 9, 0xDEADBEEF)
KMC = ROOT / "bin/kmc"
# real-reads check; skipped when KMC3_BENCH_ROOT/data/f_vesca is not available
FQ = Path(os.environ.get("KMC3_BENCH_ROOT", "/nonexistent")) / "data/f_vesca/fastq_gz/SRR072005.fastq.gz"


def dna_to_binary_key(dna: int, p: int) -> int:
    key = 0
    for i in range(p - 1, -1, -1):
        b = (dna >> (2 * i)) & 3
        key = (key << 1) | (1 if b in (1, 2) else 0)
    return key


def rank_hash(dna: int, seed: int) -> int:
    return (dna ^ seed) & ((1 << 64) - 1)


def rank_pcg(dna: int, seed: int) -> int:
    return rxs_m_xs64(dna ^ seed)


def load_gm(path: Path) -> list[int]:
    data = path.read_bytes()
    (n,) = struct.unpack_from("<Q", data, 0)
    return list(struct.unpack_from(f"<{n}Q", data, 8))


def build_rank_lut(ranks: list[int], k: int) -> list[int]:
    n = 1 << (2 * k)
    return [ranks[dna_to_binary_key(i, k)] for i in range(n)]


def table_less(i: int, j: int, rank_lut: list[int], tie_seed: int) -> bool:
    ri, rj = rank_lut[i], rank_lut[j]
    if ri != rj:
        return ri < rj
    return (i ^ tie_seed) < (j ^ tie_seed)


def build_cmp_table(ranks: list[int], k: int, seed: int) -> list[int]:
    n = 1 << (2 * k)
    cmp = [0] * n
    for i in range(n):
        key = dna_to_binary_key(i, k)
        rank = ranks[key]
        cmp[i] = (rank << 32) | (i ^ (seed & 0xFFFFFFFF))
    return cmp


def build_cmp_hash(k: int, seed: int) -> list[int]:
    n = 1 << (2 * k)
    return [rank_hash(i, seed) for i in range(n)]


def build_cmp_pcg(k: int, seed: int) -> list[int]:
    n = 1 << (2 * k)
    return [rank_pcg(i, seed) for i in range(n)]


def test_projection() -> None:
    assert dna_to_binary_key(0b0000, 2) == 0b00
    assert dna_to_binary_key(0b0010, 2) == 0b01
    assert dna_to_binary_key(0b0101, 2) == 0b11
    assert dna_to_binary_key(0b1010, 2) == 0b11
    assert dna_to_binary_key(0b1000, 2) == 0b10
    print("ok  projection AT/GC")


def test_hash_bijection() -> None:
    k, seed = 9, 42
    n = 1 << (2 * k)
    cmp = build_cmp_hash(k, seed)
    assert len(set(cmp)) == n
    assert cmp[0] == rank_hash(0, seed)
    assert cmp[100] == (100 ^ seed)
    print("ok  hash bijective, no tie-break suffix")


def test_pcg_differs_from_hash() -> None:
    k, seed = 9, 7
    h = build_cmp_hash(k, seed)
    p = build_cmp_pcg(k, seed)
    diffs = sum(1 for a, b in zip(h, p) if a != b)
    assert diffs > 0
    print(f"ok  pcg differs from hash on {diffs}/{1 << (2 * k)} entries")


def test_table_tiebreak() -> None:
    ranks = load_gm(GM)
    assert len(ranks) == 512
    k = 9
    cmp0 = build_cmp_table(ranks, k, 0)
    cmp99 = build_cmp_table(ranks, k, 99)
    changed = same_rank_diff_tie = 0
    for i in range(1 << (2 * k)):
        r0, t0 = cmp0[i] >> 32, cmp0[i] & 0xFFFFFFFF
        r99, t99 = cmp99[i] >> 32, cmp99[i] & 0xFFFFFFFF
        if cmp0[i] != cmp99[i]:
            changed += 1
        if r0 == r99 and t0 != t99:
            same_rank_diff_tie += 1
        assert t0 == (i ^ 0)
        assert t99 == (i ^ 99)
    assert changed == same_rank_diff_tie == (1 << (2 * k))
    print("ok  table tie-break seed affects all entries, rank unchanged")


def test_table_rank_lut_matches_packed() -> None:
    ranks = load_gm(GM)
    k, seed = 9, 0xDEADBEEF
    rank_lut = build_rank_lut(ranks, k)
    cmp = build_cmp_table(ranks, k, seed)
    n = 1 << (2 * k)
    import random

    rng = random.Random(0)
    for _ in range(200000):
        i = rng.randrange(n)
        j = rng.randrange(n)
        assert (cmp[i] < cmp[j]) == table_less(i, j, rank_lut, seed & 0xFFFFFFFF)
    print("ok  rank LUT + tie matches packed uint64 order (200k pairs)")


def test_table_sentinel() -> None:
    ranks = load_gm(GM)
    n = len(ranks)
    non_uhs = sum(1 for r in ranks if r >= n)
    assert non_uhs > 0
    cmp = build_cmp_table(ranks, 9, 0)
    assert any((c >> 32) >= n for c in cmp)
    print(f"ok  gm table: {n - non_uhs} UHS keys, sentinel={n}")


def run_kmc(label: str, extra: list[str], fq: Path, work: Path) -> dict:
    out = work / label
    tmp = out / "tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    json_path = out / "kmc.json"
    cmd = [str(KMC), "-ci2", "-k28", "-hp", f"-j{json_path}", *extra, str(fq), str(out / "db"), str(tmp)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"kmc failed ({label}): {proc.stderr or proc.stdout}")
    return json.loads(json_path.read_text())


def test_kmc_integration() -> None:
    if not KMC.is_file():
        print("skip kmc integration (bin/kmc missing)")
        return
    if not FQ.is_file():
        print("skip kmc integration (f_vesca sample FASTQ missing)")
        return
    if not GM_DNA.is_file():
        print("skip kmc integration (_dna.gm missing)")
        return

    with tempfile.TemporaryDirectory(prefix="kmc_order_") as td:
        work = Path(td)
        def_run = run_kmc("def", [], FQ, work)
        hash_run = run_kmc("xor", ["-orxor", "-olxor", "-os42"], FQ, work)
        pcg_run = run_kmc("xor-pcg-hash", ["-orxor-pcg-hash", "-olxor-pcg-hash", "-os7"], FQ, work)
        gm_run = run_kmc("gm", [f"-co{GM_DNA}", "-olgm"], FQ, work)
        gm0 = run_kmc("gm0", [f"-co{GM_DNA}", "-olgm0"], FQ, work)

        counts = {
            "default": def_run["Stats"]["#Unique_counted_k-mers"],
            "xor": hash_run["Stats"]["#Unique_counted_k-mers"],
            "xor-pcg-hash": pcg_run["Stats"]["#Unique_counted_k-mers"],
            "gm_s99": gm_run["Stats"]["#Unique_counted_k-mers"],
            "gm_s0": gm0["Stats"]["#Unique_counted_k-mers"],
        }
        assert len(set(counts.values())) == 1, counts
        assert hash_run["signature_order_mode"] == "xor"
        assert pcg_run["signature_order_mode"] == "xor-pcg-hash"
        assert gm_run["signature_order_mode"] == "greedymini"
        assert gm_run["signature_order_gm_lut"] == "dna"
        sk_def = def_run["Stats"]["#Total_super-k-mers"]
        sk_hash = hash_run["Stats"]["#Total_super-k-mers"]
        assert sk_hash != sk_def
        print(f"ok  kmc integration: counted={counts['default']} (all modes match)")
        print(f"     super-k-mers default={sk_def} hash={sk_hash} pcg={pcg_run['Stats']['#Total_super-k-mers']}")


def main() -> int:
    test_projection()
    test_hash_bijection()
    test_pcg_differs_from_hash()
    test_table_tiebreak()
    test_table_rank_lut_matches_packed()
    test_table_sentinel()
    test_kmc_integration()
    print("all sanity checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
