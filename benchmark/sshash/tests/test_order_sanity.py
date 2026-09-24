#!/usr/bin/env python3
"""Offline sanity checks for sshash's custom minimizer orders (no network,
no big downloads -- uses the in-repo sshash sample dataset only).

Mirrors benchmark/kmc/tests/test_custom_order_sanity.py::test_kmc_integration:
build the same input under every order, assert they all report the exact
same num_kmers, and assert every order's index answers the exact same fixed
positive/negative k-mer queries (see data/queries/) identically.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))

from scripts.config import monorepo_root, repo_root
from scripts.orders import ORDERS, build_flags, order_t_values

BENCH_DIR = Path(__file__).resolve().parents[1]
SSHASH = repo_root() / "build" / "sshash"
SAMPLE = repo_root() / "data" / "unitigs_stitched" / "salmonella_enterica_k31_ust.fa.gz"
# Plain greedymini's table (gm_k=15, Extend-k'd up to l=m=21).
GM_FILE = monorepo_root() / "orders" / "GreedyMini" / "w11_k15_canonical.gm"
# mod-minimizer-greedymini's O_t table (t=11, no extend needed: gm_k=11=t).
MOD_GM_FILE = monorepo_root() / "orders" / "GreedyMini" / "w21_k11_canonical.gm"
POS_TXT = BENCH_DIR / "data" / "queries" / "salmonella_enterica_k31_pos.txt"
NEG_TXT = BENCH_DIR / "data" / "queries" / "salmonella_enterica_k31_neg.txt"
# m=21 (up from the original 11): matches the mod-minimizer paper's own
# tested SSHash regime (k=31, m=21, w=11) -- see orders/README.md.
K, M = 31, 21
EXPECTED_NUM_KMERS = 4787534  # independent of m: m only affects minimizer construction, not the k-mer set

SEEDS = {
    "native": 1,
    "xor-hash": 42,
    "xor-pcg-hash": 7,
    "greedymini": 99,
    "mod-minimizer-native-hash": 11,
    "mod-minimizer-xor-hash": 12,
    "mod-minimizer-pcg-xor-hash": 13,
    "mod-minimizer-greedymini": 14,
}


def _gm_file_for(order: str) -> str | None:
    if order == "greedymini":
        return str(GM_FILE)
    if order == "mod-minimizer-greedymini":
        return str(MOD_GM_FILE)
    return None


def _run_label(order: str, t: int | None) -> str:
    return order if t is None else f"{order}_t{t}"


def _run(cmd: list[str]) -> tuple[int, str]:
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    return proc.returncode, (proc.stdout or "") + (proc.stderr or "")


def _json_blob(text: str) -> dict:
    m = re.search(r"\{.*\}", text, re.DOTALL)
    assert m, f"no JSON blob found in output:\n{text[-2000:]}"
    return json.loads(m.group(0))


def _build(order: str, t: int | None, out_path: Path) -> dict:
    flags = build_flags(order, SEEDS[order], _gm_file_for(order), t)
    cmd = [str(SSHASH), "build", "-i", str(SAMPLE), "-k", str(K), "-m", str(M), "-t", "4", "--verbose", *flags, "-o", str(out_path)]
    rc, out = _run(cmd)
    assert rc == 0, f"build failed for order={_run_label(order, t)}: {out[-2000:]}"
    return _json_blob(out)


def _check(index_path: Path) -> None:
    rc, out = _run([str(SSHASH), "check", "-i", str(index_path), "--verbose"])
    assert rc == 0 and "EVERYTHING OK!" in out, f"check failed for {index_path}: {out[-2000:]}"


def _kmers_fasta(txt_path: Path, out_fasta: Path) -> None:
    kmers = [ln.strip() for ln in txt_path.read_text().splitlines() if ln.strip()]
    with open(out_fasta, "w") as f:
        for i, kmer in enumerate(kmers):
            f.write(f">{i}\n{kmer}\n")


def _query(index_path: Path, query_fasta: Path) -> dict:
    rc, out = _run([str(SSHASH), "query", "-i", str(index_path), "-q", str(query_fasta), "--verbose"])
    assert rc == 0, f"query failed for {index_path}: {out[-2000:]}"
    return _json_blob(out)


def test_orders_module() -> None:
    assert set(ORDERS) == {
        "native", "xor-hash", "xor-pcg-hash", "greedymini",
        "mod-minimizer-native-hash", "mod-minimizer-xor-hash",
        "mod-minimizer-pcg-xor-hash", "mod-minimizer-greedymini",
    }
    assert build_flags("native", 1) == ["-s", "1"]
    assert build_flags("xor-hash", 42) == ["-O", "xor-hash", "-s", "42"]
    assert build_flags("mod-minimizer-xor-hash", 5, t=10) == ["-O", "mod-minimizer-xor-hash", "--t-length", "10", "-s", "5"]
    assert order_t_values("native") == (None,)
    assert order_t_values("mod-minimizer-greedymini") == (11,)
    assert order_t_values("mod-minimizer-xor-hash") == (10, 11)
    print("ok  orders module: flag mapping (incl. mod-minimizer-* + t)")


def test_lib_orders_sh_converts_hex_seed_to_decimal() -> None:
    """Regression test: sshash's CLI parser reads -s via
    strtoull(value, nullptr, 10) -- BASE 10 ONLY -- so a hex literal like
    config.yaml's "0xDEADBEEF" silently parses as 0 if passed through
    verbatim (strtoull stops at the non-decimal 'x'). lib/orders.sh's
    build_order_flags must convert to decimal ($((seed))) before handing it
    to sshash. This bug previously made every "different seed" in a real
    sweep silently build with seed=0, making all seeds produce byte-identical
    results (caught by comparing two real Cod runs' num_minimizer_positions)."""
    script = (
        f"source '{BENCH_DIR}/lib/orders.sh'\n"
        "build_order_flags xor-hash 0xDEADBEEF ''\n"
        'echo "${ORDER_FLAGS[@]}"\n'
    )
    proc = subprocess.run(["bash", "-c", script], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    assert proc.returncode == 0, f"lib/orders.sh sourcing failed: {proc.stderr}"
    out = proc.stdout.strip()
    assert out == "-s 3735928559 -O xor-hash", (
        f"build_order_flags did not decimal-convert the hex seed: got {out!r} "
        "(0xDEADBEEF = 3735928559 decimal; a literal '0xDEADBEEF' or '0' here "
        "means the strtoull(base=10) bug regressed)"
    )
    print("ok  lib/orders.sh: hex seed 0xDEADBEEF -> decimal 3735928559 for sshash's -s")


def test_sshash_integration() -> None:
    if not SSHASH.is_file():
        print("skip sshash integration (sshash/build/sshash missing; run `make -C sshash all`)")
        return
    if not SAMPLE.is_file():
        print("skip sshash integration (salmonella_enterica sample FASTA missing)")
        return
    if not POS_TXT.is_file() or not NEG_TXT.is_file():
        print("skip sshash integration (committed query sets missing under data/queries/)")
        return

    with tempfile.TemporaryDirectory(prefix="sshash_order_sanity_") as td:
        work = Path(td)
        pos_fa, neg_fa = work / "pos.fa", work / "neg.fa"
        _kmers_fasta(POS_TXT, pos_fa)
        _kmers_fasta(NEG_TXT, neg_fa)
        n_pos = sum(1 for ln in POS_TXT.read_text().splitlines() if ln.strip())
        n_neg = sum(1 for ln in NEG_TXT.read_text().splitlines() if ln.strip())

        counts: dict[str, int] = {}
        for order in ORDERS:
            for t in order_t_values(order):
                label = _run_label(order, t)
                index_path = work / f"{label}.sshash"
                build_json = _build(order, t, index_path)
                num_kmers = int(build_json["num_kmers"])
                counts[label] = num_kmers
                # `check` exercises sshash's own internal is_member/neighbour
                # consistency for every stored k-mer, including its reverse
                # complement -- this is what actually regression-tests
                # mod-minimizer's two-probe Lookup fallback (see
                # sshash/include/util.hpp's compute_mod_minimizer_single_strand
                # doc comment): it caught the "wrong orientation" bug during
                # this feature's own development.
                _check(index_path)

                pos_result = _query(index_path, pos_fa)
                assert int(pos_result["num_positive_kmers"]) == n_pos, (
                    f"{label}: expected all {n_pos} positive k-mers present, "
                    f"got {pos_result['num_positive_kmers']}"
                )
                assert int(pos_result["num_negative_kmers"]) == 0, f"{label}: unexpected absent k-mer in positive set"

                neg_result = _query(index_path, neg_fa)
                assert int(neg_result["num_negative_kmers"]) == n_neg, (
                    f"{label}: expected all {n_neg} negative k-mers absent, "
                    f"got {neg_result['num_negative_kmers']}"
                )
                assert int(neg_result["num_positive_kmers"]) == 0, f"{label}: unexpected present k-mer in negative set"

        assert len(set(counts.values())) == 1, counts
        assert counts["native"] == EXPECTED_NUM_KMERS, counts
        print(f"ok  sshash integration: num_kmers={counts['native']} identical across {sorted(counts)}")
        print(f"     positive/negative query answers identical across all orders/t ({n_pos} pos / {n_neg} neg)")


def main() -> int:
    test_orders_module()
    test_lib_orders_sh_converts_hex_seed_to_decimal()
    test_sshash_integration()
    print("all sanity checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
