#!/usr/bin/env python3
"""Offline sanity checks for sshash's custom minimizer schemes (no network,
no big downloads -- uses the in-repo sshash sample dataset only).

Mirrors benchmark/kmc/tests/test_custom_order_sanity.py::test_kmc_integration:
build the same input under every scheme, assert they all report the exact
same num_kmers, and assert every scheme's index answers the exact same fixed
positive/negative k-mer queries (see data/queries/) identically.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))

from scripts.config import monorepo_root, repo_root
from scripts.schemes import SCHEMES, build_flags, check_gm_scheme, expand_variants, scheme_dir_name

BENCH_DIR = Path(__file__).resolve().parents[1]
SSHASH = Path(os.environ.get("SSHASH_BIN") or repo_root() / "build" / "sshash")
SAMPLE = repo_root() / "data" / "unitigs_stitched" / "salmonella_enterica_k31_ust.fa.gz"
POS_TXT = BENCH_DIR / "data" / "queries" / "salmonella_enterica_k31_pos.txt"
NEG_TXT = BENCH_DIR / "data" / "queries" / "salmonella_enterica_k31_neg.txt"
# k=31, m=21 gives w=11, the window the GreedyMini tables were built for.
K, M = 31, 21
EXPECTED_NUM_KMERS = 4787534  # independent of m: m only affects minimizer construction, not the k-mer set

SEEDS = {
    "native": 1,
    "xor": 42,
    "xor-pcg-hash": 7,
    "greedymini": 98,
    "mod-minimizer-xor": 102,
    "mod-minimizer-greedymini": 103,
}


def _variants() -> list[dict]:
    """Every variant config.yaml registers, plus the flipped projection of each GreedyMini table."""
    out = expand_variants(SCHEMES, K, M)
    out += [v for v in expand_variants(["greedymini"], K, M, projections=["at1-cg0"])]
    return out


def _run_label(v: dict, regular: bool) -> str:
    label = scheme_dir_name(v["scheme"], v["gm_k"], v["projection"], v["mod_t"]).replace("/", "_")
    return f"{label}_regular" if regular else label


def _run(cmd: list[str]) -> tuple[int, str]:
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    return proc.returncode, (proc.stdout or "") + (proc.stderr or "")


def _json_blob(text: str) -> dict:
    m = re.search(r"\{.*\}", text, re.DOTALL)
    assert m, f"no JSON blob found in output:\n{text[-2000:]}"
    return json.loads(m.group(0))


def _build(v: dict, regular: bool, out_path: Path) -> dict:
    scheme = v["scheme"]
    flags = build_flags(scheme, SEEDS[scheme], v["gm_file"], K, M, v["projection"]) + (["--regular"] if regular else [])
    cmd = [str(SSHASH), "build", "-i", str(SAMPLE), "-k", str(K), "-m", str(M), "-t", "4", "--verbose", *flags, "-o", str(out_path)]
    rc, out = _run(cmd)
    assert rc == 0, f"build failed for scheme={_run_label(v, regular)}: {out[-2000:]}"
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


def test_schemes_module() -> None:
    assert SCHEMES == ("native", "xor", "xor-pcg-hash", "greedymini", "mod-minimizer-xor",
                      "mod-minimizer-greedymini")
    variants = expand_variants(SCHEMES, K, M)
    # one variant per scheme, except greedymini/mod-minimizer-greedymini (one per table x projection)
    assert [(v["scheme"], v["gm_k"], v["projection"], v["gm_table"], v["mod_t"]) for v in variants] == [
        ("native", None, None, None, None), ("xor", None, None, None, None),
        ("xor-pcg-hash", None, None, None, None),
        ("greedymini", 13, "at0-cg1", "w11_k13_canonical", None),
        ("greedymini", 13, "at1-cg0", "w11_k13_canonical", None),
        ("greedymini", 15, "at0-cg1", "w11_k15_canonical", None),
        ("greedymini", 15, "at1-cg0", "w11_k15_canonical", None),
        ("mod-minimizer-xor", None, None, None, 10),
        ("mod-minimizer-greedymini", 10, "at0-cg1", "w22_k10_canonical", 10),
        ("mod-minimizer-greedymini", 10, "at1-cg0", "w22_k10_canonical", 10),
    ]
    gm13, gm15 = variants[3], variants[5]
    flags = build_flags("greedymini", 7, gm13["gm_file"], K, M, "at1-cg0")
    assert flags[:2] == ["--scheme", "greedymini"] and flags[-2:] == ["-s", "7"]
    assert "--gm-projection" in flags and flags[flags.index("--gm-projection") + 1] == "at1-cg0"
    assert "--gm-projection" not in build_flags("greedymini", 7, gm15["gm_file"])
    assert check_gm_scheme("greedymini", 21, gm15["gm_file"]) is None
    assert check_gm_scheme("greedymini", 21, gm13["gm_file"]) is None
    assert "m - gm_k even" in (check_gm_scheme("greedymini", 20, gm13["gm_file"]) or "")
    assert scheme_dir_name("greedymini", 15, "at1-cg0", None) == "greedymini/k15_at1-cg0"
    assert scheme_dir_name("xor", None, None, None) == "xor"
    # strand-symmetric mod-minimizer: t = 4 + (m - 4) mod w for xor, the table's gm_k otherwise
    assert build_flags("mod-minimizer-xor", 7) == ["--scheme", "mod-minimizer-xor", "--mod-t", "10", "-s", "7"]
    gm10 = build_flags("mod-minimizer-greedymini", 7)
    assert gm10[:2] == ["--scheme", "mod-minimizer-greedymini"] and gm10[gm10.index("--mod-t") + 1] == "10"
    assert build_flags("mod-minimizer-xor", 7, k=31, m=20)[3] == "8"  # w = 12
    for old in ("xor-hash", "mod-xor-hash", "greedymini-k15"):
        try:
            build_flags(old, 7)
        except ValueError:
            continue
        raise AssertionError(f"old scheme spelling {old} accepted")
    print("ok  schemes module: flag mapping")


def test_lib_schemes_sh_converts_hex_seed_to_decimal() -> None:
    """Regression test: sshash's CLI parser reads -s via
    strtoull(value, nullptr, 10) -- BASE 10 ONLY -- so a hex literal like
    config.yaml's "0xDEADBEEF" silently parses as 0 if passed through
    verbatim (strtoull stops at the non-decimal 'x'). lib/schemes.sh's
    build_scheme_flags must convert to decimal ($((seed))) before handing it
    to sshash. This bug previously made every "different seed" in a real
    sweep silently build with seed=0, making all seeds produce byte-identical
    results (caught by comparing two real Cod runs' num_minimizer_positions)."""
    script = (
        f"source '{BENCH_DIR}/lib/schemes.sh'\n"
        "build_scheme_flags xor 0xDEADBEEF ''\n"
        'echo "${SCHEME_FLAGS[@]}"\n'
    )
    proc = subprocess.run(["bash", "-c", script], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    assert proc.returncode == 0, f"lib/schemes.sh sourcing failed: {proc.stderr}"
    out = proc.stdout.strip()
    assert out == "-s 3735928559 --scheme xor", (
        f"build_scheme_flags did not decimal-convert the hex seed: got {out!r} "
        "(0xDEADBEEF = 3735928559 decimal; a literal '0xDEADBEEF' or '0' here "
        "means the strtoull(base=10) bug regressed)"
    )
    print("ok  lib/schemes.sh: hex seed 0xDEADBEEF -> decimal 3735928559 for sshash's -s")


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

    with tempfile.TemporaryDirectory(prefix="sshash_scheme_sanity_") as td:
        work = Path(td)
        pos_fa, neg_fa = work / "pos.fa", work / "neg.fa"
        _kmers_fasta(POS_TXT, pos_fa)
        _kmers_fasta(NEG_TXT, neg_fa)
        n_pos = sum(1 for ln in POS_TXT.read_text().splitlines() if ln.strip())
        n_neg = sum(1 for ln in NEG_TXT.read_text().splitlines() if ln.strip())

        counts: dict[str, int] = {}
        # every scheme with canonical parsing, plus native under --regular (the two-probe
        # "regular parsing" modality)
        native = expand_variants(["native"], K, M)[0]
        for v, regular in [(v, False) for v in _variants()] + [(native, True)]:
            label = _run_label(v, regular)
            index_path = work / f"{label}.sshash"
            build_json = _build(v, regular, index_path)
            num_kmers = int(build_json["num_kmers"])
            counts[label] = num_kmers
            # `check` exercises sshash's own internal is_member/neighbour
            # consistency for every stored k-mer, including its reverse
            # complement (and, under --regular, the two-probe Lookup fallback).
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
        print(f"     positive/negative query answers identical across all schemes ({n_pos} pos / {n_neg} neg)")


def main() -> int:
    test_schemes_module()
    test_lib_schemes_sh_converts_hex_seed_to_decimal()
    test_sshash_integration()
    print("all sanity checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
