"""Offline sanity/equivalence tests for the minimap2 scheme benchmark harness.
Mirrors benchmark/sshash/tests/test_scheme_sanity.py's pattern: build under
every scheme on a small local dataset, and assert the schemes don't silently
diverge from a shared invariant -- here, that they agree on the *set of
reads mapped* (a scheme changes which minimizers get sampled, not whether
the underlying alignment problem is solvable), even though their minimizer
counts and exact identity differ.

greedymini uses orders/GreedyMini/w10_k15_canonical.gm (gm_k=15, matching
this test's -k 15 -w 10), in both AT/CG projections (--gm-projection
at0-cg1 / at1-cg0, a parameter of greedymini rather than a scheme of its own).
mod-minimizer-greedymini uses w20_k9_canonical.gm, a table of 9-mers, so it
runs at -k 19 --mod-t 9 (the benchmark's map-pb setting).

Run with: python3.11 -m pytest benchmark/minimap2/tests -q
(needs a python3 new enough for these type hints -- see lib/constants.sh's
own PYTHON fallback logic; system python3 on this host is 3.6).
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent.parent
REPO_ROOT = HERE.parent.parent
MINIMAP2_BIN = REPO_ROOT / "minimap2" / "minimap2"
DATA_DIR = HERE / "_test_data"

GM_FILE = REPO_ROOT / "orders" / "GreedyMini" / "w10_k15_canonical.gm"
MOD_GM_FILE = REPO_ROOT / "orders" / "GreedyMini" / "w20_k9_canonical.gm"
# (scheme, greedymini projection): a variant label is "<scheme>" or "<scheme>_<projection>"
SCHEMES = ["native", "xor", "xor-pcg-hash", "greedymini_at0-cg1", "greedymini_at1-cg0", "mod-minimizer-xor",
           "mod-minimizer-greedymini_at0-cg1", "mod-minimizer-greedymini_at1-cg0"]


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    out = tmp_path_factory.mktemp("minimap2_bench_data")
    subprocess.run(
        [sys.executable, str(HERE / "lib" / "make_dataset.py"),
         "--out-dir", str(out), "--seed", "1",
         "--genome-len", "50000", "--repeat-len", "4000", "--n-reads", "60"],
        check=True, cwd=HERE,
    )
    return out


def _scheme_flags(label: str) -> list[str]:
    scheme, _, projection = label.partition("_")
    if scheme == "native":
        return []
    if scheme == "greedymini":
        return ["--scheme", "greedymini", "--greedymini-file", str(GM_FILE), "--gm-projection", projection]
    if scheme == "mod-minimizer-greedymini":
        return ["--scheme", scheme, "--greedymini-file", str(MOD_GM_FILE), "--gm-projection", projection,
                "--order-seed", "42", "-k", "19", "--mod-t", "9"]
    flags = ["--scheme", scheme, "--order-seed", "42"]
    if scheme.startswith("mod-minimizer-"):
        flags += ["--mod-t", "5"]
    return flags


def _build_and_map(dataset: Path, tmp_path: Path, scheme: str) -> tuple[dict, list[str]]:
    mmi = tmp_path / f"{scheme}.mmi"
    build = subprocess.run(
        [str(MINIMAP2_BIN), "-d", str(mmi), "-k", "15", "-w", "10", *_scheme_flags(scheme), str(dataset / "genome.fa")],  # a later -k wins
        capture_output=True, text=True,
    )
    assert build.returncode == 0, build.stderr
    density = None
    for line in build.stderr.splitlines():
        if "distinct minimizers" in line:
            density = int(line.split("distinct minimizers:")[1].split()[0].rstrip(","))
    paf = subprocess.run(
        [str(MINIMAP2_BIN), str(mmi), str(dataset / "reads.fa")],
        capture_output=True, text=True,
    )
    assert paf.returncode == 0, paf.stderr
    mapped = sorted({line.split("\t")[0] for line in paf.stdout.splitlines() if line.strip()})
    return {"density": density}, mapped


@pytest.mark.skipif(not MINIMAP2_BIN.is_file(), reason="minimap2 not built (run `make -C minimap2 all` first)")
def test_all_schemes_build_and_map(dataset, tmp_path):
    results = {}
    for scheme in SCHEMES:
        info, mapped = _build_and_map(dataset, tmp_path, scheme)
        results[scheme] = (info, mapped)
        assert info["density"] and info["density"] > 0, f"{scheme}: no minimizers sampled"


@pytest.mark.skipif(not MINIMAP2_BIN.is_file(), reason="minimap2 not built (run `make -C minimap2 all` first)")
def test_mapped_read_set_agrees_across_schemes(dataset, tmp_path):
    mapped_sets = {}
    for scheme in SCHEMES:
        _info, mapped = _build_and_map(dataset, tmp_path, scheme)
        mapped_sets[scheme] = set(mapped)
    reference = mapped_sets["native"]
    for scheme, mapped in mapped_sets.items():
        # A handful of reads landing right at the repeat boundary may map
        # under one scheme and not another -- allow a small slack rather
        # than requiring byte-identical read sets.
        diff = reference.symmetric_difference(mapped)
        assert len(diff) <= max(2, int(0.05 * len(reference))), (
            f"{scheme} disagrees with native on which reads mapped: {sorted(diff)[:10]}"
        )


@pytest.mark.skipif(not MINIMAP2_BIN.is_file(), reason="minimap2 not built (run `make -C minimap2 all` first)")
def test_non_native_schemes_are_denser_or_sparser_but_never_degenerate(dataset, tmp_path):
    n_kmers = 50000  # genome_len passed to make_dataset.py above
    for scheme in SCHEMES:
        info, _mapped = _build_and_map(dataset, tmp_path, scheme)
        assert 0.01 * n_kmers < info["density"] < n_kmers, (
            f"{scheme}: density {info['density']} looks degenerate for a {n_kmers}bp genome"
        )
