"""Offline sanity/equivalence tests for the minimap2 order benchmark harness.
Mirrors benchmark/sshash/tests/test_order_sanity.py's pattern: build under
every (non-greedymini) order on a small local dataset, and assert the orders
don't silently diverge from a shared invariant -- here, that they agree on
the *set of reads mapped* (an order changes which minimizers get sampled,
not whether the underlying alignment problem is solvable), even though
their minimizer counts and exact identity differ.

greedymini / mod-minimizer-greedymini are excluded until a real
orders/GreedyMini/w10_k15_canonical.gm (and a gm_k=5 table for the
mod-minimizer case) is supplied -- see README.md's "Default dataset" and
lib/orders.sh's DEFAULT_ORDERS.

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

ORDERS = ["native", "xor-hash", "xor-pcg-hash", "mod-minimizer-xor-hash", "mod-minimizer-pcg-xor-hash"]


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


def _order_flags(order: str) -> list[str]:
    if order == "native":
        return []
    flags = ["--minimizer-order", order, "--order-seed", "42"]
    if order.startswith("mod-minimizer-"):
        flags += ["--mod-t", "5"]
    return flags


def _build_and_map(dataset: Path, tmp_path: Path, order: str) -> tuple[dict, list[str]]:
    mmi = tmp_path / f"{order}.mmi"
    build = subprocess.run(
        [str(MINIMAP2_BIN), "-d", str(mmi), "-k", "15", "-w", "10", *_order_flags(order), str(dataset / "genome.fa")],
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
def test_all_orders_build_and_map(dataset, tmp_path):
    results = {}
    for order in ORDERS:
        info, mapped = _build_and_map(dataset, tmp_path, order)
        results[order] = (info, mapped)
        assert info["density"] and info["density"] > 0, f"{order}: no minimizers sampled"


@pytest.mark.skipif(not MINIMAP2_BIN.is_file(), reason="minimap2 not built (run `make -C minimap2 all` first)")
def test_mapped_read_set_agrees_across_orders(dataset, tmp_path):
    mapped_sets = {}
    for order in ORDERS:
        _info, mapped = _build_and_map(dataset, tmp_path, order)
        mapped_sets[order] = set(mapped)
    reference = mapped_sets["native"]
    for order, mapped in mapped_sets.items():
        # A handful of reads landing right at the repeat boundary may map
        # under one order and not another -- allow a small slack rather
        # than requiring byte-identical read sets.
        diff = reference.symmetric_difference(mapped)
        assert len(diff) <= max(2, int(0.05 * len(reference))), (
            f"{order} disagrees with native on which reads mapped: {sorted(diff)[:10]}"
        )


@pytest.mark.skipif(not MINIMAP2_BIN.is_file(), reason="minimap2 not built (run `make -C minimap2 all` first)")
def test_non_native_orders_are_denser_or_sparser_but_never_degenerate(dataset, tmp_path):
    n_kmers = 50000  # genome_len passed to make_dataset.py above
    for order in ORDERS:
        info, _mapped = _build_and_map(dataset, tmp_path, order)
        assert 0.01 * n_kmers < info["density"] < n_kmers, (
            f"{order}: density {info['density']} looks degenerate for a {n_kmers}bp genome"
        )
