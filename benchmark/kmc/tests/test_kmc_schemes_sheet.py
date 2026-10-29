#!/usr/bin/env python3
"""The `schemes` sheet: one row per (scheme variant, k, w, L), unified names, no empty density cell."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.density import PARTICULAR_FIELDS, load_particular_rows, particular_csv_path, scan_fastq
from scripts.export import (
    check_schemes_complete,
    schemes_field_order,
    schemes_rows,
    unify_rows,
)
from scripts.expected_density import EXPECTED_SCHEMES_FIELDS, build_dna_keys, expected_schemes_csv_path


def _run(scheme, ds, df, seed=0, *, bins=512, threads=12, projection="", lut=""):
    gm = scheme == "greedymini"
    return {
        "scheme": scheme, "dataset": ds, "k": 28, "p": 9, "w": 20, "n_bins": bins, "threads": threads,
        "strict_memory": "", "seed_idx": seed, "density_factor": df,
        "gm_k": 9 if gm else "", "projection": projection, "gm_table": "w20_k9_canonical" if gm else "",
        "gm_lut": lut,
    }


def _write_expected(root: Path, rows):
    with expected_schemes_csv_path(root).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(EXPECTED_SCHEMES_FIELDS))
        w.writeheader()
        for scheme, proj in rows:
            w.writerow({"scheme": scheme, "projection": proj, "k": 9, "w": 20, "L": 28, "n_seeds": 1,
                        "density_factor_mean": 1.9, "density_factor_std": "" if scheme in ("default", "lex") else 0.01})


def _write_scans(root: Path, rows):
    with particular_csv_path(root).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(PARTICULAR_FIELDS), extrasaction="ignore")
        w.writeheader()
        for scheme, proj, ds, df in rows:
            w.writerow({"scheme": scheme, "projection": proj, "dataset": ds, "k": 9, "w": 20, "L": 28,
                        "density_factor_dna": df})


VARIANTS = [("default", "", ""), ("lex", "", ""), ("xor", "", ""),
            ("greedymini", "at0-cg1", "dna"), ("greedymini", "at0-cg1", "binary"),
            ("greedymini", "at1-cg0", "binary")]


def _setup(tmp_path):
    rows = []
    for scheme, proj, lut in VARIANTS:
        for ds in ("a", "b"):
            for seed in (0, 1):
                # binary variants ran on another bins/threads config only
                cfg = dict(bins=1024, threads=24) if lut == "binary" else {}
                rows.append(_run(scheme, ds, 2.0 + 0.1 * seed, seed, projection=proj, lut=lut, **cfg))
    # lex has no run on dataset "b" -> falls back to the scan rows
    rows = [r for r in rows if not (r["scheme"] == "lex" and r["dataset"] == "b")]
    _write_expected(tmp_path, [("default", ""), ("lex", ""), ("xor", ""),
                               ("greedymini", "at0-cg1"), ("greedymini", "at1-cg0")])
    _write_scans(tmp_path, [("lex", "", "b", 2.5)])
    return rows


def test_one_row_per_variant_including_non_canonical_configs(tmp_path):
    out, datasets = schemes_rows(_setup(tmp_path), tmp_path)
    assert datasets == ["a", "b"]
    assert [(r["scheme"], r["projection"], r["gm_lut"]) for r in out] == [
        ("default", None, None), ("lex", None, None), ("xor", None, None),
        ("greedymini", "at0-cg1", "binary"), ("greedymini", "at0-cg1", "dna"),
        ("greedymini", "at1-cg0", "binary"),
    ]
    assert {(r["k"], r["w"], r["L"]) for r in out} == {(9, 20, 28)}
    check_schemes_complete(out, datasets)  # no empty density cell


def test_values_and_fallbacks(tmp_path):
    out, _ = schemes_rows(_setup(tmp_path), tmp_path)
    by = {(r["scheme"], r["gm_lut"]): r for r in out}
    d = by[("default", None)]
    assert d["a_particular_density_factor_mean"] == 2.05  # mean over the two seeds
    assert d["a_particular_density_factor_std"] == pytest.approx(0.0707, abs=1e-4)
    assert d["expected_density_factor_std"] == 0.0  # deterministic: 0, not blank
    lex = by[("lex", None)]
    assert lex["b_particular_density_factor_mean"] == 2.5 and lex["b_particular_density_factor_std"] == 0.0


def test_fields_have_unified_names_and_no_mod_t(tmp_path):
    fields = schemes_field_order(["a"])
    assert fields[:9] == ["scheme", "k", "w", "L", "gm_k", "projection", "gm_table", "gm_lut", "description"]
    assert "mod_t" not in fields and "is_allowed_denylist" in fields
    assert fields[-4:] == ["expected_density_factor_mean", "expected_density_factor_std",
                           "a_particular_density_factor_mean", "a_particular_density_factor_std"]


def test_missing_density_cell_fails_loudly(tmp_path):
    rows = _setup(tmp_path)
    _write_scans(tmp_path, [])  # drop the lex fallback
    out, datasets = schemes_rows(rows, tmp_path)
    with pytest.raises(SystemExit):
        check_schemes_complete(out, datasets)


def test_unify_rows_swaps_k_and_p():
    assert unify_rows([{"k": 28, "p": 9, "w": 20, "x": 1}]) == [{"L": 28, "k": 9, "w": 20, "x": 1}]


def test_particular_cache_in_old_names_is_still_read(tmp_path):
    with particular_csv_path(tmp_path).open("w") as f:
        f.write("dataset,k,p,w,density_factor_dna\nf_vesca,28,9,20,1.8\n")
    (row,) = load_particular_rows(particular_csv_path(tmp_path))
    assert (row["k"], row["w"], row["L"]) == ("9", "20", "28")


def test_scan_with_explicit_keys_counts_windows(tmp_path):
    fq = tmp_path / "r.fastq"
    seq = "ACGTTGCAAGCTTGCATGCCGTAAGGCTAACGTGCA" * 3
    fq.write_text(f"@r\n{seq}\n+\n{'I' * len(seq)}\n")
    keys = build_dna_keys("lex", 28, 9)
    acc = scan_fastq([fq], p=9, w=20, dna_keys=keys, max_kmers=None)
    assert acc["n_kmers"] == len(seq) - 28 + 1
    assert acc["density_factor_binary"] is None and acc["density_factor_dna"] > 0
