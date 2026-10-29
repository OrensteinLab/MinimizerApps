"""Expected and particular GreedyMini density factors.

Expected density (binary alphabet) is the fraction of charged (w+p)-windows,
equal to lim |f(S)| / (|S|-p+1) on a random binary string. The density factor
is that value times (w+1) (Marçais et al.; GreedyMini).

Particular density is the same ratio on a concrete string. For DNA we use
KMC's current ranking (AT/GC projection + GM ranks + DNA tie-break). For
binary we project the same bases to bits and apply the binary GM order.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import random
import struct
import subprocess
import sys
import tempfile
from collections import deque
from pathlib import Path

from paths import bench_paths
from scripts.common import die, load_dotenv
from scripts.config import gm_dir
from xlsx_util import write_csv as _shared_write_csv
from scripts.gm import (
    GM_SERIALIZER_BIN,
    GM_SERIALIZER_DIR,
    GM_SIGMA_EXTEND_BIN,
    GM_SIGMA_EXTEND_DIR,
    gm_binary_path,
    gm_table_stem,
    make,
)
from scripts.metrics import minimizer_window
from scripts.schemes import resolve_scheme_label
from density import density_factor as _shared_density_factor

_BASE_FROM_BYTE = [-1] * 256
for _ch, _v in zip(b"ACGTNacgtn", (0, 1, 2, 3, -1, 0, 1, 2, 3, -1)):
    _BASE_FROM_BYTE[_ch] = _v


def density_factor(n_selected: int, n_kmers: int, w: int) -> float | None:
    # 4-decimal rounding is this file's own established precision (kept
    # exactly, unlike benchmark/common/density.py's 6-decimal default, so
    # existing CSV output/tests here are unaffected).
    return _shared_density_factor(n_selected, n_kmers, w, decimals=4)


def _read_raw_gm(path: Path) -> list[int]:
    """Read the classic [u64 len][len x u64] .gm format (see
    common/minimizer-schemes/greedymini_rank.h's load_gm_file_raw)."""
    data = path.read_bytes()
    (n,) = struct.unpack_from("<Q", data, 0)
    return list(struct.unpack_from(f"<{n}Q", data, 8))


def _write_raw_gm(path: Path, ranks: list[int]) -> None:
    with path.open("wb") as f:
        f.write(struct.pack("<Q", len(ranks)))
        f.write(struct.pack(f"<{len(ranks)}Q", *ranks))


def load_gm_ranks(path: Path, p: int | None = None) -> tuple[int, list[int]]:
    """Return (binary k, rank[0..2^k)). Accepts classic .gm or a rank-list .txt
    (serialized via the shared C++ gm_serializer, not reimplemented here)."""
    if path.suffix == ".gm":
        ranks = _read_raw_gm(path)
        k = (len(ranks) - 1).bit_length()
        return k, ranks
    make(GM_SERIALIZER_DIR)
    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "table.gm"
        args = [str(GM_SERIALIZER_BIN), "--canonical", str(path), "-o", str(out)]
        if p is not None:
            args += ["-k", str(p)]
        subprocess.run(args, check=True)
        ranks = _read_raw_gm(out)
    k = (len(ranks) - 1).bit_length()
    return k, ranks


def resolve_gm_ranks(k_kmc: int, signature_len: int, table_dir: Path | None = None) -> tuple[int, list[int], Path]:
    root = table_dir or gm_dir()
    binary = gm_binary_path(k_kmc, signature_len, root)
    txt = root / f"{gm_table_stem(k_kmc, signature_len)}.txt"
    if binary.is_file():
        p, ranks = load_gm_ranks(binary)
        return p, ranks, binary
    if txt.is_file():
        p, ranks = load_gm_ranks(txt, signature_len)
        return p, ranks, txt
    die(f"GM table missing: {binary} (no source {txt})")
    return 0, [], binary  # unreachable; die() exits


def _cmp_keys_binary(ranks: list[int], p: int) -> list[int]:
    del p
    inf = len(ranks)
    return [r if r < inf else inf for r in ranks]


def _cmp_keys_dna(ranks: list[int], p: int, seed: int = 0) -> list[int]:
    """DNA-space sort key from the shared C++ gm_sigma_extend tool (not
    reimplemented here) -- write ranks as a plain .gm, extend it, read back
    the sort_key array from the resulting _dna.gm."""
    make(GM_SIGMA_EXTEND_DIR)
    with tempfile.TemporaryDirectory() as td:
        binary = Path(td) / "table.gm"
        dna = Path(td) / "table_dna.gm"
        _write_raw_gm(binary, ranks)
        subprocess.run(
            [str(GM_SIGMA_EXTEND_BIN), "--seed", str(seed), "-o", str(dna), str(binary)],
            check=True,
        )
        data = dna.read_bytes()
        (_magic, _k, _seed, n) = struct.unpack_from("<4Q", data, 0)
        return list(struct.unpack_from(f"<{n}I", data, 32))


def _sliding_runs(keys: list[int], w: int) -> tuple[int, int, int]:
    """Return (n_windows, n_signature_runs, n_selected_positions).

    n_windows is the KMC k-mer count (one window of w signatures per k-mer).
    n_signature_runs matches KMC super-k-mers (runs of the same min key).
    n_selected_positions is |f(S)| (argmin index changes; paper particular).

    NOT the same function as benchmark/common/density.py's `sliding_runs`
    (used by sshash's and minimap2's own scripts/density.py): this one pops
    the deque on a *strict* `>` (rightmost tie wins), the shared one on `>=`
    (leftmost tie wins) -- a real tie-break difference, not just a naming
    one, discovered while extracting the shared module; left as its own
    KMC-specific copy rather than risk changing KMC's existing behavior.
    """
    n_win = len(keys) - w + 1
    if n_win <= 0:
        return 0, 0, 0
    dq: deque[int] = deque()
    n_super = 0
    n_sel = 0
    prev_key: int | None = None
    prev_pos: int | None = None
    for i, key in enumerate(keys):
        while dq and keys[dq[-1]] > key:
            dq.pop()
        dq.append(i)
        if dq[0] <= i - w:
            dq.popleft()
        if i + 1 < w:
            continue
        pos = dq[0]
        min_key = keys[pos]
        if prev_key is None or min_key != prev_key:
            n_super += 1
            prev_key = min_key
        if prev_pos != pos:
            n_sel += 1
            prev_pos = pos
    return n_win, n_super, n_sel


def context_is_charged(bits: int, ranks: list[int], p: int, w: int) -> bool:
    """True if sliding one step in a binary (w+p)-window selects a new position.

    A (w+p)-window holds two overlapping w-windows of p-mers. The fraction of
    such contexts where the leftmost-min index changes equals expected density
    (Marçais et al.; GreedyMini).
    """
    length = w + p
    mask = (1 << p) - 1
    kmers = [0] * (w + 1)
    for i in range(w + 1):
        kmers[i] = (bits >> (length - p - i)) & mask
    inf = len(ranks)

    def argmin(lo: int, hi: int) -> int:
        best_i = lo
        best_r = ranks[kmers[lo]] if ranks[kmers[lo]] < inf else inf
        for i in range(lo + 1, hi):
            r = ranks[kmers[i]] if ranks[kmers[i]] < inf else inf
            if r < best_r:
                best_r = r
                best_i = i
        return best_i

    return argmin(0, w) != argmin(1, w + 1)


def expected_binary_density(
    ranks: list[int],
    p: int,
    w: int,
    n_samples: int = 4_000_000,
    seed: int = 0xC0FFEE01,
) -> dict[str, float | int]:
    """Monte Carlo expected density of the binary GM order (charged windows)."""
    rng = random.Random(seed)
    length = w + p
    mask = (1 << length) - 1
    charged = 0
    for _ in range(n_samples):
        bits = rng.getrandbits(length) & mask
        if context_is_charged(bits, ranks, p, w):
            charged += 1
    dens = charged / n_samples
    return {
        "n_samples": n_samples,
        "n_charged": charged,
        "minimizer_density": round(dens, 8),
        "density_factor": round(dens * (w + 1), 4),
        "w": w,
        "p": p,
        "alphabet": "binary",
    }


def expected_binary_density_stream(
    ranks: list[int],
    p: int,
    w: int,
    n_bits: int = 8_000_000,
    seed: int = 0xC0FFEE01,
) -> dict[str, float | int]:
    """Particular density on one long random binary string (cross-check)."""
    rng = random.Random(seed)
    keys_cmp = _cmp_keys_binary(ranks, p)
    mask = (1 << p) - 1
    state = rng.getrandbits(p) & mask
    keys = [keys_cmp[state]]
    for _ in range(n_bits - p):
        state = ((state << 1) | rng.getrandbits(1)) & mask
        keys.append(keys_cmp[state])
    n_win, n_super, n_sel = _sliding_runs(keys, w)
    return {
        "n_bits": n_bits,
        "n_kmers": n_win,
        "n_super_kmers": n_super,
        "n_selected": n_sel,
        "density_factor": density_factor(n_super, n_win, w),
        "density_factor_positions": density_factor(n_sel, n_win, w),
        "w": w,
        "p": p,
        "alphabet": "binary",
    }


def expected_dna_density(
    ranks: list[int],
    p: int,
    w: int,
    n_bases: int = 4_000_000,
    seed: int = 0xC0FFEE01,
    tie_seed: int = 0,
) -> dict[str, float | int]:
    """Particular density of the KMC DNA ranking on a random DNA string."""
    rng = random.Random(seed)
    dna_keys = _cmp_keys_dna(ranks, p, tie_seed)
    mask = (1 << (2 * p)) - 1
    state = rng.getrandbits(2 * p) & mask
    keys = [dna_keys[state]]
    for _ in range(n_bases - p):
        state = ((state << 2) | rng.getrandbits(2)) & mask
        keys.append(dna_keys[state])
    n_win, n_super, n_sel = _sliding_runs(keys, w)
    return {
        "n_bases": n_bases,
        "n_kmers": n_win,
        "n_super_kmers": n_super,
        "n_selected": n_sel,
        "density_factor": density_factor(n_super, n_win, w),
        "density_factor_positions": density_factor(n_sel, n_win, w),
        "w": w,
        "p": p,
        "alphabet": "dna",
    }


def expected_at_read_length(
    ranks: list[int],
    p: int,
    w: int,
    read_len: int,
    n_reads: int = 20_000,
    seed: int = 0xC0FFEE01,
    alphabet: str = "binary",
    tie_seed: int = 0,
) -> dict[str, float | int | None]:
    """Expected particular DF on random reads of a fixed length (finite-string effect)."""
    rng = random.Random(seed)
    if alphabet == "binary":
        key_lut = _cmp_keys_binary(ranks, p)
        bits_per_base = 1
        pack = lambda s, b: _pack_bin(s, b, p)
        rand_base = lambda: rng.getrandbits(1)
        n_sym = 1 << p
    else:
        key_lut = _cmp_keys_dna(ranks, p, tie_seed)
        bits_per_base = 2
        pack = lambda s, b: _pack_dna(s, b, p)
        rand_base = lambda: rng.getrandbits(2)
        n_sym = 1 << (2 * p)
    del n_sym, bits_per_base
    n_kmers = n_super = 0
    for _ in range(n_reads):
        bases = [rand_base() for _ in range(read_len)]
        keys = _segment_keys(bases, p, pack, key_lut)
        nw, ns, _ = _sliding_runs(keys, w)
        n_kmers += nw
        n_super += ns
    return {
        "read_len": read_len,
        "n_reads": n_reads,
        "n_kmers": n_kmers,
        "n_super_kmers": n_super,
        "density_factor": density_factor(n_super, n_kmers, w),
        "alphabet": alphabet,
        "w": w,
        "p": p,
    }


def _segment_keys(bases: list[int], p: int, pack, key_lut: list[int]) -> list[int]:
    if len(bases) < p:
        return []
    state = 0
    for i in range(p):
        state = pack(state, bases[i])
    keys = [key_lut[state]]
    for b in bases[p:]:
        state = pack(state, b)
        keys.append(key_lut[state])
    return keys


def _pack_dna(state: int, base: int, p: int) -> int:
    return ((state << 2) | base) & ((1 << (2 * p)) - 1)


def _pack_bin(state: int, base: int, p: int) -> int:
    bit = 1 if base in (1, 2) else 0
    return ((state << 1) | bit) & ((1 << p) - 1)


def particular_on_bases(
    bases: list[int],
    *,
    p: int,
    w: int,
    dna_keys: list[int],
    bin_keys: list[int] | None,
) -> dict[str, int]:
    dna = _segment_keys(bases, p, lambda s, b: _pack_dna(s, b, p), dna_keys)
    n_win_d, n_super_d, n_sel_d = _sliding_runs(dna, w)
    if bin_keys is None:  # order with no binary-alphabet form (lex, xor, ...)
        n_win_b = n_super_b = n_sel_b = 0
    else:
        binary = _segment_keys(bases, p, lambda s, b: _pack_bin(s, b, p), bin_keys)
        n_win_b, n_super_b, n_sel_b = _sliding_runs(binary, w)
    gc = sum(1 for b in bases if b in (1, 2))
    return {
        "n_bases": len(bases),
        "n_gc": gc,
        "n_kmers": n_win_d,
        "n_super_dna": n_super_d,
        "n_sel_dna": n_sel_d,
        "n_super_bin": n_super_b,
        "n_sel_bin": n_sel_b,
        "n_kmers_bin": n_win_b,
    }


def _flush_segment(
    acc: dict[str, int],
    bases: list[int],
    p: int,
    w: int,
    dna_keys: list[int],
    bin_keys: list[int] | None,
) -> None:
    if len(bases) < p + w - 1:
        bases.clear()
        return
    part = particular_on_bases(bases, p=p, w=w, dna_keys=dna_keys, bin_keys=bin_keys)
    for key, val in part.items():
        acc[key] = acc.get(key, 0) + val
    bases.clear()


def scan_fastq(
    paths: list[Path],
    *,
    p: int,
    w: int,
    ranks: list[int] | None = None,
    max_kmers: int | None = 50_000_000,
    tie_seed: int = 0,
    dna_keys: list[int] | None = None,
) -> dict:
    """Density-only scan (no kmc counting). Default: the GM table `ranks` (DNA
    and binary comparators). Pass `dna_keys` (any order's 4^p comparator, e.g.
    expected_density.build_dna_keys) to scan that order instead; the binary
    columns are then 0/None."""
    if dna_keys is None:
        if ranks is None:
            raise ValueError("scan_fastq needs ranks or dna_keys")
        dna_keys = _cmp_keys_dna(ranks, p, tie_seed)
    bin_keys = _cmp_keys_binary(ranks, p) if ranks is not None else None
    acc: dict[str, int] = {
        "n_bases": 0,
        "n_gc": 0,
        "n_kmers": 0,
        "n_super_dna": 0,
        "n_sel_dna": 0,
        "n_super_bin": 0,
        "n_sel_bin": 0,
        "n_kmers_bin": 0,
        "n_reads": 0,
        "n_segments": 0,
    }
    bases: list[int] = []
    for path in paths:
        opener = gzip.open if path.suffix == ".gz" or path.name.endswith(".gz") else open
        with opener(path, "rt", errors="replace") as fh:
            for line_i, line in enumerate(fh):
                if line_i & 3 != 1:
                    continue
                acc["n_reads"] += 1
                for ch in line.encode("ascii", "replace"):
                    b = _BASE_FROM_BYTE[ch]
                    if b < 0:
                        if bases:
                            acc["n_segments"] += 1
                            _flush_segment(acc, bases, p, w, dna_keys, bin_keys)
                        continue
                    bases.append(b)
                if bases:
                    acc["n_segments"] += 1
                    _flush_segment(acc, bases, p, w, dna_keys, bin_keys)
                if max_kmers is not None and acc["n_kmers"] >= max_kmers:
                    break
        if max_kmers is not None and acc["n_kmers"] >= max_kmers:
            break
    n_win = acc["n_kmers"]
    n_win_b = acc["n_kmers_bin"] or n_win
    acc["gc_fraction"] = round(acc["n_gc"] / acc["n_bases"], 4) if acc["n_bases"] else None
    acc["density_factor_dna"] = density_factor(acc["n_super_dna"], n_win, w)
    acc["density_factor_binary"] = density_factor(acc["n_super_bin"], n_win_b, w) if bin_keys is not None else None
    acc["density_factor_dna_positions"] = density_factor(acc["n_sel_dna"], n_win, w)
    acc["density_factor_binary_positions"] = (
        density_factor(acc["n_sel_bin"], n_win_b, w) if bin_keys is not None else None
    )
    acc["w"] = w
    acc["p"] = p
    return acc


def dataset_fastq_files(data_root: Path, dataset: str) -> list[Path]:
    lst_candidates = (
        data_root / dataset / "files_gz.lst",
    )
    for lst in lst_candidates:
        if lst.is_file():
            out = []
            for line in lst.read_text().splitlines():
                line = line.strip()
                if line:
                    out.append(Path(line))
            if out:
                return out
    dir_candidates = (
        data_root / dataset / "fastq_gz",
        data_root / dataset,
    )
    out: list[Path] = []
    for d in dir_candidates:
        if not d.is_dir():
            continue
        out.extend(sorted(d.glob("*.fastq.gz")))
        out.extend(sorted(d.glob("*.fq.gz")))
        out.extend(sorted(d.glob("*.fastq")))
        if out:
            return out
    return []


def particular_csv_path(results: Path) -> Path:
    return results / "kmc3_density_particular.csv"


def expected_csv_path(results: Path) -> Path:
    return results / "kmc3_density_expected.csv"


# Rows with a blank `scheme` are the legacy per-dataset GreedyMini DNA-vs-binary
# diagnostic (keyed by dataset/k/p). Rows with `scheme` (+ `projection` for
# greedymini) are per-variant density-only scans (`scripts density
# particular-schemes`) that fill the schemes sheet's particular-density cells for
# (variant, dataset) pairs with no kmc run.
PARTICULAR_FIELDS = (
    "scheme",
    "projection",
    "dataset",
    "k",
    "w",
    "L",
    "n_bases",
    "n_kmers",
    "n_reads",
    "gc_fraction",
    "density_factor_dna",
    "density_factor_binary",
    "density_factor_dna_positions",
    "density_factor_binary_positions",
    "n_super_dna",
    "n_super_bin",
    "n_sel_dna",
    "n_sel_bin",
)


EXPECTED_FIELDS = (
    "k",
    "w",
    "L",
    "alphabet",
    "density_factor",
    "density_factor_positions",
    "minimizer_density",
    "n_samples",
    "n_kmers",
    "method",
)


def unify_param_row(row: dict) -> dict:
    """Rows in the unified names (k = minimizer length, w, L = k + w - 1). A cache
    written before the rename has k = k-mer length and p = minimizer length: map it."""
    if "L" in row or "p" not in row:
        return row
    out = {("L" if f == "k" else "k" if f == "p" else f): v for f, v in row.items()}
    return out


def load_particular_rows(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    return [unify_param_row(r) for r in csv.DictReader(path.open())]


def particular_index(rows: list[dict]) -> dict[tuple, dict]:
    out: dict[tuple, dict] = {}
    for r in rows:
        if r.get("scheme"):  # per-variant row, see particular_variant_index
            continue
        try:
            key = (r["dataset"], int(r["L"]), int(r["k"]))  # (dataset, k-mer length L, minimizer length k)
        except (KeyError, TypeError, ValueError):
            continue
        out[key] = r
    return out


def particular_variant_index(rows: list[dict]) -> dict[tuple, dict]:
    """Per-variant scan rows keyed (scheme, projection, dataset, L, k)."""
    out: dict[tuple, dict] = {}
    for r in rows:
        if not r.get("scheme"):
            continue
        try:
            key = (r["scheme"], r.get("projection") or "", r["dataset"], int(r["L"]), int(r["k"]))
        except (KeyError, TypeError, ValueError):
            continue
        out[key] = r
    return out


def _float_or_none(v) -> float | None:
    if v is None or v == "":
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def attach_density_refs(rows: list[dict], *, particular: dict[tuple, dict]) -> None:
    """Attach binary particular DF on greedymini rows with the dna LUT (gm_lut=dna) (compare measured DNA vs binary scan)."""
    for row in rows:
        row["binary_density_factor"] = None
        if resolve_scheme_label(row) != "greedymini" or (row.get("gm_lut") or "") != "dna":
            continue
        ds = row.get("dataset")
        try:
            k = int(row.get("k") or 0)
            p = int(row.get("p") or row.get("signature_len") or 0)
        except (TypeError, ValueError):
            continue
        part = particular.get((ds, k, p)) if ds else None
        if not part:
            continue
        row["binary_density_factor"] = _float_or_none(
            part.get("density_factor_binary_positions") or part.get("density_factor_binary")
        )


def _write_csv(path: Path, rows: list[dict], fields: tuple[str, ...]) -> None:
    _shared_write_csv(path, list(fields), rows)
    print(f"wrote {path} ({len(rows)} rows)")


def _bench_roots() -> tuple[Path, Path]:
    """(data dir, csvs dir) from the shared layout (see benchmark/common/paths.py)."""
    load_dotenv()
    p = bench_paths("kmc")
    return p["data"], p["csvs"]


def cmd_expected(args: argparse.Namespace) -> int:
    p = args.p
    w = minimizer_window(args.k, p)
    kk, ranks, src = resolve_gm_ranks(args.k, p, Path(args.gm_dir) if args.gm_dir else None)
    if kk != p:
        die(f"table k={kk} does not match -p {p} ({src})")
    charged = expected_binary_density(ranks, p, w, n_samples=args.samples, seed=args.seed)
    stream = expected_binary_density_stream(ranks, p, w, n_bits=args.stream_bases, seed=args.seed)
    dna = expected_dna_density(ranks, p, w, n_bases=args.stream_bases, seed=args.seed, tie_seed=args.tie_seed)
    print(
        f"table {src.name}  p={p} w={w}  UHS={sum(1 for r in ranks if r < (1 << p))}",
        file=sys.stderr,
    )
    print(
        f"expected binary (charged windows, n={charged['n_samples']}): "
        f"df={charged['density_factor']}"
    )
    print(
        f"random binary stream (n_bits={stream['n_bits']}): "
        f"df_super={stream['density_factor']}  df_pos={stream['density_factor_positions']}"
    )
    print(
        f"random DNA stream, KMC ranking (n_bases={dna['n_bases']}): "
        f"df_super={dna['density_factor']}  df_pos={dna['density_factor_positions']}"
    )
    for L in (64, 100, 351):
        b = expected_at_read_length(ranks, p, w, L, n_reads=8000, seed=args.seed, alphabet="binary")
        d = expected_at_read_length(ranks, p, w, L, n_reads=8000, seed=args.seed, alphabet="dna", tie_seed=args.tie_seed)
        print(f"random reads L={L}: binary_df={b['density_factor']}  dna_df={d['density_factor']}")
    if args.output:
        rows = [
            {
                "k": p,
                "w": w,
                "L": args.k,
                "alphabet": "binary",
                "density_factor": charged["density_factor"],
                "density_factor_positions": stream["density_factor_positions"],
                "minimizer_density": charged["minimizer_density"],
                "n_samples": charged["n_samples"],
                "n_kmers": stream["n_kmers"],
                "method": "charged_windows",
            },
            {
                "k": p,
                "w": w,
                "L": args.k,
                "alphabet": "dna",
                "density_factor": dna["density_factor"],
                "density_factor_positions": dna["density_factor_positions"],
                "minimizer_density": None,
                "n_samples": dna["n_kmers"],
                "n_kmers": dna["n_kmers"],
                "method": "random_dna_kmc_rank",
            },
        ]
        for L in (64, 100, 351):
            b = expected_at_read_length(ranks, p, w, L, n_reads=8000, seed=args.seed, alphabet="binary")
            d = expected_at_read_length(ranks, p, w, L, n_reads=8000, seed=args.seed, alphabet="dna", tie_seed=args.tie_seed)
            rows.append(
                {
                    "k": p,
                    "w": w,
                    "L": args.k,
                    "alphabet": "binary",
                    "density_factor": b["density_factor"],
                    "density_factor_positions": None,
                    "minimizer_density": None,
                    "n_samples": b["n_reads"],
                    "n_kmers": b["n_kmers"],
                    "method": f"random_reads_len{L}",
                }
            )
            rows.append(
                {
                    "k": p,
                    "w": w,
                    "L": args.k,
                    "alphabet": "dna",
                    "density_factor": d["density_factor"],
                    "density_factor_positions": None,
                    "minimizer_density": None,
                    "n_samples": d["n_reads"],
                    "n_kmers": d["n_kmers"],
                    "method": f"random_reads_len{L}",
                }
            )
        _write_csv(Path(args.output), rows, EXPECTED_FIELDS)
    return 0


def cmd_particular(args: argparse.Namespace) -> int:
    data_root, results = _bench_roots()
    p = args.p
    w = minimizer_window(args.k, p)
    _kk, ranks, src = resolve_gm_ranks(args.k, p, Path(args.gm_dir) if args.gm_dir else None)
    datasets = [d.strip() for d in args.datasets.split(",") if d.strip()]
    out_rows: list[dict] = []
    for ds in datasets:
        files = dataset_fastq_files(data_root, ds)
        if not files:
            print(f"warning: no FASTQ for {ds} under {data_root}", file=sys.stderr)
            continue
        print(f"scanning {ds} ({len(files)} files, max_kmers={args.max_kmers}) table={src.name}", file=sys.stderr)
        acc = scan_fastq(
            files,
            p=p,
            w=w,
            ranks=ranks,
            max_kmers=args.max_kmers,
            tie_seed=args.tie_seed,
        )
        row = {**acc, "dataset": ds, "k": p, "w": w, "L": args.k}
        out_rows.append(row)
        print(
            f"{ds}: dna_df={acc['density_factor_dna']}  "
            f"bin_df={acc['density_factor_binary']}  "
            f"gc={acc['gc_fraction']}  kmers={acc['n_kmers']}",
            file=sys.stderr,
        )
    out = Path(args.output) if args.output else particular_csv_path(results)
    # keep the per-variant scan rows (particular-schemes) -- only the legacy rows are rebuilt here
    keep = [r for r in load_particular_rows(out) if r.get("scheme")]
    _write_csv(out, out_rows + keep, PARTICULAR_FIELDS)
    return 0


def cmd_particular_schemes(args: argparse.Namespace) -> int:
    """Density-only scan (no kmc run) of one scheme variant on each dataset;
    upserts rows keyed (scheme, projection, dataset, k, p) into
    kmc3_density_particular.csv, leaving every other row untouched."""
    from scripts.expected_density import PROJECTION_AT0, build_dna_keys

    data_root, results = _bench_roots()
    p = args.p
    w = minimizer_window(args.k, p)
    projection = (args.projection or PROJECTION_AT0) if args.scheme == "greedymini" else ""
    from scripts.config import seeds as bench_seeds

    tie_seed = args.tie_seed if args.tie_seed is not None else bench_seeds()[0]
    keys = build_dna_keys(args.scheme, args.k, p, seed=None if args.scheme in ("default", "lex") else tie_seed, projection=projection or PROJECTION_AT0,
                          gm_dir=Path(args.gm_dir) if args.gm_dir else None)
    out = Path(args.output) if args.output else particular_csv_path(results)
    rows = load_particular_rows(out)
    index = {}
    for i, r in enumerate(rows):
        index[(r.get("scheme") or "", r.get("projection") or "", r.get("dataset"), r.get("L"), r.get("k"))] = i
    for ds in [d.strip() for d in args.datasets.split(",") if d.strip()]:
        files = dataset_fastq_files(data_root, ds)
        if not files:
            die(f"no FASTQ for {ds} under {data_root}")
        print(f"scanning {ds} ({len(files)} files, max_kmers={args.max_kmers}) "
              f"scheme={args.scheme} {projection}", file=sys.stderr)
        acc = scan_fastq(files, p=p, w=w, dna_keys=keys, max_kmers=args.max_kmers)
        row = {"scheme": args.scheme, "projection": projection, "dataset": ds, **acc, "k": p, "w": w, "L": args.k}
        print(f"{ds}: dna_df={acc['density_factor_dna']} gc={acc['gc_fraction']} kmers={acc['n_kmers']}",
              file=sys.stderr)
        key = (args.scheme, projection, ds, str(args.k), str(p))
        if key in index:
            rows[index[key]] = row
        else:
            index[key] = len(rows)
            rows.append(row)
        _write_csv(out, rows, PARTICULAR_FIELDS)  # after every dataset: scans are slow
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="scripts density", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_exp = sub.add_parser("expected", help="expected binary (and random-DNA) density factor")
    p_exp.add_argument("--k", type=int, default=28)
    p_exp.add_argument("-p", type=int, default=9)
    p_exp.add_argument("--samples", type=int, default=4_000_000)
    p_exp.add_argument("--stream-bases", type=int, default=4_000_000)
    p_exp.add_argument("--seed", type=lambda x: int(x, 0), default=0xC0FFEE01)
    p_exp.add_argument("--tie-seed", type=lambda x: int(x, 0), default=0)
    p_exp.add_argument("--gm-dir", type=str, default=None)
    p_exp.add_argument("-o", "--output", type=str, default=None)

    p_part = sub.add_parser("particular", help="particular DNA vs binary density on FASTQ")
    p_part.add_argument("--k", type=int, default=28)
    p_part.add_argument("-p", type=int, default=9)
    p_part.add_argument("--datasets", default="f_vesca,g_gallus,m_balbisiana,h_sapiens_2,h_sapiens_3")
    p_part.add_argument("--max-kmers", type=int, default=50_000_000)
    p_part.add_argument("--tie-seed", type=lambda x: int(x, 0), default=0)
    p_part.add_argument("--gm-dir", type=str, default=None)
    p_part.add_argument("-o", "--output", type=str, default=None)

    p_sch = sub.add_parser(
        "particular-schemes",
        help="density-only scan of one scheme variant (no kmc run); upserts into kmc3_density_particular.csv",
    )
    p_sch.add_argument("--scheme", required=True, choices=("default", "lex", "xor", "xor-pcg-hash", "greedymini"))
    p_sch.add_argument("--projection", default=None, help="greedymini only: at0-cg1 (default) | at1-cg0")
    p_sch.add_argument("--k", type=int, default=28)
    p_sch.add_argument("-p", type=int, default=9)
    p_sch.add_argument("--datasets", required=True)
    p_sch.add_argument("--max-kmers", type=int, default=12_000_000)
    p_sch.add_argument("--tie-seed", type=lambda x: int(x, 0), default=None,
                       help="xor/xor-pcg-hash/greedymini tie-break seed (default: config seeds[0])")
    p_sch.add_argument("--gm-dir", type=str, default=None)
    p_sch.add_argument("-o", "--output", type=str, default=None)

    args = ap.parse_args(argv)
    if args.cmd == "particular-schemes":
        return cmd_particular_schemes(args)
    if args.cmd == "expected":
        return cmd_expected(args)
    if args.cmd == "particular":
        return cmd_particular(args)
    die(f"unknown command: {args.cmd}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
