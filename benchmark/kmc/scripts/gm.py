"""GreedyMini table paths and build/ensure logic.

The actual GreedyMini serialization/extension logic is C++ (shared with
kraken2 where possible), not reimplemented in Python: this module only
resolves paths and runs the compiled tools (built via their Makefiles).
"""

from __future__ import annotations

import argparse
import struct
import subprocess
import sys
from pathlib import Path

import order_names
from scripts.common import die
from scripts.config import default_signature_len, gm_dir

KMC_DNA_GM_MAGIC = 0x314D444E414D4B43  # "KMCANDM1"; see custom_order.cpp

_MONOREPO_ROOT = Path(__file__).resolve().parents[3]
GM_SERIALIZER_DIR = _MONOREPO_ROOT / "common/minimizer-orders"
GM_SERIALIZER_BIN = GM_SERIALIZER_DIR / "gm_serializer"
GM_SIGMA_EXTEND_DIR = _MONOREPO_ROOT / "benchmark/kmc/tools"
GM_SIGMA_EXTEND_BIN = GM_SIGMA_EXTEND_DIR / "gm_sigma_extend"


def make(build_dir: Path) -> None:
    # make's output goes to stderr: callers (lib/gm.sh's resolve_gm_path) capture this script's stdout
    # as the table path, so any build chatter on stdout would end up inside the path
    subprocess.run(["make", "-C", str(build_dir)], check=True, stdout=sys.stderr)


def signature_len_effective(signature_len: int | None = None) -> int:
    if signature_len is not None:
        return signature_len
    import os

    env = os.environ.get("SIGNATURE_LEN", "").strip()
    return int(env) if env else default_signature_len()


def gm_window(k: int, signature_len: int | None = None) -> int:
    p = signature_len_effective(signature_len)
    return k - p + 1


def gm_table_stem(k: int, signature_len: int | None = None) -> str:
    p = signature_len_effective(signature_len)
    w = gm_window(k, p)
    return f"w{w}_k{p}_canonical"


def gm_label(k: int, signature_len: int | None = None) -> str:
    return f"gm_{gm_table_stem(k, signature_len)}"


def gm_binary_path(
    k: int,
    signature_len: int | None = None,
    table_dir: Path | None = None,
) -> Path:
    root = table_dir or gm_dir()
    return root / f"{gm_table_stem(k, signature_len)}.gm"


def txt_path(k: int, signature_len: int | None = None, table_dir: Path | None = None) -> Path:
    root = table_dir or gm_dir()
    return root / f"{gm_table_stem(k, signature_len)}.txt"


def dna_gm_path(
    k: int,
    signature_len: int | None = None,
    seed: int | None = None,
    table_dir: Path | None = None,
) -> Path:
    root = table_dir or gm_dir()
    stem = root / gm_table_stem(k, signature_len)
    if seed is None:
        return Path(f"{stem}_dna.gm")
    return Path(f"{stem}_dna_s{seed:#x}.gm")


def resolve_gm_table(
    k: int,
    *,
    lut_mode: str = "dna",
    signature_len: int | None = None,
    seed: int | None = None,
    table_dir: Path | None = None,
) -> Path:
    """Return GreedyMini table path for the requested LUT backend."""
    mode = lut_mode.strip().lower()
    if mode == "binary":
        binary = gm_binary_path(k, signature_len, table_dir)
        if not binary.is_file():
            ensure_gm_tables(k, signature_len, seed, table_dir)
        if not binary.is_file():
            die(f"GM binary table missing: {binary}")
        return binary.resolve()
    return Path(ensure_gm_tables(k, signature_len, seed, table_dir)).resolve()


def read_dna_gm_seed(path: Path) -> int | None:
    if not path.is_file() or path.stat().st_size < 24:
        return None
    magic, _k, seed, _n = struct.unpack_from("<4Q", path.read_bytes(), 0)
    if magic != KMC_DNA_GM_MAGIC:
        return None
    return int(seed)


def is_dna_gm(path: Path) -> bool:
    return read_dna_gm_seed(path) is not None


def _needs_rebuild(out: Path, *deps: Path) -> bool:
    if not out.is_file():
        return True
    out_m = out.stat().st_mtime
    return any(dep.is_file() and dep.stat().st_mtime > out_m for dep in deps)


def ensure_gm_tables(
    k: int,
    signature_len: int | None = None,
    seed: int | None = None,
    table_dir: Path | None = None,
) -> Path:
    """Ensure binary .gm and seed-specific _dna.gm exist; return DNA table path."""
    p = signature_len_effective(signature_len)
    root = table_dir or gm_dir()
    txt = txt_path(k, p, root)
    binary = gm_binary_path(k, p, root)
    dna = dna_gm_path(k, p, seed, root)

    if txt.is_file() and _needs_rebuild(binary, txt):
        print(f"building {binary.name} from {txt.name}", file=sys.stderr)
        make(GM_SERIALIZER_DIR)
        subprocess.run(
            [str(GM_SERIALIZER_BIN), "--canonical", str(txt), "-o", str(binary),
             "-k", str(p)],
            check=True,
        )
    elif not binary.is_file():
        die(f"GM binary table missing: {binary} (no source {txt})")

    legacy = dna_gm_path(k, p, None, root)
    if seed is not None:
        if read_dna_gm_seed(legacy) == seed and legacy.is_file():
            return legacy
        if _needs_rebuild(dna, binary):
            print(f"building {dna.name} (seed {seed:#x})", file=sys.stderr)
            make(GM_SIGMA_EXTEND_DIR)
            subprocess.run(
                [str(GM_SIGMA_EXTEND_BIN), "--seed", str(seed), "-o", str(dna),
                 str(binary)],
                check=True,
            )
        return dna

    if _needs_rebuild(legacy, binary):
        use_seed = read_dna_gm_seed(legacy) or 0
        print(f"building {legacy.name} (seed {use_seed:#x})", file=sys.stderr)
        make(GM_SIGMA_EXTEND_DIR)
        subprocess.run(
            [str(GM_SIGMA_EXTEND_BIN), "--seed", str(use_seed), "-o", str(legacy),
             str(binary)],
            check=True,
        )
    return legacy


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="scripts gm", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_ensure = sub.add_parser("ensure", help="build if needed; print DNA .gm path")
    p_ensure.add_argument("--k", type=int, required=True)
    p_ensure.add_argument("-p", type=int, default=None)
    p_ensure.add_argument("--seed", type=lambda x: int(x, 0), default=None)
    p_ensure.add_argument("--gm-dir", type=Path, default=None)

    p_dna = sub.add_parser("dna-seed", help="print embedded seed for _dna.gm (exit 1 if not DNA table)")
    p_dna.add_argument("path", type=Path)

    p_resolve = sub.add_parser("resolve", help="print table path for dna or binary LUT mode")
    p_resolve.add_argument("--k", type=int, required=True)
    p_resolve.add_argument("-p", type=int, default=None)
    p_resolve.add_argument("--seed", type=lambda x: int(x, 0), default=None)
    p_resolve.add_argument("--lut-mode", choices=("dna", "binary"), default="dna")
    p_resolve.add_argument("--gm-dir", type=Path, default=None)

    p_dir = sub.add_parser("order-dir", help="print the run-dir path for an order (order_names.order_dir_name)")
    p_dir.add_argument("--order", required=True)
    p_dir.add_argument("--gm-k", type=int, default=None)
    p_dir.add_argument("--projection", default=None)

    args = ap.parse_args(argv)
    if args.cmd == "order-dir":
        order_names.validate_order(args.order, ("default", "lex"))
        print(order_names.order_dir_name(args.order, args.gm_k, args.projection))
        return 0
    if args.cmd == "ensure":
        print(ensure_gm_tables(args.k, args.p, args.seed, args.gm_dir))
        return 0
    if args.cmd == "dna-seed":
        s = read_dna_gm_seed(args.path)
        if s is None:
            return 1
        print(f"{s:#x}")
        return 0
    if args.cmd == "resolve":
        print(resolve_gm_table(args.k, lut_mode=args.lut_mode, signature_len=args.p, seed=args.seed, table_dir=args.gm_dir))
        return 0
    die(f"unknown command: {args.cmd}")
    return 1
