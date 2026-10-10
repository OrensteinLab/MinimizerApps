"""Shared helpers for the sshash benchmark."""

from __future__ import annotations

import os
import sys
from pathlib import Path

from scripts.config import benchmark_dir, repo_root


def load_env_file(path: Path) -> None:
    if not path.is_file():
        return
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, val = line.partition("=")
        key, val = key.strip(), val.strip()
        if len(val) >= 2 and val[0] == val[-1] and val[0] in "\"'":
            val = val[1:-1]
        if key and key not in os.environ:
            os.environ[key] = val


def load_dotenv() -> None:
    load_env_file(benchmark_dir() / ".env")
    load_env_file(repo_root().parent / ".env")


def die(msg: str, code: int = 1) -> None:
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(code)


def bench_paths() -> tuple[Path, Path]:
    """(data dir, results dir): $ROOT/data and $ROOT/results unless overridden."""
    from paths import bench_paths as _layout  # benchmark/common (on sys.path via scripts/__init__)

    if not os.environ.get("SSHASH_BENCH_ROOT", "").strip():
        die("SSHASH_BENCH_ROOT is not set (see benchmark/sshash/.env)")
    p = _layout("sshash")
    return p["data"], p["results"]
