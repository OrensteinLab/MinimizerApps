"""Load shared benchmark config (single source for seeds, paths, defaults)."""

from __future__ import annotations

import os
import sys
from functools import lru_cache
from pathlib import Path
from typing import Any

from yaml_util import load_yaml  # benchmark/common (see scripts/__init__.py)

_CONFIG_NAME = "config.yaml"


def _die(msg: str, code: int = 1) -> None:
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(code)


def benchmark_dir() -> Path:
    env = os.environ.get("BENCHMARK_DIR", "").strip()
    if env:
        return Path(env).resolve()
    cwd = Path.cwd()
    if (cwd / _CONFIG_NAME).is_file():
        return cwd.resolve()
    root = Path(__file__).resolve().parent.parent
    if (root / _CONFIG_NAME).is_file():
        return root
    _die("BENCHMARK_DIR is not set and config.yaml was not found")


def repo_root() -> Path:
    # The kmc source tree. Derived from this module's own location (not
    # benchmark_dir(), which can be redirected via $BENCHMARK_DIR to an
    # unrelated data/results path): benchmark/kmc/scripts/ -> repo root -> kmc/.
    return Path(__file__).resolve().parents[3] / "kmc"


@lru_cache(maxsize=1)
def _cfg() -> dict[str, Any]:
    path = benchmark_dir() / _CONFIG_NAME
    if not path.is_file():
        _die(f"missing {path}")
    return load_yaml(path.read_text())


def seeds() -> tuple[int, ...]:
    raw = _cfg().get("seeds")
    if not isinstance(raw, list) or not raw:
        _die("config.yaml: seeds must be a non-empty list")
    out: list[int] = []
    for i, item in enumerate(raw):
        try:
            out.append(int(item, 0) if isinstance(item, str) else int(item))
        except (TypeError, ValueError):
            _die(f"config.yaml: seeds[{i}] is not an integer: {item!r}")
    return tuple(out)


def gm_dir() -> Path:
    rel = _cfg().get("paths", {}).get("gm_dir", "minimizers/GreedyMini")
    return benchmark_dir() / rel


def experiments_dir() -> Path:
    rel = _cfg().get("paths", {}).get("experiments", "experiments")
    return benchmark_dir() / rel


def data_dir() -> Path:
    rel = _cfg().get("paths", {}).get("data", "data")
    return benchmark_dir() / rel


def default_datasets() -> str:
    return str(_cfg().get("defaults", {}).get("datasets", ""))


def default_signature_len() -> int:
    return int(_cfg().get("defaults", {}).get("signature_len", 9))
