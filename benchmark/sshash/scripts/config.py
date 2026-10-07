"""Load shared benchmark config (single source for seeds, paths, dataset list, GreedyMini labels)."""

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
    # The sshash source tree. Derived from this module's own location (not
    # benchmark_dir(), which can be redirected via $BENCHMARK_DIR to an
    # unrelated data/results path): benchmark/sshash/scripts/ -> repo root -> sshash/.
    return Path(__file__).resolve().parents[3] / "sshash"


def monorepo_root() -> Path:
    return Path(__file__).resolve().parents[3]


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
    rel = _cfg().get("paths", {}).get("gm_dir", "../../orders/GreedyMini")
    return benchmark_dir() / rel


def experiments_dir() -> Path:
    rel = _cfg().get("paths", {}).get("experiments", "experiments")
    return benchmark_dir() / rel


def greedymini_labels() -> dict[str, str]:
    """Labeled plain-greedymini table registry (config.yaml's top-level
    greedymini_labels map): label -> filename, for orders "greedymini-<label>"
    (and their "-at1" flipped-polarity twins)."""
    return dict(_cfg().get("greedymini_labels", {}) or {})


def mod_greedymini_labels() -> dict[str, str]:
    """config.yaml's mod_greedymini_labels: label -> GreedyMini table of t-mers, for orders
    "mod-minimizer-greedymini-<label>"."""
    return dict(_cfg().get("mod_greedymini_labels", {}) or {})


def mod_greedymini_label_file(label: str) -> Path:
    name = mod_greedymini_labels().get(label)
    if name is None:
        _die(f"unknown mod-minimizer GreedyMini label '{label}' (add it to config.yaml's mod_greedymini_labels)")
    return gm_dir() / str(name)


def greedymini_label_file(label: str) -> Path:
    labels = _cfg().get("greedymini_labels", {}) or {}
    name = labels.get(label)
    if name is None:
        _die(f"unknown GreedyMini label '{label}' (add it to config.yaml's greedymini_labels)")
    return gm_dir() / str(name)
