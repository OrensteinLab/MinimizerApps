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


def data_dir() -> Path:
    rel = _cfg().get("paths", {}).get("data", "data")
    return benchmark_dir() / rel


def default_datasets() -> str:
    return str(_cfg().get("defaults", {}).get("datasets", ""))


def default_k() -> int:
    return int(_cfg().get("defaults", {}).get("k", 31))


def default_m() -> int:
    return int(_cfg().get("defaults", {}).get("m", 21))


def default_orders() -> str:
    return str(
        _cfg().get("defaults", {}).get(
            "orders",
            "native,xor-hash,xor-pcg-hash,greedymini,"
            "mod-minimizer-native-hash,mod-minimizer-xor-hash,"
            "mod-minimizer-pcg-xor-hash,mod-minimizer-greedymini,"
            "mod-minimizer-greedymini-v2",
        )
    )


def default_n_seeds() -> int:
    return int(_cfg().get("defaults", {}).get("n_seeds", 5))


def default_gm_file() -> Path:
    name = str(_cfg().get("defaults", {}).get("gm_file", "w11_k15_canonical.gm"))
    return gm_dir() / name


def mod_t_gm_file() -> Path:
    """mod-minimizer-greedymini's O_t table for t=11 (no CLI override -- see
    lib/orders.sh::order_gm_file)."""
    name = str(_cfg().get("defaults", {}).get("mod_t_gm_file", "w21_k11_canonical.gm"))
    return gm_dir() / name


def mod_t_gm_file_t10() -> Path:
    """mod-minimizer-greedymini's O_t table for t=10."""
    name = str(_cfg().get("defaults", {}).get("mod_t_gm_file_t10", "w22_k10_canonical.gm"))
    return gm_dir() / name


def mod_t_gm_labels() -> dict[str, str]:
    """The labeled GreedyMini table-variant registry (config.yaml's
    top-level mod_t_gm_labels map): label -> filename. Lets scripts.orders'
    ORDERS enumerate registered variants directly, instead of needing a
    matching manual edit alongside every new config.yaml line."""
    return dict(_cfg().get("mod_t_gm_labels", {}) or {})


def mod_t_gm_label_file(label: str) -> Path:
    """gm file for a labeled mod-minimizer-greedymini variant (e.g. "v2" for
    order "mod-minimizer-greedymini-v2"), from config.yaml's top-level
    mod_t_gm_labels map -- register a new variant by adding a line there, not
    by touching this function (mirrors lib/orders.sh's mod_t_gm_label_file)."""
    labels = _cfg().get("mod_t_gm_labels", {}) or {}
    name = labels.get(label)
    if name is None:
        _die(f"unknown GreedyMini label '{label}' (add it to config.yaml's mod_t_gm_labels)")
    return gm_dir() / str(name)
