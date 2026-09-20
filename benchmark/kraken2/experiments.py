#!/usr/bin/env python3
"""Load and expand benchmark experiment YAML files.

kraken2's own config field names and shell-variable mapping; the actual
YAML loading/merging/rendering is the shared engine in
benchmark/common/experiments_lib.py (also used by kmc's own
benchmark/kmc/scripts/experiments.py).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

_MONOREPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_MONOREPO_ROOT / "benchmark" / "common"))

import experiments_lib  # noqa: E402

EXPERIMENTS_DIR = Path(__file__).resolve().parent / "experiments"

RUN_KEYS = (
    "scheme",
    "n_seeds",
    "n_classify",
    "threads",
    "parallel_jobs",
    "library",
)

SHELL_VAR_MAPPING = {
    "scheme": "SCHEME",
    "n_seeds": "N_SEEDS",
    "n_classify": "N_CLASSIFY",
    "threads": "SCHEME_THREADS",
    "parallel_jobs": "PARALLEL_JOBS",
    "library": "LIBRARY",
}


def expand_configs(spec: dict[str, Any]) -> list[dict[str, str]]:
    return experiments_lib.expand_configs(spec, RUN_KEYS)


def list_experiments() -> list[tuple[str, str, int]]:
    return experiments_lib.list_experiments(EXPERIMENTS_DIR)


def shell_env_line(cfg: dict[str, str]) -> str:
    return experiments_lib.shell_env_line(cfg, SHELL_VAR_MAPPING)


def shell_env_lines(spec: dict[str, Any]) -> list[str]:
    return experiments_lib.shell_env_lines(spec, RUN_KEYS, SHELL_VAR_MAPPING)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="experiments.py", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    sub.add_parser("list")
    p_expand = sub.add_parser("expand")
    p_expand.add_argument("experiment")
    p_env = sub.add_parser("shell-env")
    p_env.add_argument("experiment")
    p_env.add_argument("index", type=int)
    p_env_all = sub.add_parser("shell-env-all")
    p_env_all.add_argument("experiment")
    p_count = sub.add_parser("count")
    p_count.add_argument("experiment")

    args = ap.parse_args(argv)

    if args.cmd == "list":
        for name, desc, n in list_experiments():
            suffix = f" ({desc})" if desc else ""
            print(f"{name}: {n} config(s){suffix}")
        return 0

    spec = experiments_lib.load_experiment(
        experiments_lib.resolve_path(EXPERIMENTS_DIR, args.experiment)
    )

    if args.cmd == "shell-env-all":
        for line in shell_env_lines(spec):
            print(line)
        return 0

    configs = expand_configs(spec)

    if args.cmd == "count":
        print(len(configs))
        return 0
    if args.cmd == "expand":
        for cfg in configs:
            print(json.dumps(cfg, sort_keys=True))
        return 0
    if args.cmd == "shell-env":
        if args.index < 0 or args.index >= len(configs):
            experiments_lib.die(f"config index out of range: {args.index} (have {len(configs)})")
        print(shell_env_line(configs[args.index]))
        return 0

    experiments_lib.die(f"unknown command: {args.cmd}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
