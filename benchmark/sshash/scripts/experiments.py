"""Load and expand benchmark experiment YAML files.

sshash's own config field names and shell-variable mapping; the actual YAML
loading/merging/rendering is the shared engine in
benchmark/common/experiments_lib.py (also used by kraken2's and kmc's own
experiments.py).
"""

from __future__ import annotations

import argparse
import json
from typing import Any

import experiments_lib
from scripts.common import die
from scripts.config import experiments_dir

RUN_KEYS = (
    "datasets",
    "k",
    "m",
    "schemes",
    "gm_tables",
    "gm_projections",
    "n_seeds",
    "threads",
    "force",
    "run_tag",
)

SHELL_VAR_MAPPING = {
    "datasets": "DATASETS",
    "k": "K",
    "m": "M",
    "schemes": "SCHEMES",
    "gm_tables": "GM_TABLES",
    "gm_projections": "GM_PROJECTIONS",
    "n_seeds": "N_SEEDS",
    "threads": "THREADS",
    "force": "FORCE",
    "run_tag": "RUN_TAG",
}

BOOL_VARS = ("FORCE",)


def expand_configs(spec: dict[str, Any]) -> list[dict[str, str]]:
    return experiments_lib.expand_configs(spec, RUN_KEYS)


def list_experiments() -> list[tuple[str, str, int]]:
    return experiments_lib.list_experiments(experiments_dir())


def shell_env_line(cfg: dict[str, str]) -> str:
    return experiments_lib.shell_env_line(cfg, SHELL_VAR_MAPPING, BOOL_VARS)


def shell_env_lines(spec: dict[str, Any]) -> list[str]:
    return experiments_lib.shell_env_lines(spec, RUN_KEYS, SHELL_VAR_MAPPING, BOOL_VARS)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="scripts experiment", description=__doc__)
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
        experiments_lib.resolve_path(experiments_dir(), args.experiment)
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
            die(f"config index out of range: {args.index} (have {len(configs)})")
        print(shell_env_line(configs[args.index]))
        return 0

    die(f"unknown command: {args.cmd}")
    return 1
