"""Load and expand benchmark experiment YAML files.

kmc's own config field names and shell-variable mapping; the actual YAML
loading/merging/rendering is the shared engine in
benchmark/common/experiments_lib.py (also used by kraken2's own
benchmark/kraken2/experiments.py).
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
    "input",
    "threads",
    "mem",
    "bins",
    "strict_mem",
    "signature_order",
    "order_label",
    "custom_order",
    "n_seeds",
    "signature_len",
    "force",
    "run_tag",
    "gm_lut_mode",
)

SHELL_VAR_MAPPING = {
    "datasets": "DATASETS",
    "k": "KS",
    "input": "INPUTS",
    "threads": "THREADS",
    "mem": "MEM_OVERRIDE",
    "bins": "BINS",
    "strict_mem": "STRICT_MEM",
    "signature_order": "SIGNATURE_ORDER",
    "order_label": "ORDER_LABEL",
    "custom_order": "CUSTOM_ORDER",
    "n_seeds": "N_SEEDS",
    "signature_len": "SIGNATURE_LEN",
    "force": "FORCE",
    "run_tag": "RUN_TAG",
    "gm_lut_mode": "GM_LUT_MODE",
}

BOOL_VARS = ("STRICT_MEM", "FORCE")


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
