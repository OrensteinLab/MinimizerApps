#!/usr/bin/env python3
"""Unified CLI for KMC3 benchmark Python tools."""

from __future__ import annotations

import sys

COMMANDS = {
    "export": ("scripts.export", "Export runs to CSV/XLSX"),
    "experiment": ("scripts.experiments", "List/expand experiment YAML"),
    "gm": ("scripts.gm", "Build/ensure GreedyMini tables"),
    "density": ("scripts.density", "Expected / particular GreedyMini density"),
    "expected-density": (
        "scripts.expected_density",
        "Theoretical expected density_factor per scheme (any scheme)",
    ),
    "meta": ("scripts.meta", "Write run meta.json from META_* env"),
    "make-readlen-dataset": (
        "scripts.make_readlen_dataset",
        "Build a joined-read derived dataset (read-length confound check)",
    ),
}


def usage() -> None:
    print("usage: python -m scripts <command> [args...]")
    print("\ncommands:")
    for name, (_, desc) in COMMANDS.items():
        print(f"  {name:<12} {desc}")


def main() -> int:
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        usage()
        return 0

    cmd = sys.argv[1]
    if cmd not in COMMANDS:
        print(f"error: unknown command: {cmd}", file=sys.stderr)
        usage()
        return 1

    module_name, _ = COMMANDS[cmd]
    mod = __import__(module_name, fromlist=["main"])
    return int(mod.main(sys.argv[2:]))


if __name__ == "__main__":
    raise SystemExit(main())
