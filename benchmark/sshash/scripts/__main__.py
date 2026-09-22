#!/usr/bin/env python3
"""Unified CLI for the sshash benchmark Python tools."""

from __future__ import annotations

import sys

COMMANDS = {
    "export": ("scripts.export", "Export runs to CSV/XLSX"),
    "experiment": ("scripts.experiments", "List/expand experiment YAML"),
    "orders": ("scripts.orders", "Order name -> sshash CLI flags"),
    "equivalence": ("scripts.equivalence", "Fixed pos/neg query sets + per-run equivalence checks"),
    "zenodo": ("scripts.zenodo", "Zenodo record 17582116 lookup (maintainers only)"),
    "meta": ("scripts.meta", "Write run meta.json from META_* env"),
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
