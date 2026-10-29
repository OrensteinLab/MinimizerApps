#!/usr/bin/env python3
"""Migrate a result set made before the order -> scheme renaming, so the current exporters read it.

    migrate_scheme_names.py APP RESULTS_DIR [--dry-run]

APP is kraken2, minimap2, kmc or sshash. The identity column `order` held a scheme (an order plus the
selection rule: xor, greedymini, mod-minimizer-xor, ...), so it became `scheme`; the apps' flags that
select a scheme were renamed with it. In RESULTS_DIR this:
  * renames the key "order" to "scheme" in every runs/**/meta.json (other keys, e.g. order_seed or kmc's
    signature_order, are kept: they name a true order);
  * sshash: rewrites the recorded build command in meta.json, `-O <scheme>` -> `--scheme <scheme>`;
  * minimap2: rewrites the `[M::main] CMD:` line of build.log / map.log, `--minimizer-order` (or
    the interim `--minimizer-scheme`) -> `--scheme` (experiments/profile_run.sh replays it);
  * renames the CSVs of csvs/ (`*_orders.csv` -> `*_schemes.csv`, `*_density_expected_orders.csv` ->
    `*_density_expected_schemes.csv`, `analysis_accuracy_order.csv` -> `analysis_accuracy_scheme.csv`,
    `order_tests_*.csv` -> `scheme_tests_*.csv`) and renames their header columns (`order` -> `scheme`,
    `n_order` -> `n_scheme`, `n_orders` -> `n_schemes`, `orders_run` -> `schemes_run`, `c_order_only` ->
    `c_scheme_only`); the caches among them (expected and particular density) are inputs of the export.
Re-export afterwards: the workbook and the other tables are rewritten from the migrated files. Running it
again changes nothing. Standard library only.
"""
import argparse
import csv
import io
import json
import os
import re
import sys
from collections import OrderedDict
from pathlib import Path

APPS = ("kraken2", "minimap2", "kmc", "sshash")
COLUMNS = {"order": "scheme", "n_order": "n_scheme", "n_orders": "n_schemes", "orders_run": "schemes_run",
           "c_order_only": "c_scheme_only"}
FILES = [(re.compile(r"^(.*)_density_expected_orders\.csv$"), r"\1_density_expected_schemes.csv"),
         (re.compile(r"^(.*)_orders\.csv$"), r"\1_schemes.csv"),
         (re.compile(r"^analysis_accuracy_order\.csv$"), "analysis_accuracy_scheme.csv"),
         (re.compile(r"^order_tests_(.*)\.csv$"), r"scheme_tests_\1.csv")]
SSHASH_FLAG = re.compile(r"(?<=\s)-O(?=\s)")
MINIMAP2_FLAG = re.compile(r"--minimizer-(?:order|scheme)(?=\s)")
MINIMAP2_CMD = re.compile(r"^(\[M::main\] CMD: .*)$", re.M)


def write(path, text, dry, log):
    log.append(str(path))
    if not dry:
        tmp = path.with_name(path.name + ".migrating")
        tmp.write_text(text)
        os.replace(tmp, path)


def migrate_meta(path, app, dry, log):
    meta = json.loads(path.read_text(), object_pairs_hook=OrderedDict)
    if "order" not in meta and not (app == "sshash" and any(" -O " in c for c in meta.get("command") or [])):
        return
    out = OrderedDict((("scheme" if k == "order" else k), v) for k, v in meta.items())
    if app == "sshash" and isinstance(out.get("command"), list):
        out["command"] = [SSHASH_FLAG.sub("--scheme", c) if " build " in c else c for c in out["command"]]
    write(path, json.dumps(out, indent=2) + "\n", dry, log)


def migrate_minimap2_log(path, dry, log):
    text = path.read_text(errors="surrogateescape")
    new = MINIMAP2_CMD.sub(lambda m: MINIMAP2_FLAG.sub("--scheme", m.group(1)), text)
    if new != text:
        write(path, new, dry, log)


def migrate_csv(path, dry, log):
    text = path.read_text()
    rows = list(csv.reader(io.StringIO(text)))
    if rows and any(c in COLUMNS for c in rows[0]):
        rows[0] = [COLUMNS.get(c, c) for c in rows[0]]
        buf = io.StringIO()
        csv.writer(buf, lineterminator="\n").writerows(rows)
        write(path, buf.getvalue(), dry, log)
    for pat, repl in FILES:
        if pat.match(path.name):
            dest = path.with_name(pat.sub(repl, path.name))
            log.append("%s -> %s" % (path, dest.name))
            if not dry:
                if dest.exists():
                    dest.unlink()   # a stale copy the export would overwrite anyway
                os.replace(path, dest)
            break


def migrate(app, results, dry=False):
    log = []
    for meta in sorted(results.glob("runs/**/meta.json")):
        migrate_meta(meta, app, dry, log)
    if app == "minimap2":
        for name in ("build.log", "map.log"):
            for p in sorted(results.glob("runs/**/" + name)):
                migrate_minimap2_log(p, dry, log)
    for p in sorted(results.glob("**/csvs/*.csv")) + sorted(results.glob("*.csv")):
        migrate_csv(p, dry, log)
    return log


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("app", choices=APPS)
    ap.add_argument("results", type=Path)
    ap.add_argument("--dry-run", action="store_true", help="list what would change, change nothing")
    a = ap.parse_args(argv)
    if not (a.results / "runs").is_dir() and not (a.results / "csvs").is_dir():
        ap.error("%s is not a result set (no runs/ or csvs/)" % a.results)
    log = migrate(a.app, a.results, a.dry_run)
    for line in log:
        print(("would migrate " if a.dry_run else "migrated ") + line)
    print("%s %s: %d changes%s" % (a.app, a.results, len(log), " (dry run)" if a.dry_run else ""), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
