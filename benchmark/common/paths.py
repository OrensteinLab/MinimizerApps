"""Directory layout shared by all four benchmark harnesses.

    <ROOT>/data/<dataset>/...      inputs (kraken2 also keeps <ROOT>/db, <ROOT>/zenodo)
    <ROOT>/results/runs/...        per-run outputs (meta.json, time.txt, logs, ...)
    <ROOT>/results/csvs/*.csv      exported tables
    <ROOT>/results/*.xlsx          exported workbooks

<ROOT> comes from <APP>_BENCH_ROOT (KRAKEN2_, KMC3_, SSHASH_, MINIMAP2_). The
optional <APP>_DATA_DIR / <APP>_RESULTS_DIR override the two sub-roots.
There is no extra `benchmark/` folder under <ROOT>.
"""
import os
from pathlib import Path

ENV_PREFIX = {"kraken2": "KRAKEN2", "kmc": "KMC3", "sshash": "SSHASH", "minimap2": "MINIMAP2"}


def bench_paths(app, env=None):
    """Resolve the layout for `app` ('kraken2'|'kmc'|'sshash'|'minimap2') from the environment."""
    env = os.environ if env is None else env
    p = ENV_PREFIX[app]
    root = env.get(p + "_BENCH_ROOT")
    if not root:
        raise SystemExit("error: %s_BENCH_ROOT is not set (see benchmark/%s/.env.example)" % (p, app))
    root = Path(root)
    data = Path(env.get(p + "_DATA_DIR") or root / "data")
    results = Path(env.get(p + "_RESULTS_DIR") or root / "results")
    return {"root": root, "data": data, "results": results,
            "runs": results / "runs", "csvs": results / "csvs", "xlsx": results}


def ensure_dirs(paths):
    for k in ("runs", "csvs"):
        paths[k].mkdir(parents=True, exist_ok=True)
