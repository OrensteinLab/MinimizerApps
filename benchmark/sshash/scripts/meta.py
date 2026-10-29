"""Write per-run meta.json (reads META_* env vars from sshash_bench.sh)."""

from __future__ import annotations

import json
import os
import socket
import sys
from datetime import datetime, timezone
from pathlib import Path

import scheme_names  # benchmark/common (see scripts/__init__.py)


def main(_argv: list[str] | None = None) -> int:
    cmd_file = Path(os.environ["META_CMD_FILE"])
    command = cmd_file.read_text().splitlines()

    seed_s = os.environ.get("META_SEED", "").strip()
    gm_k_s = os.environ.get("META_GM_K", "").strip()
    mod_t_s = os.environ.get("META_MOD_T", "").strip()
    # identity fields (benchmark/common/scheme_names.py): validated, canonical scheme name
    ident = scheme_names.make_id(
        os.environ["META_SCHEME"],
        int(gm_k_s) if gm_k_s else None,
        os.environ.get("META_PROJECTION", "").strip() or None,
        os.environ.get("META_GM_TABLE", "").strip() or None,
        int(mod_t_s) if mod_t_s else None,
        baselines=("native",),
        default_projection=False,
    )

    meta = {
        "run_id": os.environ["META_RUN_ID"],
        "dataset": os.environ["META_DATASET"],
        "k": int(os.environ["META_K"]),
        "m": int(os.environ["META_M"]),
        **ident,
        "modality": os.environ.get("META_MODALITY", "") or "canonical",
        "order_seed": int(seed_s, 0) if seed_s else None,
        "seed_idx": int(os.environ.get("META_SEED_IDX", "0")),
        "gm_file": os.environ.get("META_GM_FILE", "") or "",
        "threads": int(os.environ.get("META_THREADS") or 1),
        "sshash_path": os.environ["META_SSHASH_PATH"],
        "sshash_version": os.environ.get("META_SSHASH_VERSION", ""),
        "host": socket.gethostname(),
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "build_exit_code": int(os.environ.get("META_BUILD_EXIT_CODE", "-1")),
        "check_exit_code": int(os.environ.get("META_CHECK_EXIT_CODE", "-1")),
        "check_ok": os.environ.get("META_CHECK_OK", "0") == "1",
        "bench_exit_code": int(os.environ.get("META_BENCH_EXIT_CODE", "-1")),
        "query_exit_code": int(os.environ.get("META_QUERY_EXIT_CODE", "-1")),
        "query_pos_match": os.environ.get("META_QUERY_POS_MATCH", "") or "",
        "query_neg_match": os.environ.get("META_QUERY_NEG_MATCH", "") or "",
        "index_bytes_on_disk": int(os.environ["META_INDEX_BYTES"]) if os.environ.get("META_INDEX_BYTES") else None,
        "python_wall_s": float(os.environ["META_WALL"]),
        "command": command,
    }
    Path(os.environ["META_PATH"]).write_text(json.dumps(meta, indent=2) + "\n")
    cmd_file.unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
