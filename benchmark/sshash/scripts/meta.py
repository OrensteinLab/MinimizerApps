"""Write per-run meta.json (reads META_* env vars from sshash_bench.sh)."""

from __future__ import annotations

import json
import os
import socket
import sys
from datetime import datetime, timezone
from pathlib import Path


def main(_argv: list[str] | None = None) -> int:
    cmd_file = Path(os.environ["META_CMD_FILE"])
    command = cmd_file.read_text().splitlines()

    seed_s = os.environ.get("META_SEED", "").strip()
    t_s = os.environ.get("META_T", "").strip()

    meta = {
        "run_id": os.environ["META_RUN_ID"],
        "dataset": os.environ["META_DATASET"],
        "k": int(os.environ["META_K"]),
        "m": int(os.environ["META_M"]),
        "order": os.environ["META_ORDER"],
        "order_seed": int(seed_s, 0) if seed_s else None,
        "seed_idx": int(os.environ.get("META_SEED_IDX", "0")),
        "t": int(t_s) if t_s else None,
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
