"""Write per-run meta.json (reads META_* env vars from kmc3_bench.sh)."""

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

    sig_len = os.environ.get("META_SIGNATURE_LEN", "").strip()
    seed_s = os.environ.get("META_SEED", "").strip()
    order = os.environ["META_ORDER"]

    meta = {
        "run_id": os.environ["META_RUN_ID"],
        "dataset": os.environ["META_DATASET"],
        "label": os.environ["META_LABEL"],
        "k": int(os.environ["META_K"]),
        "input_format": os.environ["META_KIND"],
        "threads": int(os.environ["META_THREADS"]),
        "mem_limit_gb": int(os.environ["META_MEM"]),
        "n_bins": int(os.environ.get("META_BINS") or 512),
        "strict_memory": os.environ.get("META_STRICT_MEM", "0") == "1",
        "config_tag": os.environ.get("META_CONFIG_TAG", ""),
        "signature_len": int(sig_len) if sig_len else None,
        "signature_order": "default" if order == "default" else "custom",
        "signature_order_mode": os.environ.get("META_SIG_MODE", ""),
        "gm_lut_mode": os.environ.get("META_GM_LUT_MODE", ""),
        "signature_order_label": os.environ.get("META_ORDER_LABEL", "") if order != "default" else "",
        "signature_order_seed": int(seed_s, 0) if seed_s else None,
        "seed_idx": int(os.environ.get("META_SEED_IDX", "0")),
        "custom_order_file": os.environ.get("META_CUSTOM_ORDER", "") if order == "gm" else "",
        "kmc_path": os.environ["META_KMC_PATH"],
        "kmc_version": os.environ["META_KMC_VERSION"],
        "host": socket.gethostname(),
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "exit_code": int(os.environ["META_EXIT_CODE"]),
        "python_wall_s": float(os.environ["META_WALL"]),
        "command": command,
    }
    Path(os.environ["META_PATH"]).write_text(json.dumps(meta, indent=2) + "\n")
    cmd_file.unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
