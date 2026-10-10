#!/usr/bin/env bash
# Loads MINIMAP2_BENCH_ROOT (and optional overrides) from .env, mirroring
# benchmark/sshash/lib/env.sh: never overrides a variable already set in the
# shell, and dies loudly if the root ends up unset (no default -- data/
# results must never land inside the git working tree).

load_bench_env() {
  local here
  here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
  local repo_root
  repo_root="$(cd "$here/../.." && pwd)"
  for f in "$here/.env" "$repo_root/.env"; do
    if [ -f "$f" ]; then
      while IFS='=' read -r k v; do
        [[ "$k" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] || continue
        [ -n "${!k:-}" ] && continue
        export "$k=$v"
      done < <(grep -v '^\s*#' "$f" | grep '=')
    fi
  done
  if [ -z "${MINIMAP2_BENCH_ROOT:-}" ]; then
    echo "error: MINIMAP2_BENCH_ROOT is not set (copy benchmark/minimap2/.env.example to .env and edit it, or export it yourself)" >&2
    exit 1
  fi
  MINIMAP2_DATA_DIR="${MINIMAP2_DATA_DIR:-$MINIMAP2_BENCH_ROOT/data}"
  MINIMAP2_RESULTS_DIR="${MINIMAP2_RESULTS_DIR:-$MINIMAP2_BENCH_ROOT/results}"
  export MINIMAP2_BENCH_ROOT MINIMAP2_DATA_DIR MINIMAP2_RESULTS_DIR
  mkdir -p "$MINIMAP2_DATA_DIR" "$MINIMAP2_RESULTS_DIR/runs" "$MINIMAP2_RESULTS_DIR/csvs"
}
