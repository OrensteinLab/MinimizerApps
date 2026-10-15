#!/usr/bin/env bash
# Source this file, then call `wait_low_load` before each timed step.
#
# The benchmarks are timing measurements on a shared host: block until the
# machine is quiet. "Quiet" = fewer than BENCH_LOAD_MAX busy cores (default 4)
# measured from /proc/stat over BENCH_LOAD_WINDOW seconds (default 5), on
# BENCH_LOAD_STABLE consecutive samples (default 1). The call is made while the
# benchmark itself is idle, so what is measured is other people's load -- unlike
# the 1-minute load average, which lags and would count the previous run.
# Prints the busy-core figure it started under (so the run log records it).
# Disabled when BENCH_WAIT_LOW_LOAD=0.
_bench_busy_cores() {  # busy cores over a $1-second window
  local w="${1:-5}" a b c d e f g h u1 t1 u2 t2
  read -r _ a b c d e f g h _ < /proc/stat; t1=$((a+b+c+d+e+f+g+h)); u1=$((t1-d-e))
  sleep "$w"
  read -r _ a b c d e f g h _ < /proc/stat; t2=$((a+b+c+d+e+f+g+h)); u2=$((t2-d-e))
  awk -v u=$((u2-u1)) -v t=$((t2-t1)) -v n="$(nproc)" 'BEGIN{printf "%.1f", (t > 0 ? u/t*n : 0)}'
}

wait_low_load() {
  [ "${BENCH_WAIT_LOW_LOAD:-1}" = "0" ] && return 0
  local max="${BENCH_LOAD_MAX:-4}" need="${BENCH_LOAD_STABLE:-1}" win="${BENCH_LOAD_WINDOW:-5}"
  local ok=0 busy waited=0 first=""
  while :; do
    busy=$(_bench_busy_cores "$win")
    [ -z "$first" ] && first="$busy"
    if awk -v b="$busy" -v m="$max" 'BEGIN{exit !(b < m)}'; then
      ok=$((ok + 1))
      [ "$ok" -ge "$need" ] && break
    else
      ok=0
    fi
    waited=$((waited + win))
  done
  echo "[load] $busy busy cores (first sample $first, max $max), waited ${waited}s" >&2
}
