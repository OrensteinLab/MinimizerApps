#!/usr/bin/env bash
# Source this file, then call `bench_pin_cpus <n_threads>` before each timed step and put
# "${BENCH_PIN[@]}" directly in front of the timed program, inside /usr/bin/time:
#
#     bench_pin_cpus "$THREADS"
#     "$TIME_BIN" -v -o "$t" ${BENCH_PIN[@]+"${BENCH_PIN[@]}"} prog args...
#
# Pins the step with `taskset -c` (no privileges needed) to ONE NUMA node and one logical CPU per
# physical core, so a timed run cannot spread onto cores another user is loading, never shares a
# core with a second hardware thread of its own, and keeps its memory node-local; the step's helper
# threads/processes (readers, gzip, IO threads) stay inside the pinned set (Wood et al., Kraken 2,
# Genome Biol 2019 pin each classifier the same way). Selection logic: pin_cpus.py (stdlib python3);
# call after wait_low_load. It takes N + BENCH_PIN_SPARE (default 2) cores from the quietest node.
#   BENCH_CPUS=<list>       use exactly this CPU list (e.g. 0-23 or 4,6,8) instead of choosing
#   BENCH_PIN_SPARE=<n>     spare cores for helper threads (default 2)
#   BENCH_PIN_WINDOW=<s>    load sampling window (default 1)
# Pinning is not optional: every timed step of every benchmark is pinned, so all runs are made under
# the same conditions; if taskset or the CPU selection is unavailable the benchmark stops with an error
# instead of timing unpinned. Sets BENCH_PIN (array), BENCH_PIN_LIST (e.g. "2,6,10") and BENCH_PIN_NODE
# (NUMA node; empty for BENCH_CPUS), and prints "[pin] node N: taskset -c <list>" to stderr.
_BENCH_PIN_PY="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/pin_cpus.py"

bench_pin_cpus() {
  local n="${1:?bench_pin_cpus: thread count}" out
  BENCH_PIN=(); BENCH_PIN_LIST=""; BENCH_PIN_NODE=""
  if ! command -v taskset >/dev/null 2>&1; then echo "error: taskset not found; the benchmark always pins its timed steps" >&2; exit 1; fi
  if [ -n "${BENCH_CPUS:-}" ]; then
    BENCH_PIN_LIST="$BENCH_CPUS"
  else
    if ! out="$(python3 "$_BENCH_PIN_PY" "$n" --spare "${BENCH_PIN_SPARE:-2}" --window "${BENCH_PIN_WINDOW:-1}")" || [ -z "$out" ]; then
      echo "error: CPU selection failed (pin_cpus.py); the benchmark always pins its timed steps" >&2; exit 1
    fi
    BENCH_PIN_NODE="${out%% *}"; BENCH_PIN_LIST="${out#* }"
  fi
  BENCH_PIN=(taskset -c "$BENCH_PIN_LIST")
  echo "[pin] node ${BENCH_PIN_NODE:-?}: taskset -c $BENCH_PIN_LIST" >&2
}
