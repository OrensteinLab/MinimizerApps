#!/usr/bin/env bash
# Source this file, then call `bench_pin_cpus <n_threads> [STEP_RECORD.json]` before each timed step
# and put "${BENCH_PIN[@]}" directly in front of the timed program, inside /usr/bin/time:
#
#     bench_pin_cpus "$THREADS" "$run_dir/step_build.json"
#     "$TIME_BIN" -v -o "$t" ${BENCH_PIN[@]+"${BENCH_PIN[@]}"} prog args...
#
# Pins the step to ONE NUMA node and one logical CPU per physical core, so a timed run cannot spread
# onto cores another user is loading and never shares a core with a second hardware thread of its
# own; the step's helper threads/processes (readers, gzip, IO threads) stay inside the pinned set
# (Wood et al., Kraken 2, Genome Biol 2019 pin each classifier the same way). It takes N +
# BENCH_PIN_SPARE (default 2) cores, chosen by pin_cpus.py among the least busy (call after
# wait_low_load). The program runs under bench_step.py, which sets the CPU affinity, binds the step's
# memory to the same node (set_mempolicy MPOL_BIND; BENCH_MEMBIND=0 disables) and, given
# STEP_RECORD.json, records the node's memory state and the step's page-size, remote-memory and
# foreign-load figures (see bench_step.py).
#   BENCH_NODE=<k>          the slot's NUMA node: CPUs are chosen on node k only (full_run_all.sh sets
#                           it per slot); unset: the node whose best cores are least busy
#   BENCH_CPUS=<list>       use exactly this CPU list (e.g. 0-23 or 4,6,8) instead of choosing
#   BENCH_PIN_SPARE=<n>     spare cores for helper threads (default 2)
#   BENCH_PIN_WINDOW=<s>    load sampling window (default 1)
# Pinning is not optional: every timed step of every benchmark is pinned, so all runs are made under
# the same conditions; if the CPU selection fails the benchmark stops with an error instead of timing
# unpinned. Sets BENCH_PIN (array), BENCH_PIN_LIST (e.g. "2,6,10") and BENCH_PIN_NODE (NUMA node),
# and prints "[pin] node N: <list>" to stderr.
_BENCH_PIN_PY="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/pin_cpus.py"

_BENCH_STEP_PY="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/bench_step.py"

bench_pin_cpus() {
  local n="${1:?bench_pin_cpus: thread count}" record="${2:-}" out
  BENCH_PIN=(); BENCH_PIN_LIST=""; BENCH_PIN_NODE=""
  if [ -n "${BENCH_CPUS:-}" ]; then
    BENCH_PIN_LIST="$BENCH_CPUS"
    BENCH_PIN_NODE="${BENCH_NODE:-}"
  else
    if ! out="$(python3 "$_BENCH_PIN_PY" "$n" --spare "${BENCH_PIN_SPARE:-2}" --window "${BENCH_PIN_WINDOW:-1}" \
                ${BENCH_NODE:+--node "$BENCH_NODE"})" || [ -z "$out" ]; then
      echo "error: CPU selection failed (pin_cpus.py); the benchmark always pins its timed steps" >&2; exit 1
    fi
    BENCH_PIN_NODE="${out%% *}"; BENCH_PIN_LIST="${out#* }"
  fi
  BENCH_PIN=(python3 "$_BENCH_STEP_PY" --cpus "$BENCH_PIN_LIST")
  [ -n "$BENCH_PIN_NODE" ] && BENCH_PIN+=(--node "$BENCH_PIN_NODE")
  [ "${BENCH_MEMBIND:-1}" = "0" ] && BENCH_PIN+=(--no-membind)
  [ -n "$record" ] && BENCH_PIN+=(--record "$record")
  BENCH_PIN+=(--)
  echo "[pin] node ${BENCH_PIN_NODE:-?}: $BENCH_PIN_LIST$([ "${BENCH_MEMBIND:-1}" = "0" ] || echo ", memory bound to the node")" >&2
}
