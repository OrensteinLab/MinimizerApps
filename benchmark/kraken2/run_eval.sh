#!/usr/bin/env bash
# run_eval.sh — alias for paper-mode run_benchmark.sh (labeled + ARGOS).
#
# Kept so older docs/commands still work. Prefer:
#   ./benchmark/kraken2/run_benchmark.sh -o "$KRAKEN2_BENCHMARK_DIR" -S xor greedymini
#
# Former flags: -A is the default (ARGOS on); -a skips ARGOS; -T is ignored
# (seed #1 is the paper toggle). -j is parallel builds; -t is k2 build threads.

set -euo pipefail
EVAL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "${EVAL_DIR}/run_benchmark.sh" "$@"
