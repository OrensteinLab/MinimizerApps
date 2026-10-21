#!/usr/bin/env bash
# Full kraken2 benchmark (strain-exclusion library -E; query sets strain-exclusion, fda-argos, cami;
# 7 order variants x 5 seeds; --minimum-hit-groups 0 and 2; build 16 threads, classify 1 thread)
# into a fresh results dir, each timed step waiting for a quiet machine
# (benchmark/common/wait_low_load.sh), then density microbench + comparison + export (one
# run_benchmark.sh invocation does all of these).
#
#   RESULTS_NAME=results experiments/full_run.sh
#
# Re-running resumes: databases with a hash.k2d are not rebuilt; classification of a seed is
# redone for that seed. Re-uses db/strain_exclusion (library + taxonomy) and data/ as they are.
# Per-run databases (~9.4 GB each, 35 of them) are all kept: export reads hash.k2d (header ->
# distinct_minimizers, size -> index_size_gb) and opts.k2d, so they cannot be deleted.
set -uo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
cd "$HERE"
if [ -z "${KRAKEN2_BENCH_ROOT:-}" ] && [ -f .env ]; then
  KRAKEN2_BENCH_ROOT="$(sed -n 's/^export KRAKEN2_BENCH_ROOT=//p;s/^KRAKEN2_BENCH_ROOT=//p' .env | head -1 | tr -d "\"' ")"
fi
: "${KRAKEN2_BENCH_ROOT:?KRAKEN2_BENCH_ROOT is not set (benchmark/kraken2/.env)}"
export KRAKEN2_BENCH_ROOT
export KRAKEN2_RESULTS_DIR="$KRAKEN2_BENCH_ROOT/${RESULTS_NAME:-results}"
mkdir -p "$KRAKEN2_RESULTS_DIR"
echo "results: $KRAKEN2_RESULTS_DIR  $(date -Is)"
CMD=(./run_benchmark.sh -E -n 5 -t 16 -j 1 -P at0-cg1,at1-cg0
     -S xor xor-pcg-hash greedymini mod-minimizer-xor mod-minimizer-greedymini)
echo "KRAKEN2_RESULTS_DIR=$KRAKEN2_RESULTS_DIR ${CMD[*]}"
"${CMD[@]}" || echo "run failed (re-run resumes)" >&2
echo "=== done $(date -Is)"
