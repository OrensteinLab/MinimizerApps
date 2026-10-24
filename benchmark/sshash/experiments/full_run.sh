#!/usr/bin/env bash
# Full sshash benchmark (all datasets x all orders/variants x 5 seeds, threads 16, canonical
# modality) into a fresh results dir, each timed step waiting for a quiet machine
# (benchmark/common/wait_low_load.sh), then export.
#
#   RESULTS_NAME=results experiments/full_run.sh
#
# Re-running skips runs that already have a meta.json. The sshash binary must be built with the
# default nucleotide encoding (GreedyMini needs it): SSHASH_BIN, default <repo>/sshash/build/sshash.
set -uo pipefail
fail=0   # a failed step is reported and the run goes on; the exit status says so (full_run_all.sh)
HERE="$(cd "$(dirname "$0")/.." && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$HERE"
if [ -z "${SSHASH_BENCH_ROOT:-}" ] && [ -f .env ]; then
  SSHASH_BENCH_ROOT="$(sed -n 's/^SSHASH_BENCH_ROOT=//p' .env | tr -d "\"' ")"
fi
: "${SSHASH_BENCH_ROOT:?SSHASH_BENCH_ROOT is not set (benchmark/sshash/.env)}"
export SSHASH_BENCH_ROOT
RES="$SSHASH_BENCH_ROOT/${RESULTS_NAME:-results}"
BIN="${SSHASH_BIN:-$REPO/sshash/build/sshash}"
[ -x "$BIN" ] || { echo "error: sshash binary $BIN not found" >&2; exit 1; }
mkdir -p "$RES"
echo "results: $RES  binary: $BIN  $(date -Is)"
./sshash_bench.sh run -d "${DATASETS:-cod,kestrel,human,ncbi-virus,se,hprc}" -O all -n 5 -t 16 \
  --sshash "$BIN" --results-dir "$RES" || { echo "run failed (re-run resumes)" >&2; fail=1; }
./sshash_bench.sh export --results-dir "$RES" || { echo "export failed" >&2; fail=1; }
echo "=== done $(date -Is)$([ "$fail" = 0 ] || echo "  (some step failed, see above)")"
exit "$fail"
