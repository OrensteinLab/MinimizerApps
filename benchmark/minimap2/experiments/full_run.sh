#!/usr/bin/env bash
# Full minimap2 benchmark (all datasets, both repeat-filter conditions, scheme analysis, export)
# into a fresh results dir, one step at a time under low machine load.
#
#   RESULTS_NAME=results THREADS=16 experiments/full_run.sh
#
# Per dataset: `sweep` with minimap2's own default mid_occ, then `sweep` with MID_OCC fixed to
# native's effective threshold (every order filters repeats by the same rule), then
# analysis/run_seed_experiments.sh (seed-0 indexes, fixed random read sample, 5 timing replicates).
# Finally accuracy_paired.py, stats.py and `export`. Every timed step first waits for a quiet
# machine (benchmark/common/wait_low_load.sh; BENCH_LOAD_MAX busy cores, default 4).
# Re-running skips runs that already have a meta.json.
set -uo pipefail
fail=0   # a failed step is reported and the run goes on; the exit status says so (full_run_all.sh)
HERE="$(cd "$(dirname "$0")/.." && pwd)"
cd "$HERE"
source lib/env.sh
load_bench_env
export MINIMAP2_RESULTS_DIR="$MINIMAP2_BENCH_ROOT/${RESULTS_NAME:-results}"
export THREADS="${THREADS:-16}"
mkdir -p "$MINIMAP2_RESULTS_DIR/runs" "$MINIMAP2_RESULTS_DIR/csvs"
source lib/constants.sh
echo "results: $MINIMAP2_RESULTS_DIR  threads: $THREADS  $(date -Is)"

# dataset:n_reads of the scheme-analysis read sample
DATASETS="${DATASETS:-chm13_chrX:3000 grch38_wg:500 chm13_rel3_ont:1000}"
REPS="${REPS:-5}"
for entry in $DATASETS; do
  ds="${entry%%:*}"; n="${entry##*:}"
  echo "=== $ds: sweep (default mid_occ) $(date -Is)"
  DATASET="$ds" ./minimap2_bench.sh sweep || { echo "sweep failed for $ds (continuing)" >&2; fail=1; }
  cap=$(ls "$MINIMAP2_RESULTS_DIR"/runs/"$ds"/k*_w*/native/0/meta.json 2>/dev/null | grep -v '_f[0-9]' | head -1 \
        | xargs -r "$PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1]))["mid_occ"] or "")')
  if [ -n "$cap" ]; then
    echo "=== $ds: sweep (fixed mid_occ $cap) $(date -Is)"
    DATASET="$ds" MID_OCC="$cap" ./minimap2_bench.sh sweep || { echo "capped sweep failed for $ds (continuing)" >&2; fail=1; }
    echo "=== $ds: scheme analysis $(date -Is)"
    analysis/run_seed_experiments.sh "$ds" "$n" "$REPS" default "$cap" || { echo "analysis failed for $ds (continuing)" >&2; fail=1; }
  else
    echo "no native mid_occ for $ds: skipping the fixed-cap sweep and analysis" >&2
  fi
done
echo "=== paired accuracy + statistics $(date -Is)"
"$PYTHON" analysis/accuracy_paired.py --results-dir "$MINIMAP2_RESULTS_DIR" --data-dir "$MINIMAP2_DATA_DIR" || { echo "accuracy_paired failed" >&2; fail=1; }
"$PYTHON" analysis/stats.py --results-dir "$MINIMAP2_RESULTS_DIR" || { echo "stats failed" >&2; fail=1; }
./minimap2_bench.sh export || { echo "export failed" >&2; fail=1; }
echo "=== done $(date -Is)$([ "$fail" = 0 ] || echo "  (some step failed, see above)")"
exit "$fail"
