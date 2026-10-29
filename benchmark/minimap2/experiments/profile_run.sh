#!/usr/bin/env bash
# Mapping profile of a finished full run (experiments/full_run.sh): for every run of the
# fixed-threshold blocks (mid_occ fixed to native's) of the seeds SEEDS (default the full run's five,
# 0-4), map the first N_READS reads (default 3000) of its dataset with one thread, using a profiling
# build of minimap2 (-DMM_PHASE_PROF, in build-prof/minimap2) on the run's index and threshold. The
# full run keeps the index of seed 0 only; for the other seeds the index is rebuilt with the command
# recorded in the run's build.log ([M::main] CMD:) into a temporary directory, and deleted after. The build times the
# phases of mapping a read -- computing its selected k-mers, collecting their hits in the index,
# chaining, aligning -- and prints their shares, written to
#   <run>/map_profile.log
# (no timed file of the run is touched). Then it re-exports (minimap2_bench.sh export), which adds the
# map_prof_*_pct columns to the runs sheet. The runs follow the measurement protocol (node gate,
# pinning, memory binding).
#
# It works on a mirror of the full run, RESULTS_NAME (default results-prof), made from SOURCE_RESULTS
# (default results) by common/mirror_results.py: large files hard-linked, small ones copied, so the
# profiles and the re-export are written there and the source result set is never modified.
#
#   experiments/profile_run.sh
#   SEEDS="0" N_READS=1000 DATASETS="chm13_chrX" ...     FORCE=1 ...     NO_EXPORT=1 ...
set -uo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$HERE"
source lib/env.sh
load_bench_env
source lib/constants.sh
: "${MINIMAP2_DATA_DIR:?MINIMAP2_DATA_DIR is not set (benchmark/minimap2/.env)}"
export MINIMAP2_RESULTS_DIR="$MINIMAP2_BENCH_ROOT/${RESULTS_NAME:-results-prof}"
"$PYTHON" "$HERE/../common/mirror_results.py" "$MINIMAP2_BENCH_ROOT/${SOURCE_RESULTS:-results}" "$MINIMAP2_RESULTS_DIR" || exit 1
SEEDS="${SEEDS:-0 1 2 3 4}"
N_READS="${N_READS:-3000}"
DATASETS="${DATASETS:-chm13_chrX,grch38_wg,chm13_rel3_ont}"

# shellcheck source=../common/wait_low_load.sh
source "$HERE/../common/wait_low_load.sh"
# shellcheck source=../common/pin_cpus.sh
source "$HERE/../common/pin_cpus.sh"

# Profiling build, separate from the timed binary (minimap2/minimap2); common/ beside it for order.cpp
B="$REPO/build-prof"
mkdir -p "$B"
rsync -a --delete --exclude '*.o' --exclude '*.a' --exclude minimap2 --exclude sdust --exclude minimap2-lite \
  "$REPO/minimap2/" "$B/minimap2/"
rsync -a --delete "$REPO/common/" "$B/common/"
make -C "$B/minimap2" minimap2 CFLAGS="-g -Wall -O2 -Wc++-compat -DMM_PHASE_PROF" >&2 \
  || { echo "error: profiling build failed" >&2; exit 1; }
BIN="$B/minimap2/minimap2"

tmp="$(mktemp -d "${TMPDIR:-/tmp}/minimap2_profile.XXXXXX")"
trap 'rm -rf "$tmp"' EXIT
fail=0
for ds in ${DATASETS//,/ }; do
  reads="$tmp/$ds.fa"
  awk -v n="$N_READS" '/^>/{c++} c>n{exit} {print}' "$MINIMAP2_DATA_DIR/$ds/reads.fa" > "$reads"
  for seed in $SEEDS; do
    for meta in $(find "$MINIMAP2_RESULTS_DIR/runs/$ds" -path "*_f*/$seed/meta.json" | sort); do
      run="$(dirname "$meta")"
      [ -s "$run/map_profile.log" ] && grep -q "M::PHASE" "$run/map_profile.log" && [ "${FORCE:-0}" != 1 ] && continue
      read -r preset cap threads < <("$PYTHON" -c "import json,sys;m=json.load(open(sys.argv[1]));print(m['preset'],m['mid_occ_cap'],m.get('threads') or 16)" "$meta")
      mmi="$run/index.mmi"
      if [ ! -f "$mmi" ]; then   # seeds other than 0: rebuild the index with the recorded command
        cmd="$(sed -n 's/^\[M::main\] CMD: //p' "$run/build.log")"
        [ -n "$cmd" ] || { echo "skip (no index, no CMD in build.log): $run" >&2; fail=1; continue; }
        read -r -a args <<< "$cmd"
        for i in "${!args[@]}"; do [ "${args[$i]}" = -d ] && args[$((i + 1))]="$tmp/index.mmi"; done
        args[0]="$BIN"; mmi="$tmp/index.mmi"
        wait_low_load 16 0 $((threads + 2))
        bench_pin_cpus "$threads"
        ${BENCH_PIN[@]+"${BENCH_PIN[@]}"} "${args[@]}" > /dev/null 2> "$tmp/build.log" \
          || { echo "index rebuild failed: $run" >&2; fail=1; continue; }
      fi
      wait_low_load 16 0 2
      BENCH_PIN_SPARE=1 bench_pin_cpus 1
      ${BENCH_PIN[@]+"${BENCH_PIN[@]}"} "$BIN" -x "$preset" -c -f "$cap" -t 1 "$mmi" "$reads" \
        > /dev/null 2> "$run/map_profile.log" || fail=1
      rm -f "$tmp/index.mmi"
      echo "$(date +%T) ${run#"$MINIMAP2_RESULTS_DIR"/runs/}" >&2
    done
  done
done

if [ "${NO_EXPORT:-0}" != 1 ]; then
  ./minimap2_bench.sh export >&2 || fail=1
fi
exit "$fail"
