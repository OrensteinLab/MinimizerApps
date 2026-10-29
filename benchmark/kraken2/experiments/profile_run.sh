#!/usr/bin/env bash
# Classification profile of a finished full run (experiments/full_run.sh): for every database
# (order variant x seed) and query set, one classify run with a classify built with
# ENABLE_MINIMIZER_PROF=1, which times the three phases of each read pair (computing the
# selected k-mers into a token stream, the hash-table lookups, resolving the hits) and counts the
# lookups, the lookups that find no taxon, and the GreedyMini ranks taken on the finite (UHS) and
# the infinite path. It writes
#   <run>/results/classify_profile_<dataset>_mhg<m>.txt
# next to the timed results (no timed file is touched), then re-exports (run_benchmark.sh -c),
# which turns those files into
# the `profile` sheet of kraken2_results.xlsx (mean / std over seeds, ratio to xor). These runs are
# not the timed ones: the profiling build is separate (kraken2/src-prof) and the timed classify
# binary is left as it is. Every run follows the measurement protocol (node gate, pinning, memory
# binding: benchmark/common/wait_low_load.sh, pin_cpus.sh).
#
# It works on a mirror of the full run, RESULTS_NAME (default results-prof), made from SOURCE_RESULTS
# (default results) by common/mirror_results.py: large files hard-linked, small ones copied, so the
# profiles and the re-export are written there and the source result set is never modified.
#
#   experiments/profile_run.sh                                  # all seeds, minimum hit groups 2
#   SOURCE_RESULTS=results RESULTS_NAME=results-prof ...
#   FORCE=1 ...                                                 # redo runs that have a profile file
#   NO_EXPORT=1 ...                                             # skip the re-export
set -uo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
cd "$HERE"
if [ -z "${KRAKEN2_BENCH_ROOT:-}" ] && [ -f .env ]; then
  KRAKEN2_BENCH_ROOT="$(sed -n 's/^export KRAKEN2_BENCH_ROOT=//p;s/^KRAKEN2_BENCH_ROOT=//p' .env | head -1 | tr -d "\"' ")"
fi
: "${KRAKEN2_BENCH_ROOT:?KRAKEN2_BENCH_ROOT is not set (benchmark/kraken2/.env)}"
export KRAKEN2_BENCH_ROOT
# the harness's conda env (CONDA_ENV in .env, default kraken_env; setup.sh)
if [ -z "${PYTHON:-}" ]; then
  CONDA_ENV="${CONDA_ENV:-$(sed -n 's/^export CONDA_ENV="\${CONDA_ENV:-\(.*\)}"$/\1/p;s/^CONDA_ENV=//p' .env 2>/dev/null | head -1 | tr -d "\"' ")}"
  CONDA_ENV="${CONDA_ENV:-kraken_env}"
  env_dir="$(command -v conda >/dev/null 2>&1 && conda env list 2>/dev/null | awk -v e="$CONDA_ENV" '$1==e{print $NF}')"
  [ -z "$env_dir" ] && env_dir="$HOME/.conda/envs/$CONDA_ENV"
  PYTHON="$env_dir/bin/python"
fi
[ -x "$PYTHON" ] || { echo "error: no python in the conda env ($PYTHON)" >&2; exit 1; }
SRC_RESULTS="$KRAKEN2_BENCH_ROOT/${SOURCE_RESULTS:-results}"
export KRAKEN2_RESULTS_DIR="$KRAKEN2_BENCH_ROOT/${RESULTS_NAME:-results-prof}"
"$PYTHON" "$HERE/../common/mirror_results.py" "$SRC_RESULTS" "$KRAKEN2_RESULTS_DIR" || exit 1
# the run manifest names each run's results directory by absolute path: point it at the mirror (a copied file)
for m in "$KRAKEN2_RESULTS_DIR"/runs/*/.benchmark_run_manifest; do
  [ -f "$m" ] && sed -i "s#$SRC_RESULTS/#$KRAKEN2_RESULTS_DIR/#g" "$m"
done
RUNS="$KRAKEN2_RESULTS_DIR/runs/${LIBRARY:-strain_exclusion}"
DATA="$KRAKEN2_BENCH_ROOT/data"
MHG="${MHG:-2}"
DATASETS=(${DATASETS:-strain-exclusion fda-argos cami})
[ -d "$RUNS" ] || { echo "error: no run directory $RUNS" >&2; exit 1; }

# shellcheck source=../common/wait_low_load.sh
source "$HERE/../common/wait_low_load.sh"
# shellcheck source=../common/pin_cpus.sh
source "$HERE/../common/pin_cpus.sh"

# Profiling build in a sibling of kraken2/src (its relative includes still resolve)
SRC="$HERE/../../kraken2/src"; PROF="$HERE/../../kraken2/src-prof"
rsync -a --delete --exclude '*.o' --exclude '.minimizer_prof_mode' --exclude classify \
  --exclude build_db --exclude estimate_capacity --exclude dump_table --exclude lookup_accession_numbers \
  --exclude k2mask --exclude blast_to_fasta --exclude merge --exclude libtax.so \
  --exclude benchmark_minimizers --exclude 'test_*' "$SRC/" "$PROF/"
make -C "$PROF" classify ENABLE_MINIMIZER_PROF=1 >&2 || { echo "error: profiling build failed" >&2; exit 1; }
BIN="$PROF/classify"

fail=0
# seed-major, all order variants of a seed one after the other (as the timed classify)
mapfile -t dbs < <(find "$RUNS" -path '*/db/hash.k2d' | sed 's#/db/hash.k2d##' | sort)
seeds=$(printf '%s\n' "${dbs[@]}" | sed 's#.*/##' | sort -un)
for seed in $seeds; do
  for ds in "${DATASETS[@]}"; do
    for run in "${dbs[@]}"; do
      [ "${run##*/}" = "$seed" ] || continue
      res="$run/results"; mkdir -p "$res"
      prof="$res/classify_profile_${ds}_mhg${MHG}.txt"
      [ -s "$prof" ] && [ "${FORCE:-0}" != 1 ] && continue
      wait_low_load 16 12 4
      bench_pin_cpus 2 "$res/profile_${ds}_mhg${MHG}_step.json"
      if ! env KRAKEN2_MINIMIZER_PROF="$prof" ${BENCH_PIN[@]+"${BENCH_PIN[@]}"} "$BIN" \
          -H "$run/db/hash.k2d" -t "$run/db/taxo.k2d" -o "$run/db/opts.k2d" -p 1 -P -g "$MHG" \
          -O /dev/null "$DATA/$ds/R1.fastq" "$DATA/$ds/R2.fastq" 2> "$res/profile_${ds}_mhg${MHG}.log"; then
        echo "profile failed: $run $ds (see $res/profile_${ds}_mhg${MHG}.log)" >&2; fail=1
      fi
      echo "$(date +%T) ${run#"$RUNS"/} $ds" >&2
    done
  done
done

if [ "${NO_EXPORT:-0}" != 1 ]; then
  PATH="$(dirname "$PYTHON"):$PATH" ./run_benchmark.sh -E -c >&2 || fail=1
fi
exit "$fail"
