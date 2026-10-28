#!/usr/bin/env bash
# Streaming-query profile of a finished full run (experiments/full_run.sh): for every run of the
# seeds SEEDS (default the full run's five, 0-4), rebuild its index with the build command recorded in its meta.json, using a
# profiling build of sshash (-DSSHASH_PROFILE=On, in sshash/build-prof), into a temporary directory,
# and stream its query FASTQ once. The profiling build counts the lookups that decode a bucket and the
# positions they decode, and the lookups that verify candidate positions against the stored strings
# and the positions they verify; `sshash query` adds the counts to its JSON line, written to
#   <run>/query_profile.log
# (no timed file of the run is touched; the temporary index is deleted). Then it re-exports
# (sshash_bench.sh export), which adds the query_prof_* columns to the runs and summary sheets of the workbook. The counts
# depend only on the index and the reads, not on timing; the runs still follow the measurement
# protocol (node gate, pinning, memory binding).
#
# It works on a mirror of the full run, RESULTS_NAME (default results-prof), made from SOURCE_RESULTS
# (default results) by common/mirror_results.py: large files hard-linked, small ones copied, so the
# profiles and the re-export are written there and the source result set is never modified.
#
#   experiments/profile_run.sh
#   SEEDS="0" DATASETS="cod,kestrel" SKIP="" ...     FORCE=1 ...     NO_EXPORT=1 ...
set -uo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$HERE"
if [ -z "${SSHASH_BENCH_ROOT:-}" ] && [ -f .env ]; then
  SSHASH_BENCH_ROOT="$(sed -n 's/^SSHASH_BENCH_ROOT=//p' .env | tr -d "\"' ")"
fi
: "${SSHASH_BENCH_ROOT:?SSHASH_BENCH_ROOT is not set (benchmark/sshash/.env)}"
SEEDS="${SEEDS:-0 1 2 3 4}"
DATASETS="${DATASETS:-cod,kestrel,human,ncbi-virus,se,hprc}"
# run directories to skip (regex on the path below the dataset): the 13-mer GreedyMini tables, not used
SKIP="${SKIP:-/greedymini/k13_}"

# the harness venv (benchmark/sshash/.venv), as sshash_bench.sh
PYTHON="${PYTHON:-$HERE/.venv/bin/python}"
[ -x "$PYTHON" ] || { echo "error: no harness venv $PYTHON (make venv)" >&2; exit 1; }
RES="$SSHASH_BENCH_ROOT/${RESULTS_NAME:-results-prof}"
"$PYTHON" "$HERE/../common/mirror_results.py" "$SSHASH_BENCH_ROOT/${SOURCE_RESULTS:-results}" "$RES" || exit 1

# shellcheck source=../common/wait_low_load.sh
source "$HERE/../common/wait_low_load.sh"
# shellcheck source=../common/pin_cpus.sh
source "$HERE/../common/pin_cpus.sh"

# Profiling build, separate from the timed binary (sshash/build)
BIN="$REPO/sshash/build-prof/sshash"
mkdir -p "$REPO/sshash/build-prof"
( cd "$REPO/sshash/build-prof" && cmake .. -DCMAKE_BUILD_TYPE=Release -DSSHASH_PROFILE=On >/dev/null \
  && make -j16 sshash >&2 ) || { echo "error: profiling build failed" >&2; exit 1; }

fail=0
for ds in ${DATASETS//,/ }; do
  for seed in $SEEDS; do
    for meta in $(find "$RES/runs/$ds" -path "*/$seed/meta.json" | sort); do
      run="$(dirname "$meta")"
      [ -n "$SKIP" ] && [[ "$run" =~ $SKIP ]] && continue
      [ -s "$run/query_profile.log" ] && grep -q prof_full_lookups "$run/query_profile.log" \
        && [ "${FORCE:-0}" != 1 ] && continue
      tmp="$(mktemp -d "${TMPDIR:-/tmp}/sshash_profile.XXXXXX")"
      # the recorded commands, with the profiling binary and the temporary index
      mapfile -t cmds < <("$PYTHON" - "$meta" "$BIN" "$run" "$tmp" <<'PY'
import json, shlex, sys
meta, bin_, run, tmp = sys.argv[1:]
cmds = json.load(open(meta))["command"]
build = next(c for c in cmds if " build " in c)
query = next(c for c in cmds if " query " in c)
# the recorded commands hold the paths of the run that made them: set every output and the index to tmp
redirect = {"-o": tmp + "/index.sshash", "-d": tmp + "/tmp", "-i": tmp + "/index.sshash"}
for c, flags in ((build, ("-o", "-d")), (query, ("-i",))):
    a = shlex.split(c)
    a[0] = bin_
    for i in range(1, len(a) - 1):
        if a[i] in flags:
            a[i + 1] = redirect[a[i]]
    assert all(redirect[f] in a for f in flags), c
    print(shlex.join(a))
PY
      )
      threads="$("$PYTHON" -c "import json,sys;print(json.load(open(sys.argv[1])).get('threads',16))" "$meta")"
      mkdir -p "$tmp/tmp"
      wait_low_load 32 0 $((threads + 2))
      bench_pin_cpus "$threads"
      if ! eval ${BENCH_PIN[@]+"${BENCH_PIN[*]}"} "${cmds[0]}" > "$tmp/build.log" 2>&1; then
        echo "profile build failed: $run (log kept in $tmp)" >&2; fail=1; continue
      fi
      wait_low_load 8 6 2
      BENCH_PIN_SPARE=1 bench_pin_cpus 1
      eval ${BENCH_PIN[@]+"${BENCH_PIN[*]}"} "${cmds[1]}" > "$run/query_profile.log" 2>&1 || fail=1
      echo "$(date +%T) ${run#"$RES"/runs/}" >&2
      rm -rf "$tmp"
    done
  done
done

if [ "${NO_EXPORT:-0}" != 1 ]; then
  ./sshash_bench.sh export --results-dir "$RES" >&2 || fail=1
fi
exit "$fail"
