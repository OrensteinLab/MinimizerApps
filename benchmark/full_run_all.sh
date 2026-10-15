#!/usr/bin/env bash
# Full benchmark of all four apps, one after the other (never two at once: the timings need a quiet
# machine), every timed step pinned (common/pin_cpus.sh) and gated on low load (common/wait_low_load.sh):
#
#   benchmark/full_run_all.sh                       # results into <app root>/results_r2
#   RESULTS_NAME=results_r3 APPS="kmc sshash" benchmark/full_run_all.sh
#
# Refuses to start into a results directory that already holds runs: a rerun must not mix runs made
# under other conditions (e.g. unpinned) with the new ones. After an interruption set RESUME=1 to
# continue the same directory (finished runs are skipped, the rest are made pinned like the others).
# Per-app logs: <app root>/<RESULTS_NAME>.log; each app's experiments/full_run.sh ends with its export.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
export RESULTS_NAME="${RESULTS_NAME:-results_r2}"
APPS="${APPS:-minimap2 sshash kraken2 kmc}"

root_of() {  # <app>_BENCH_ROOT from the app's .env
  local var
  case "$1" in minimap2) var=MINIMAP2_BENCH_ROOT ;; sshash) var=SSHASH_BENCH_ROOT ;; kraken2) var=KRAKEN2_BENCH_ROOT ;; kmc) var=KMC3_BENCH_ROOT ;; *) return 1 ;; esac
  sed -n "s/^\(export \)\?$var=//p" "$HERE/$1/.env" | head -1 | tr -d "\"' "
}

for app in $APPS; do
  root="$(root_of "$app")"
  [ -n "$root" ] || { echo "error: no bench root for $app (benchmark/$app/.env)" >&2; exit 1; }
  if [ -d "$root/$RESULTS_NAME/runs" ] && [ -n "$(ls -A "$root/$RESULTS_NAME/runs" 2>/dev/null)" ] && [ "${RESUME:-0}" != 1 ]; then
    echo "error: $root/$RESULTS_NAME/runs is not empty; clean it, pick another RESULTS_NAME, or set RESUME=1" >&2
    exit 1
  fi
done

for app in $APPS; do
  root="$(root_of "$app")"
  echo "=== $app -> $root/$RESULTS_NAME  $(date -Is)"
  ( cd "$HERE/$app" && bash experiments/full_run.sh ) > "$root/$RESULTS_NAME.log" 2>&1 \
    || echo "warning: $app full_run.sh exited non-zero (see $root/$RESULTS_NAME.log)" >&2
  echo "=== $app finished $(date -Is)"
done
echo "=== all done $(date -Is)"
