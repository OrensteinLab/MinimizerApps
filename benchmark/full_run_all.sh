#!/usr/bin/env bash
# Full benchmark of all four apps (measurement protocol: benchmark/README.md). The apps run in SLOTS,
# one slot per NUMA node, the slots in parallel and the apps of a slot one after the other; every
# timed step of an app runs on its slot's node (CPUs chosen there, memory bound to it, node-local
# gate: common/pin_cpus.sh, bench_step.py, wait_low_load.sh), so all schemes of an app are measured
# under the same node. The two slots share no cores, L3 cache, memory controller or memory; they do
# share the disk (/localdata) and the interconnect, which the step records track (disk_busy_frac,
# start/end epochs):
#
#   benchmark/full_run_all.sh                                        # results into <app root>/$RESULTS_NAME
#   SLOTS="0:minimap2,kmc 1:kraken2,sshash" RESULTS_NAME=<name> benchmark/full_run_all.sh   # the default
#   SLOTS="0:kmc" benchmark/full_run_all.sh                           # one app on node 0
#
# Default assignment: node 1, which usually has more free 2 MB blocks for the huge-page-sensitive
# single-thread steps, runs kraken2 then sshash; node 0 runs minimap2 then kmc. The two chains take
# about the same time, and kraken2's disk-heavy builds overlap only kmc's small datasets.
#
# Refuses to start into a results directory that already holds runs: a rerun must not mix runs made
# under other conditions (another protocol or another version of a scheme) with the new ones. After an
# interruption set RESUME=1 to continue the same directory (finished runs are skipped; the app logs
# are appended to). The kernel settings the timings depend on (THP, NUMA balancing, zone reclaim)
# and the commit are written to <results>/host_settings.txt. Exits 1 if any app failed.
# Per-app logs: <app root>/<RESULTS_NAME>.log; each app's experiments/full_run.sh ends with its export.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
export RESULTS_NAME="${RESULTS_NAME:-results}"
SLOTS="${SLOTS:-0:minimap2,kmc 1:kraken2,sshash}"

root_of() {  # <app>_BENCH_ROOT from the app's .env
  local var
  case "$1" in minimap2) var=MINIMAP2_BENCH_ROOT ;; sshash) var=SSHASH_BENCH_ROOT ;; kraken2) var=KRAKEN2_BENCH_ROOT ;; kmc) var=KMC3_BENCH_ROOT ;; *) return 1 ;; esac
  sed -n "s/^\(export \)\?$var=//p" "$HERE/$1/.env" | head -1 | tr -d "\"' "
}

nodes=(); seen_apps=" "
for slot in $SLOTS; do
  node="${slot%%:*}"; apps="${slot#*:}"
  [ -d "/sys/devices/system/node/node$node" ] || { echo "error: slot '$slot': no NUMA node $node" >&2; exit 1; }
  case " ${nodes[*]:-} " in *" $node "*) echo "error: node $node is in two slots" >&2; exit 1 ;; esac
  nodes+=("$node")
  for app in ${apps//,/ }; do
    case "$seen_apps" in *" $app "*) echo "error: $app is in two slots" >&2; exit 1 ;; esac
    seen_apps+="$app "
    root="$(root_of "$app")"
    [ -n "$root" ] || { echo "error: no bench root for $app (benchmark/$app/.env)" >&2; exit 1; }
    if [ -d "$root/$RESULTS_NAME/runs" ] && [ -n "$(ls -A "$root/$RESULTS_NAME/runs" 2>/dev/null)" ] && [ "${RESUME:-0}" != 1 ]; then
      echo "error: $root/$RESULTS_NAME/runs is not empty; clean it, pick another RESULTS_NAME, or set RESUME=1" >&2
      exit 1
    fi
  done
done

host_settings() {  # the kernel settings the timings depend on, once per results directory
  local f
  echo "date $(date -Is)  host $(hostname)  commit $(git -C "$HERE/.." describe --always --dirty 2>/dev/null)"
  for f in /sys/kernel/mm/transparent_hugepage/enabled /sys/kernel/mm/transparent_hugepage/defrag \
           /proc/sys/kernel/numa_balancing /proc/sys/vm/zone_reclaim_mode /proc/sys/vm/swappiness; do
    echo "$f: $(cat "$f" 2>/dev/null)"
  done
  grep -m1 "model name" /proc/cpuinfo
  for f in /sys/devices/system/node/node*/meminfo; do grep -E "MemTotal|MemFree" "$f"; done
}

run_slot() {  # <node> <app,app,...>; returns the number of failed apps
  export BENCH_NODE="$1"
  local app root failed=0
  for app in ${2//,/ }; do
    root="$(root_of "$app")"
    mkdir -p "$root/$RESULTS_NAME"
    host_settings >> "$root/$RESULTS_NAME/host_settings.txt"   # appended on RESUME
    echo "=== node $BENCH_NODE: $app -> $root/$RESULTS_NAME  $(date -Is)"
    { [ "${RESUME:-0}" = 1 ] && echo "=== RESUME $(date -Is) (node $BENCH_NODE)"; } >> "$root/$RESULTS_NAME.log"
    if ! ( cd "$HERE/$app" && bash experiments/full_run.sh ) >> "$root/$RESULTS_NAME.log" 2>&1; then
      echo "warning: $app full_run.sh exited non-zero (see $root/$RESULTS_NAME.log)" >&2
      failed=$((failed + 1))
    fi
    echo "=== node $BENCH_NODE: $app finished $(date -Is)"
  done
  return "$failed"
}

echo "=== $RESULTS_NAME, slots: $SLOTS  $(git -C "$HERE/.." describe --always --dirty 2>/dev/null)  $(date -Is)"
pids=()
for slot in $SLOTS; do
  run_slot "${slot%%:*}" "${slot#*:}" &
  pids+=($!)
done
status=0
for pid in "${pids[@]}"; do
  wait "$pid" || status=1
done
echo "=== all done $(date -Is)$([ "$status" = 0 ] || echo "  (some app failed, see the warnings above)")"
exit "$status"
