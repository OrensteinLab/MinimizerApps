# Run a multi-config experiment from experiments/*.yaml

reset_run_options() {
  DATASETS="$ALL_DATASETS"
  K=31
  M=11
  ORDERS="all"
  GM_TABLES=""
  GM_PROJECTIONS=""
  THREADS=8
  N_SEEDS=5
  FORCE=0
  RUN_TAG=""
}

run_experiment() {
  local spec="${EXPERIMENT:?experiment name or path required}"
  local -a env_lines
  mapfile -t env_lines < <(bench_py experiment shell-env-all "$spec")
  local n="${#env_lines[@]}" i

  echo "=== experiment $spec ($n configs) ==="
  SSHASH_BIN="$(resolve_sshash)"
  SSHASH_VERSION="$(sshash_version "$SSHASH_BIN")"

  for ((i = 0; i < n; i++)); do
    reset_run_options
    eval "${env_lines[$i]}"
    echo "=== [$((i + 1))/$n] datasets=$DATASETS k=$K m=$M orders=$ORDERS n=$N_SEEDS $(date) ==="
    run_sweep
  done
  export_results
}

list_experiments() {
  bench_py experiment list
}
