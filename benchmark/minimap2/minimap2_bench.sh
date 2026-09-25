#!/usr/bin/env bash
# minimap2 order-comparison benchmark harness. See README.md.
set -uo pipefail
HERE_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE_DIR"

source lib/env.sh
source lib/constants.sh
source lib/orders.sh
source lib/run.sh
load_bench_env

usage() {
  cat <<'EOF'
Usage: minimap2_bench.sh <command> [args]

Commands:
  make-dataset [--seed N]      (Re)generate the default synthetic dataset
  pbsim-reads [model] [depth]  Replace reads.fa/truth.tsv with real pbsim2
                                simulated reads (needs $PBSIM; see
                                environment.yml/setup.sh). model defaults to
                                R94.model (ONT); depth defaults to 8.
  run <order> [seed]           Run one order (see lib/orders.sh ALL_ORDERS)
  sweep [orders...]            Run every default order (or the given list)
  export                       Aggregate runs/**/meta.json -> minimap2_runs.csv
  list                         List orders sweep would run
EOF
}

cmd_make_dataset() {
  mkdir -p "$MINIMAP2_DATA_DIR/$DATASET"
  "$PYTHON" lib/make_dataset.py --out-dir "$MINIMAP2_DATA_DIR/$DATASET" \
    --seed "${1:-1}" --genome-len "${GENOME_LEN:-300000}" --repeat-len "${REPEAT_LEN:-20000}" \
    --n-reads "${N_READS:-400}"
}

cmd_pbsim_reads() {
  local model="${1:-R94.model}"
  local depth="${2:-8}"
  if [ -z "${PBSIM:-}" ] || [ ! -x "$PBSIM" ]; then
    echo "error: pbsim not found -- run ./setup.sh (conda env 'minimap2_bench') first" >&2
    return 1
  fi
  local model_path="$model"
  [ -f "$model_path" ] || model_path="$(dirname "$PBSIM")/../data/$model"
  "$PYTHON" lib/pbsim_reads.py --dataset-dir "$MINIMAP2_DATA_DIR/$DATASET" \
    --pbsim-bin "$PBSIM" --model "$model_path" --depth "$depth"
}

cmd_sweep() {
  local orders=("$@")
  if [ ${#orders[@]} -eq 0 ]; then
    IFS=',' read -r -a orders <<< "$DEFAULT_ORDERS"
  fi
  local seed=${SWEEP_SEED:-0x2545F4914F6CDD1D}
  local fail=0
  for o in "${orders[@]}"; do
    run_one "$o" "$seed" || fail=1
  done
  return $fail
}

cmd_export() {
  "$PYTHON" lib/export.py --results-dir "$MINIMAP2_RESULTS_DIR" \
    --out "$MINIMAP2_RESULTS_DIR/minimap2_runs.csv"
}

cmd_list() {
  echo "$DEFAULT_ORDERS" | tr ',' '\n'
}

case "${1:-}" in
  make-dataset) shift; cmd_make_dataset "$@" ;;
  pbsim-reads) shift; cmd_pbsim_reads "$@" ;;
  run) shift; run_one "$@" ;;
  sweep) shift; cmd_sweep "$@" ;;
  export) cmd_export ;;
  list) cmd_list ;;
  *) usage; exit 1 ;;
esac
