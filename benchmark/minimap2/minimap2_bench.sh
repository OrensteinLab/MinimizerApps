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
load_bench_config

usage() {
  cat <<'EOF'
Usage: minimap2_bench.sh <command> [args]

Commands:
  make-dataset [--seed N]      (Re)generate the default synthetic dataset
  pbsim-reads [depth]          Replace reads.fa/truth.tsv with PBSIM 1.0.3
                                CLR reads as in the Winnowmap paper (10%
                                error, 15 kb mean, 100-25,000 bp; needs
                                $PBSIM1, built by setup.sh). depth defaults
                                to 5 (D1; D2 uses 0.5).
  run <order> [seed] [seed_idx] Run one order (see lib/orders.sh ALL_ORDERS)
  sweep [orders...]            Run every default order (or the given list);
                                seed-dependent orders (not native/greedymini)
                                are swept across config.yaml's full seed list
  export [--datasets a,b]      Aggregate runs/**/meta.json -> CSVs + minimap2_results.xlsx
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
  local depth="${1:-5}"
  if [ -z "${PBSIM1:-}" ] || [ ! -x "$PBSIM1" ]; then
    echo "error: PBSIM 1.0.3 not found at '$PBSIM1' -- run ./setup.sh first" >&2
    return 1
  fi
  "$PYTHON" lib/pbsim_reads.py --simulator pbsim1 --dataset-dir "$MINIMAP2_DATA_DIR/$DATASET" \
    --pbsim-bin "$PBSIM1" --model "$(dirname "$PBSIM1")/../data/model_qc_clr" \
    --depth "$depth" --length-mean 15000 --accuracy-mean 0.90 --seed 1
}

cmd_sweep() {
  local orders=("$@")
  if [ ${#orders[@]} -eq 0 ]; then
    IFS=',' read -r -a orders <<< "$DEFAULT_ORDERS"
  fi
  local fail=0
  for o in "${orders[@]}"; do
    if order_is_seed_dependent "$o" && [ "${#ALL_SEEDS[@]}" -gt 0 ]; then
      local i
      for i in "${!ALL_SEEDS[@]}"; do
        run_one "$o" "${ALL_SEEDS[$i]}" "$i" || fail=1
      done
    else
      run_one "$o" "${ALL_SEEDS[0]:-0x2545F4914F6CDD1D}" 0 || fail=1
    fi
  done
  return $fail
}

cmd_export() {
  "$PYTHON" lib/export.py --results-dir "$MINIMAP2_RESULTS_DIR" \
    --data-dir "$MINIMAP2_DATA_DIR" --samtools "$SAMTOOLS" "$@"
}

cmd_list() {
  echo "$DEFAULT_ORDERS" | tr ',' '\n'
}

case "${1:-}" in
  make-dataset) shift; cmd_make_dataset "$@" ;;
  pbsim-reads) shift; cmd_pbsim_reads "$@" ;;
  run) shift; run_one "$@" ;;
  sweep) shift; cmd_sweep "$@" ;;
  export) shift; cmd_export "$@" ;;
  list) cmd_list ;;
  *) usage; exit 1 ;;
esac
