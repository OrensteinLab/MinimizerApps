#!/usr/bin/env bash
# minimap2 scheme-comparison benchmark harness. See README.md.
set -uo pipefail
HERE_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE_DIR"

source lib/env.sh
source lib/constants.sh
source lib/schemes.sh
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
  run <scheme> [seed] [seed_idx] Run one scheme (see lib/schemes.sh ALL_SCHEMES);
                                the greedymini schemes use PROJECTION (at0-cg1
                                default or at1-cg0) and GM_TABLE (greedymini,
                                w10_k15_canonical) / MOD_GM_TABLE
                                (mod-minimizer-greedymini, w20_k9_canonical)
  sweep [schemes...]           Run every default scheme (or the given list);
                                every scheme (native too: its repeats
                                give the baseline replicates) is swept across
                                config.yaml's full seed list;
                                the greedymini schemes run once per projection
                                in PROJECTIONS (at0-cg1,at1-cg0), each its own run
  export [--datasets a,b]      Aggregate runs/**/meta.json -> results/csvs/*.csv + results/minimap2_results.xlsx
  list                         List schemes sweep would run
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
  local schemes=("$@")
  if [ ${#schemes[@]} -eq 0 ]; then
    IFS=',' read -r -a schemes <<< "$DEFAULT_SCHEMES"
  fi
  # Seed-major: every scheme (and greedymini projection) runs for seed 0, then for seed 1, ...,
  # so time-varying conditions (load, the other benchmark slot's disk use) hit all schemes alike.
  local fail=0 i s p
  for i in "${!ALL_SEEDS[@]}"; do
    for s in "${schemes[@]}"; do
      local projs=("")
      if scheme_needs_gm "$s"; then
        IFS=',' read -r -a projs <<< "$PROJECTIONS"
      fi
      for p in "${projs[@]}"; do
        [ -n "$p" ] && PROJECTION="$p"
        run_one "$s" "${ALL_SEEDS[$i]}" "$i" || fail=1
      done
    done
  done
  return $fail
}

cmd_export() {
  "$PYTHON" lib/export.py --results-dir "$MINIMAP2_RESULTS_DIR" \
    --data-dir "$MINIMAP2_DATA_DIR" --samtools "$SAMTOOLS" "$@"
}

cmd_list() {
  local s p
  for s in ${DEFAULT_SCHEMES//,/ }; do
    if scheme_needs_gm "$s"; then
      for p in ${PROJECTIONS//,/ }; do echo "$s --gm-projection $p ($(scheme_gm_table "$s"))"; done
    else
      echo "$s"
    fi
  done
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
