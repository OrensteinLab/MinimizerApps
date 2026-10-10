# Run a multi-config experiment from experiments/*.yaml

reset_run_options() {
  DATASETS="$ALL_DATASETS"
  KS="28,55"
  INPUTS="gzip"
  THREADS=12
  MEM_OVERRIDE=""
  BINS=""
  STRICT_MEM=0
  RUN_TAG=""
  SIGNATURE_ORDER="all"
  ORDER_LABEL=""
  CUSTOM_ORDER=""
  N_SEEDS=1
  SIGNATURE_LEN=""
  GM_LUT=""
  FORCE=0
}

normalize_signature_order() {
  # canonical order names only (benchmark/common/order_names.py), no aliases
  SIGNATURE_ORDER="${SIGNATURE_ORDER,,}"
  case "$SIGNATURE_ORDER" in
    default|lex|xor|xor-pcg-hash|greedymini|all) ;;
    *) die "unknown signature order '$SIGNATURE_ORDER' (use default, lex, xor, xor-pcg-hash, greedymini or all)" ;;
  esac
  case "${GM_LUT:-dna}" in dna|binary) ;; *) die "unknown gm_lut '$GM_LUT' (use dna or binary)" ;; esac
  case "${GM_PROJECTION:-at0-cg1}" in at0-cg1|at1-cg0) ;; *) die "unknown projection '$GM_PROJECTION' (use at0-cg1 or at1-cg0)" ;; esac
}

run_experiment() {
  local spec="${EXPERIMENT:?experiment name or path required}"
  local -a env_lines
  mapfile -t env_lines < <(bench_py experiment shell-env-all "$spec")
  local n="${#env_lines[@]}" i

  echo "=== experiment $spec ($n configs) ==="
  KMC_BIN="$(resolve_kmc)"
  KMC_VERSION="$(kmc_version "$KMC_BIN")"

  for ((i = 0; i < n; i++)); do
    reset_run_options
    eval "${env_lines[$i]}"
    normalize_signature_order
    echo "=== [$((i + 1))/$n] datasets=$DATASETS k=$KS t=$THREADS bins=${BINS:-512} mem=${MEM_OVERRIDE:-csv} sm=$STRICT_MEM -O $SIGNATURE_ORDER $(date) ==="
    run_schemes
  done
  export_results
}

list_experiments() {
  bench_py experiment list
}
