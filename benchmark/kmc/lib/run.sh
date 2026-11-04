# kmc binary resolution, benchmark runs, and export dispatch

# wait_low_load: block until the (shared) machine is quiet before each timed kmc invocation
# (BENCH_LOAD_MAX; BENCH_WAIT_LOW_LOAD=0 disables; benchmark/common/wait_low_load.sh).
source "${BASH_SOURCE[0]%/*}/../../common/wait_low_load.sh"
# bench_pin_cpus: pin the timed kmc invocation to $THREADS least-busy CPUs of one node, memory bound to it
# (benchmark/common/pin_cpus.sh; pinning is not optional); step_count.json records the node's memory state
source "${BASH_SOURCE[0]%/*}/../../common/pin_cpus.sh"

# Elapsed seconds between two $EPOCHREALTIME samples ("sec.usec" strings),
# formatted like the `%.6f` this replaced -- pure bash, no subprocess.
wall_time_delta() {
  local start="$1" end="$2"
  local s0="${start%%.*}" u0="${start#*.}" s1="${end%%.*}" u1="${end#*.}"
  u0=$((10#$u0)); u1=$((10#$u1))
  local sec=$((s1 - s0)) usec=$((u1 - u0))
  if ((usec < 0)); then
    usec=$((usec + 1000000))
    sec=$((sec - 1))
  fi
  printf '%d.%06d' "$sec" "$usec"
}

resolve_kmc() {
  if [[ -n "$KMC" ]]; then
    [[ -f "$KMC" ]] || die "kmc binary not found: $KMC"
    echo "$(cd "$(dirname "$KMC")" && pwd)/$(basename "$KMC")"
    return
  fi
  if [[ -f "$REPO_ROOT/bin/kmc" ]]; then
    echo "$REPO_ROOT/bin/kmc"
    return
  fi
  local w
  w="$(command -v kmc 2>/dev/null || true)"
  [[ -n "$w" ]] || die "kmc not found; pass --kmc PATH (after make)"
  echo "$w"
}

kmc_version() {
  local kmc="$1"
  "$kmc" --version 2>&1 | head -1 || echo "unknown"
}

list_path_for() {
  local root="$1" ds="$2" kind="$3"
  if [[ "$kind" == "gzip" ]]; then
    echo "$root/$ds/files_gz.lst"
  else
    echo "$root/$ds/files_plain.lst"
  fi
}

params_tag() {
  local k="$1"
  echo "k${k}_p$(signature_len_effective)"
}

config_dir_name() {
  local kind="$1" mem="$2"
  if [[ -n "$RUN_TAG" ]]; then
    echo "$RUN_TAG"
    return
  fi
  local name="${kind}_n${BINS:-512}_t${THREADS}_m${mem}"
  [[ "$STRICT_MEM" -eq 1 ]] && name="${name}_sm"
  echo "$name"
}

# Run-dir segment for a scheme: scheme_names.scheme_dir_name (e.g. "xor",
# "greedymini/k9_at1-cg0"), plus the kmc-specific table LUT backend for
# greedymini ("/lut_dna" | "/lut_binary": same order, different storage).
scheme_dir_name() {
  case "$1" in
    default|lex|xor|xor-pcg-hash) bench_py gm scheme-dir --scheme "$1" ;;
    greedymini)
      echo "$(bench_py gm scheme-dir --scheme greedymini --gm-k "$(signature_len_effective)" \
        --projection "${GM_PROJECTION:-at0-cg1}")/lut_${GM_LUT:-binary}"
      ;;
    *) die "unknown scheme: $1" ;;
  esac
}

run_rel_path() {
  local ds="$1" k="$2" kind="$3" scheme="$4" seed_idx="$5" mem="$6"
  echo "${ds}/$(params_tag "$k")/$(config_dir_name "$kind" "$mem")/$(scheme_dir_name "$scheme")/${seed_idx}"
}

run_one() {
  local ds="$1" k="$2" kind="$3" scheme="$4" custom_order="$5" order_label="$6" seed="${7:-}" seed_idx="${8:-0}"

  local mem label lst kmc_cmd wrapped p cfg_tag rid
  mem="$(dataset_mem "$ds")"
  label="$(csv_field "$ds" label)"
  p="$(signature_len_effective)"
  cfg_tag="$(config_dir_name "$kind" "$mem")"
  rid="$(run_rel_path "$ds" "$k" "$kind" "$scheme" "$seed_idx" "$mem")"

  local run_dir="$RESULTS_DIR/runs/$rid"
  local tmp="$run_dir/tmp"
  local db_prefix="$run_dir/out"
  local json_path="$run_dir/kmc.json"
  local time_path="$run_dir/time.txt"
  local meta_path="$run_dir/meta.json"

  if [[ -f "$json_path" && "$FORCE" -eq 0 ]]; then
    echo "skip existing run $rid"
    return 0
  fi

  mkdir -p "$run_dir"
  rm -rf "$tmp"
  mkdir -p "$tmp"

  lst="$(list_path_for "$DATA_DIR" "$ds" "$kind")"
  [[ -f "$lst" ]] || die "missing $lst; run: $0 download -d $ds"

  kmc_cmd=(
    "$KMC_BIN"
    -ci2 -v -hp
    "-t$THREADS"
    "-m$mem"
    "-k$k"
    "-j$json_path"
  )
  if [[ -n "$BINS" ]]; then
    kmc_cmd+=("-n$BINS")
  fi
  if [[ "$STRICT_MEM" -eq 1 ]]; then
    kmc_cmd+=(-sm)
  fi
  if [[ -n "$SIGNATURE_LEN" ]]; then
    kmc_cmd+=("-p$SIGNATURE_LEN")
  elif [[ "$scheme" == "greedymini" ]]; then
    kmc_cmd+=("-p$(signature_len_effective)")
  fi
  local gm_lut="${GM_LUT:-binary}" projection="${GM_PROJECTION:-at0-cg1}"
  if [[ "$scheme" == "greedymini" ]]; then
    local gm_file="$custom_order" embedded="$custom_order" embedded
    [[ -n "$gm_file" ]] || gm_file="$(resolve_gm_path "$k" "${seed:-0}")"
    [[ -f "$gm_file" ]] || die "GM table not found for k=$k p=$(signature_len_effective) seed=${seed:-0} gm_lut=$gm_lut: $gm_file"
    if [[ "$projection" == "at1-cg0" ]]; then
      # -oat1 permutes a binary-alphabet table: a prebuilt _dna.gm cannot be flipped
      [[ "$gm_lut" == "binary" ]] || die "projection at1-cg0 needs gm_lut=binary (kmc -oat1 permutes the binary .gm)"
      [[ -n "$seed" ]] || die "GM binary .gm requires a seed"
      kmc_cmd+=("-co$gm_file" "-ol$order_label" "-os$((seed))" "-oat1")
    elif [[ "$gm_lut" == "dna" ]] && embedded="$(bench_py gm dna-seed "$gm_file" 2>/dev/null)"; then
      seed="$embedded"
      kmc_cmd+=("-co$gm_file" "-ol$order_label")
    else
      [[ -n "$seed" ]] || die "GM binary .gm requires a seed; use a _dna.gm table or -n with predefined seeds"
      kmc_cmd+=("-co$gm_file" "-ol$order_label" "-os$((seed))")
      # kmc uses the binary rank LUT for a binary .gm (-ode would DNA-expand it).
    fi
    custom_order="$gm_file"
  elif [[ "$scheme" == "lex" ]]; then
    # True canonical lex (no is_allowed denylist): -or xor -os 0 is order-
    # equivalent to MIN(i, revcomp(i)) with no exclusions (verified offline;
    # see README). Always seed 0 -- deterministic, not a real random scheme.
    kmc_cmd+=("-orxor" "-ol$order_label" "-os0")
  elif [[ "$scheme" == "xor" || "$scheme" == "xor-pcg-hash" ]]; then
    kmc_cmd+=("-or$scheme" "-ol$order_label" "-os$((seed))")
  fi
  kmc_cmd+=("@$lst" "$db_prefix" "$tmp")

  local t0 wall rc cmd_file
  cmd_file="$(mktemp)"
  wait_low_load 40 0 $((THREADS + 2))
  bench_pin_cpus "$THREADS" "$run_dir/step_count.json"
  printf '%s\n' "$BENCH_PIN_LIST" >"$run_dir/cpus.txt"
  if [[ -x "$TIME_BIN" ]]; then
    wrapped=("$TIME_BIN" -v -o "$time_path" ${BENCH_PIN[@]+"${BENCH_PIN[@]}"} "${kmc_cmd[@]}")
  else
    echo "warning: $TIME_BIN missing; wall time from KMC JSON only" >&2
    wrapped=(${BENCH_PIN[@]+"${BENCH_PIN[@]}"} "${kmc_cmd[@]}")
  fi

  printf '%s\n' "${wrapped[*]}"
  printf '%s\n' "${wrapped[@]}" >"$cmd_file"
  t0="$EPOCHREALTIME"
  set +e
  "${wrapped[@]}"
  rc=$?
  set -e
  wall="$(wall_time_delta "$t0" "$EPOCHREALTIME")"

  # kmc.json's self-reported mode (xor / xor-pcg-hash / greedymini); lex runs
  # report xor. gm_lut / projection are the harness settings that built the
  # kmc command line (kmc.json also records signature_order_gm_lut).
  local sig_mode=""
  if [[ -f "$json_path" ]]; then
    sig_mode="$("$PYTHON" -c "import json; print(json.load(open('$json_path')).get('signature_order_mode','') or '')")"
  fi
  if [[ -z "$sig_mode" ]]; then
    case "$scheme" in
      greedymini|xor|xor-pcg-hash) sig_mode="$scheme" ;;
      lex) sig_mode="xor" ;;
    esac
  fi
  local meta_gm_k="" meta_projection="" meta_gm_lut="" meta_gm_table=""
  if [[ "$scheme" == "greedymini" ]]; then
    meta_gm_k="$p"; meta_projection="$projection"; meta_gm_lut="$gm_lut"
    meta_gm_table="$(gm_table_stem "$k")"
  fi

  META_RUN_ID="$rid" \
  META_DATASET="$ds" \
  META_LABEL="$label" \
  META_K="$k" \
  META_KIND="$kind" \
  META_THREADS="$THREADS" \
  META_MEM="$mem" \
  META_BINS="${BINS:-512}" \
  META_STRICT_MEM="$STRICT_MEM" \
  META_CONFIG_TAG="$cfg_tag" \
  META_SIGNATURE_LEN="$p" \
  META_SCHEME="$scheme" \
  META_SIG_MODE="$sig_mode" \
  META_GM_K="$meta_gm_k" \
  META_PROJECTION="$meta_projection" \
  META_GM_TABLE="$meta_gm_table" \
  META_GM_LUT="$meta_gm_lut" \
  META_ORDER_LABEL="${order_label:-}" \
  META_SEED="${seed:-}" \
  META_SEED_IDX="$seed_idx" \
  META_CUSTOM_ORDER="${custom_order:-}" \
  META_KMC_PATH="$KMC_BIN" \
  META_KMC_VERSION="$KMC_VERSION" \
  META_EXIT_CODE="$rc" \
  META_WALL="$wall" \
  META_PATH="$meta_path" \
  META_CMD_FILE="$cmd_file" \
  bench_py meta

  rm -f "$cmd_file"
  # kmc temp dir ($run_dir/tmp, on the results filesystem, up to ~32 GB): always removed after the run
  if [[ "$KEEP_TMP" -eq 0 && -d "$tmp" ]]; then
    rm -rf "$tmp"
  fi
  if [[ "$KEEP_DB" -eq 0 ]]; then
    rm -f "${db_prefix}.kmc_pre" "${db_prefix}.kmc_suf" "${db_prefix}.kmc_co" "${db_prefix}.kmc_col"
  fi

  [[ "$rc" -eq 0 ]] || die "kmc failed for $rid (exit $rc)"
  printf 'finished %s in %.1fs\n' "$rid" "$wall"
}

run_schemes() {
  local n="${#BENCH_SEEDS[@]}"
  [[ "$N_SEEDS" -ge 1 && "$N_SEEDS" -le "$n" ]] || die "-n must be between 1 and $n"

  IFS=',' read -r -a ds_arr <<<"$DATASETS"
  IFS=',' read -r -a k_arr <<<"$KS"
  IFS=',' read -r -a in_arr <<<"$INPUTS"

  local ds kind k i seed order_label

  case "$SCHEME" in
    default)
      # -n repeats the default scheme the way the seeded schemes repeat seeds:
      # it has no seed, but repeating it gives it a mean/std.
      [[ -z "$CUSTOM_ORDER" && -z "$ORDER_LABEL" ]] || die "--custom-order/--order-label require -S xor|xor-pcg-hash|greedymini|all"
      for ds in "${ds_arr[@]}"; do
        for kind in "${in_arr[@]}"; do
          for k in "${k_arr[@]}"; do
            for ((i = 0; i < N_SEEDS; i++)); do
              run_one "$(trim "$ds")" "$k" "$kind" default "" "" "" "$i"
            done
          done
        done
      done
      ;;
    lex)
      # True canonical lex, no is_allowed denylist (see README). Deterministic
      # (seed always 0, no real randomness) -- repeated -n times purely for
      # wall-time noise stats, same rationale as -S default.
      [[ -z "$CUSTOM_ORDER" && -z "$ORDER_LABEL" ]] || die "--custom-order/--order-label require -S xor|xor-pcg-hash|greedymini|all"
      for ds in "${ds_arr[@]}"; do
        for kind in "${in_arr[@]}"; do
          for k in "${k_arr[@]}"; do
            for ((i = 0; i < N_SEEDS; i++)); do
              run_one "$(trim "$ds")" "$k" "$kind" lex "" "lex" "0" "$i"
            done
          done
        done
      done
      ;;
    xor|xor-pcg-hash)
      [[ -z "$CUSTOM_ORDER" ]] || die "--custom-order is only valid with -S greedymini or all"
      for ((i = 0; i < N_SEEDS; i++)); do
        seed="${BENCH_SEEDS[$i]}"
        for ds in "${ds_arr[@]}"; do
          for kind in "${in_arr[@]}"; do
            for k in "${k_arr[@]}"; do
              order_label="${ORDER_LABEL:-$SCHEME}_s${seed}"
              run_one "$(trim "$ds")" "$k" "$kind" "$SCHEME" "" "$order_label" "$seed" "$i"
            done
          done
        done
      done
      ;;
    greedymini)
      for ((i = 0; i < N_SEEDS; i++)); do
        seed="${BENCH_SEEDS[$i]}"
        for ds in "${ds_arr[@]}"; do
          for kind in "${in_arr[@]}"; do
            for k in "${k_arr[@]}"; do
              order_label="${ORDER_LABEL}_s${seed}"
              [[ -n "$ORDER_LABEL" ]] || order_label="$(gm_label_for_k "$k")_s${seed}"
              run_one "$(trim "$ds")" "$k" "$kind" greedymini "$CUSTOM_ORDER" "$order_label" "$seed" "$i"
            done
          done
        done
      done
      ;;
    all)
      [[ -z "$ORDER_LABEL" ]] || die "--order-label is not used with -S all"
      # Interleave default/xor/xor-pcg-hash/greedymini within each repetition `i`, so
      # time-of-day/machine drift spreads evenly over the schemes instead of
      # landing as a step change between them; default is repeated like the
      # others, so the ratio to it has an error bar (README: kmc3_bench.sh).
      for ds in "${ds_arr[@]}"; do
        for kind in "${in_arr[@]}"; do
          for k in "${k_arr[@]}"; do
            for ((i = 0; i < N_SEEDS; i++)); do
              seed="${BENCH_SEEDS[$i]}"
              # DEFAULT_REPS (env, optional): repeat the default baseline only that many times
              [[ "$i" -lt "${DEFAULT_REPS:-$N_SEEDS}" ]] && run_one "$(trim "$ds")" "$k" "$kind" default "" "" "" "$i"
              run_one "$(trim "$ds")" "$k" "$kind" xor "" "xor_s${seed}" "$seed" "$i"
              run_one "$(trim "$ds")" "$k" "$kind" xor-pcg-hash "" "xor-pcg-hash_s${seed}" "$seed" "$i"
              run_one "$(trim "$ds")" "$k" "$kind" greedymini "$CUSTOM_ORDER" "$(gm_label_for_k "$k")_s${seed}" "$seed" "$i"
            done
          done
        done
      done
      ;;
    *)
      die "unsupported -S $SCHEME (use default, lex, xor, xor-pcg-hash, greedymini, or all)"
      ;;
  esac
}

export_results() {
  bench_py export \
    --results-dir "$RESULTS_DIR" \
    --paper-csv "$PAPER_CSV" \
    --paper-platform "$PAPER_PLATFORM"
}
