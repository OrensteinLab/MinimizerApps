# kmc binary resolution, benchmark runs, and export dispatch

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
    echo "$root/data/$ds/files_gz.lst"
  else
    echo "$root/data/$ds/files_plain.lst"
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

scheme_dir_name() {
  case "$1" in
    default) echo "default" ;;
    lex) echo "lexicographic" ;;
    hash) echo "hash" ;;
    pcg) echo "pcg" ;;
    gm)
      case "${GM_LUT_MODE:-dna}" in
        binary) echo "greedymini_binary" ;;
        *) echo "greedymini_dna" ;;
      esac
      ;;
    *) die "unknown signature order: $1" ;;
  esac
}

run_rel_path() {
  local ds="$1" k="$2" kind="$3" order="$4" seed_idx="$5" mem="$6"
  echo "${ds}/$(params_tag "$k")/$(config_dir_name "$kind" "$mem")/$(scheme_dir_name "$order")/${seed_idx}"
}

run_one() {
  local ds="$1" k="$2" kind="$3" order="$4" custom_order="$5" order_label="$6" seed="${7:-}" seed_idx="${8:-0}"

  local mem label lst kmc_cmd wrapped p cfg_tag rid
  mem="$(dataset_mem "$ds")"
  label="$(csv_field "$ds" label)"
  p="$(signature_len_effective)"
  cfg_tag="$(config_dir_name "$kind" "$mem")"
  rid="$(run_rel_path "$ds" "$k" "$kind" "$order" "$seed_idx" "$mem")"

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
  elif [[ "$order" == "gm" ]]; then
    kmc_cmd+=("-p$(signature_len_effective)")
  fi
  if [[ "$order" == "gm" ]]; then
    local gm_lut="${GM_LUT_MODE:-dna}" gm_file="$custom_order" embedded
    [[ -n "$gm_file" ]] || gm_file="$(resolve_gm_path "$k" "${seed:-0}")"
    [[ -f "$gm_file" ]] || die "GM table not found for k=$k p=$(signature_len_effective) seed=${seed:-0} lut=$gm_lut: $gm_file"
    if [[ "$gm_lut" == "dna" ]] && embedded="$(bench_py gm dna-seed "$gm_file" 2>/dev/null)"; then
      seed="$embedded"
      kmc_cmd+=("-co$gm_file" "-ol$order_label")
    else
      [[ -n "$seed" ]] || die "GM binary .gm requires a seed; use a _dna.gm table or -n with predefined seeds"
      kmc_cmd+=("-co$gm_file" "-ol$order_label" "-os$((seed))")
      # -ode (DNA-expanded) is opt-in now that the binary rank LUT + runtime
      # projection is kmc's default backend for a binary-alphabet .gm source.
    fi
    custom_order="$gm_file"
  elif [[ "$order" == "lex" ]]; then
    # True canonical lex (no is_allowed denylist): -or hash -os 0 is order-
    # equivalent to MIN(i, revcomp(i)) with no exclusions (verified offline;
    # see README). Always seed 0 -- deterministic, not a real random scheme.
    kmc_cmd+=("-orhash" "-ol$order_label" "-os0")
  elif [[ "$order" == "hash" || "$order" == "pcg" ]]; then
    kmc_cmd+=("-or$order" "-ol$order_label" "-os$((seed))")
  fi
  kmc_cmd+=("@$lst" "$db_prefix" "$tmp")

  if [[ -x "$TIME_BIN" ]]; then
    wrapped=("$TIME_BIN" -v -o "$time_path" "${kmc_cmd[@]}")
  else
    echo "warning: $TIME_BIN missing; wall time from KMC JSON only" >&2
    wrapped=("${kmc_cmd[@]}")
  fi

  printf '%s\n' "${wrapped[*]}"
  local t0 wall rc cmd_file
  cmd_file="$(mktemp)"
  printf '%s\n' "${wrapped[@]}" >"$cmd_file"
  t0="$EPOCHREALTIME"
  set +e
  "${wrapped[@]}"
  rc=$?
  set -e
  wall="$(wall_time_delta "$t0" "$EPOCHREALTIME")"

  local sig_mode="" gm_lut="${GM_LUT_MODE:-}"
  if [[ -f "$json_path" ]]; then
    sig_mode="$("$PYTHON" -c "import json; print(json.load(open('$json_path')).get('signature_order_mode','') or '')")"
  fi
  if [[ -z "$sig_mode" ]]; then
    case "$order" in
      gm) sig_mode="table_dna" ;;
      hash) sig_mode="hash" ;;
      pcg) sig_mode="pcg" ;;
    esac
  fi
  if [[ "$order" == "gm" && -z "$gm_lut" ]]; then
    case "$sig_mode" in
      table_binary) gm_lut="binary" ;;
      table_dna|table) gm_lut="dna" ;;
      *) gm_lut="dna" ;;
    esac
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
  META_ORDER="$order" \
  META_SIG_MODE="$sig_mode" \
  META_GM_LUT_MODE="${gm_lut:-}" \
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

  case "$SIGNATURE_ORDER" in
    default)
      # -n repeats lex the same way custom schemes repeat seeds: it has no
      # real "seed", but repeating it gives it a mean/std instead of a
      # single unrepeatable sample (see README: lex baseline noise).
      [[ -z "$CUSTOM_ORDER" && -z "$ORDER_LABEL" ]] || die "--custom-order/--order-label require -O hash|pcg|gm|all"
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
      # wall-time noise stats, same rationale as -O default.
      [[ -z "$CUSTOM_ORDER" && -z "$ORDER_LABEL" ]] || die "--custom-order/--order-label require -O hash|pcg|gm|all"
      for ds in "${ds_arr[@]}"; do
        for kind in "${in_arr[@]}"; do
          for k in "${k_arr[@]}"; do
            for ((i = 0; i < N_SEEDS; i++)); do
              run_one "$(trim "$ds")" "$k" "$kind" lex "" "lexicographic" "0" "$i"
            done
          done
        done
      done
      ;;
    hash|pcg)
      [[ -n "$ORDER_LABEL" ]] || die "--order-label is required with -O hash|pcg|gm"
      [[ -z "$CUSTOM_ORDER" ]] || die "--custom-order is only valid with -O gm or all"
      for ((i = 0; i < N_SEEDS; i++)); do
        seed="${BENCH_SEEDS[$i]}"
        for ds in "${ds_arr[@]}"; do
          for kind in "${in_arr[@]}"; do
            for k in "${k_arr[@]}"; do
              order_label="${ORDER_LABEL}_s${seed}"
              run_one "$(trim "$ds")" "$k" "$kind" "$SIGNATURE_ORDER" "" "$order_label" "$seed" "$i"
            done
          done
        done
      done
      ;;
    gm)
      [[ -n "$ORDER_LABEL" ]] || ORDER_LABEL="gm"
      for ((i = 0; i < N_SEEDS; i++)); do
        seed="${BENCH_SEEDS[$i]}"
        for ds in "${ds_arr[@]}"; do
          for kind in "${in_arr[@]}"; do
            for k in "${k_arr[@]}"; do
              order_label="${ORDER_LABEL}_s${seed}"
              [[ "$ORDER_LABEL" == "gm" ]] && order_label="$(gm_label_for_k "$k")_s${seed}"
              run_one "$(trim "$ds")" "$k" "$kind" gm "$CUSTOM_ORDER" "$order_label" "$seed" "$i"
            done
          done
        done
      done
      ;;
    all)
      [[ -z "$ORDER_LABEL" ]] || die "--order-label is not used with -O all"
      # Interleave lex/hash/pcg/gm within each repetition `i` instead of
      # running all of one scheme, then all of the next: back-to-back
      # scheme blocks let any time-of-day/machine-state drift land as a
      # step change between schemes, which reads as a fake ordering
      # effect. Interleaving spreads that drift evenly across all four
      # schemes, and repeating lex (instead of running it once, first)
      # gives it the same mean/std as the others so the /lex ratio has an
      # actual error bar. See README: "Run order and the lex baseline".
      for ds in "${ds_arr[@]}"; do
        for kind in "${in_arr[@]}"; do
          for k in "${k_arr[@]}"; do
            for ((i = 0; i < N_SEEDS; i++)); do
              seed="${BENCH_SEEDS[$i]}"
              run_one "$(trim "$ds")" "$k" "$kind" default "" "" "" "$i"
              run_one "$(trim "$ds")" "$k" "$kind" hash "" "rand_hash_s${seed}" "$seed" "$i"
              run_one "$(trim "$ds")" "$k" "$kind" pcg "" "rand_pcg_s${seed}" "$seed" "$i"
              run_one "$(trim "$ds")" "$k" "$kind" gm "$CUSTOM_ORDER" "$(gm_label_for_k "$k")_s${seed}" "$seed" "$i"
            done
          done
        done
      done
      ;;
    *)
      die "unsupported -O $SIGNATURE_ORDER (use default, lex, hash, pcg, gm, or all)"
      ;;
  esac
}

export_results() {
  bench_py export \
    --results-dir "$RESULTS_DIR" \
    --paper-csv "$PAPER_CSV" \
    --paper-platform "$PAPER_PLATFORM"
}
