# sshash binary resolution, benchmark runs (build/check/bench/query +
# equivalence gate), and export dispatch.

# Elapsed seconds between two $EPOCHREALTIME samples ("sec.usec" strings) --
# pure bash, no subprocess (same helper as benchmark/kmc/lib/run.sh).
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

resolve_sshash() {
  if [[ -n "$SSHASH" ]]; then
    [[ -f "$SSHASH" ]] || die "sshash binary not found: $SSHASH"
    echo "$(cd "$(dirname "$SSHASH")" && pwd)/$(basename "$SSHASH")"
    return
  fi
  if [[ -x "$REPO_ROOT/build/sshash" ]]; then
    echo "$REPO_ROOT/build/sshash"
    return
  fi
  local w
  w="$(command -v sshash 2>/dev/null || true)"
  [[ -n "$w" ]] || die "sshash not found; pass --sshash PATH (after 'make -C $REPO_ROOT all')"
  echo "$w"
}

sshash_version() {
  local bin="$1"
  # `sshash` with no subcommand prints the "== SSHash ... (vX.Y.Z) ==" banner
  # to stdout (exit 0); every subcommand's own --help does not.
  "$bin" 2>&1 | grep -oE 'v[0-9]+\.[0-9]+\.[0-9]+' | head -1 || echo "unknown"
}

# Idempotently ensures $DATA_DIR/data/<ds>/queries/{pos,neg}.txt (+cached
# .fa) exist for dataset `$1`, generating them from a throwaway native-order
# reference build the first time (scripts/equivalence.py). Sets
# QUERY_POS_FASTA / QUERY_NEG_FASTA on return.
ensure_query_sets() {
  local ds="$1" k="$2" m="$3" fasta out_dir json
  fasta="$(dataset_fasta_path "$ds")"
  [[ -f "$fasta" ]] || die "missing $fasta; run: $0 download -d $ds"
  out_dir="$DATA_DIR/data/$ds/queries"
  json="$(bench_py equivalence ensure-queries \
    --sshash "$SSHASH_BIN" --fasta "$fasta" --k "$k" --m "$m" --threads "$THREADS" \
    --out-dir "$out_dir")"
  local pos_txt neg_txt
  pos_txt="$("$PYTHON" -c "import json,sys; print(json.loads(sys.argv[1])['pos'])" "$json")"
  neg_txt="$("$PYTHON" -c "import json,sys; print(json.loads(sys.argv[1])['neg'])" "$json")"
  QUERY_POS_FASTA="${pos_txt%.txt}.fa"
  QUERY_NEG_FASTA="${neg_txt%.txt}.fa"
  bench_py equivalence to-fasta --txt "$pos_txt" --fasta "$QUERY_POS_FASTA" >/dev/null
  bench_py equivalence to-fasta --txt "$neg_txt" --fasta "$QUERY_NEG_FASTA" >/dev/null
}

# runs/{dataset}/{modality}/{order}/{seed_idx} -- k/m aren't in the
# path (fixed for the whole benchmark, recorded in meta.json instead).
run_rel_path() {
  local ds="$1" order="$2" seed_idx="$3" modality="${4:-canonical}"
  echo "${ds}/${modality}/$(order_dir_name "$order")/${seed_idx}"
}

# sshash --gm-projection value a plain greedymini order runs with ("" for any
# other order), recorded in meta.json.
gm_projection_of() {
  case "$1" in
    greedymini-*) order_gm_split "$1"; echo "$GM_PROJECTION" ;;
    *) echo "" ;;
  esac
}

run_one() {
  local ds="$1" k="$2" m="$3" order="$4" gm_file="$5" seed="${6:-}" seed_idx="${7:-0}" modality="${8:-canonical}"

  local rid run_dir index_path meta_path
  rid="$(run_rel_path "$ds" "$order" "$seed_idx" "$modality")"
  run_dir="$RESULTS_DIR/runs/$rid"
  meta_path="$run_dir/meta.json"

  if [[ -f "$meta_path" && "$FORCE" -eq 0 ]]; then
    echo "skip existing run $rid"
    return 0
  fi

  mkdir -p "$run_dir"
  index_path="$run_dir/index.sshash"
  rm -f "$index_path"

  # Per-run tmp dir on the same SSD as $RESULTS_DIR, isolated from other runs.
  local tmp_dir="$run_dir/tmp"
  mkdir -p "$tmp_dir"

  local fasta
  fasta="$(dataset_fasta_path "$ds")"
  [[ -f "$fasta" ]] || die "missing $fasta; run: $0 download -d $ds"

  build_order_flags "$order" "$seed" "$gm_file" "$modality"

  local build_log="$run_dir/build.log" build_time="$run_dir/build_time.txt"
  local build_cmd=("$SSHASH_BIN" build -i "$fasta" -k "$k" -m "$m" -t "$THREADS" --verbose --bucket-stats "${ORDER_FLAGS[@]}" -d "$tmp_dir" -o "$index_path")
  local wrapped=("${build_cmd[@]}")
  [[ -x "$TIME_BIN" ]] && wrapped=("$TIME_BIN" -v -o "$build_time" "${build_cmd[@]}")

  local cmd_file
  cmd_file="$(mktemp)"
  {
    printf '%s\n' "${build_cmd[*]}"
  } >"$cmd_file"

  local t0 wall build_rc
  t0="$EPOCHREALTIME"
  set +e
  "${wrapped[@]}" >"$build_log" 2>&1
  build_rc=$?
  set -e
  wall="$(wall_time_delta "$t0" "$EPOCHREALTIME")"
  rm -rf "$tmp_dir"

  local check_rc=-1 check_ok=0 bench_rc=-1 query_rc=-1
  local index_bytes="" pos_match="" neg_match=""

  if [[ "$build_rc" -eq 0 && -f "$index_path" ]]; then
    index_bytes="$(stat -c%s "$index_path" 2>/dev/null || stat -f%z "$index_path")"

    local check_log="$run_dir/check.log"
    set +e
    "$SSHASH_BIN" check -i "$index_path" --verbose >"$check_log" 2>&1
    check_rc=$?
    set -e
    printf '%s\n' "$SSHASH_BIN check -i $index_path --verbose" >>"$cmd_file"
    grep -q "EVERYTHING OK!" "$check_log" && check_ok=1

    local bench_log="$run_dir/bench.log"
    set +e
    "$SSHASH_BIN" bench -i "$index_path" >"$bench_log" 2>&1
    bench_rc=$?
    set -e
    printf '%s\n' "$SSHASH_BIN bench -i $index_path" >>"$cmd_file"

    local qpath
    qpath="$(dataset_query_path "$ds")"
    if [[ -n "$qpath" && -f "$qpath" ]]; then
      local query_log="$run_dir/query.log" query_time="$run_dir/query_time.txt"
      local query_cmd=("$SSHASH_BIN" query -i "$index_path" -q "$qpath" --verbose)
      local wrapped_q=("${query_cmd[@]}")
      [[ -x "$TIME_BIN" ]] && wrapped_q=("$TIME_BIN" -v -o "$query_time" "${query_cmd[@]}")
      set +e
      "${wrapped_q[@]}" >"$query_log" 2>&1
      query_rc=$?
      set -e
      printf '%s\n' "${query_cmd[*]}" >>"$cmd_file"
    fi

    ensure_query_sets "$ds" "$k" "$m"
    local pos_json neg_json
    pos_json="$(bench_py equivalence check --sshash "$SSHASH_BIN" --index "$index_path" \
      --query-fasta "$QUERY_POS_FASTA" --expect positive)"
    neg_json="$(bench_py equivalence check --sshash "$SSHASH_BIN" --index "$index_path" \
      --query-fasta "$QUERY_NEG_FASTA" --expect negative)"
    printf '{"pos": %s, "neg": %s}\n' "$pos_json" "$neg_json" >"$run_dir/equivalence.json"
    pos_match="$("$PYTHON" -c "import json,sys; print('true' if json.loads(sys.argv[1])['match'] else 'false')" "$pos_json")"
    neg_match="$("$PYTHON" -c "import json,sys; print('true' if json.loads(sys.argv[1])['match'] else 'false')" "$neg_json")"

    if [[ "$KEEP_INDEX" -eq 0 ]]; then
      rm -f "$index_path"
    fi
  fi

  META_RUN_ID="$rid" \
  META_DATASET="$ds" \
  META_K="$k" \
  META_M="$m" \
  META_ORDER="$order" \
  META_MODALITY="$modality" \
  META_SEED="${seed:-}" \
  META_SEED_IDX="$seed_idx" \
  META_GM_FILE="${gm_file:-}" \
  META_GM_PROJECTION="$(gm_projection_of "$order")" \
  META_THREADS="$THREADS" \
  META_SSHASH_PATH="$SSHASH_BIN" \
  META_SSHASH_VERSION="$SSHASH_VERSION" \
  META_BUILD_EXIT_CODE="$build_rc" \
  META_CHECK_EXIT_CODE="$check_rc" \
  META_CHECK_OK="$check_ok" \
  META_BENCH_EXIT_CODE="$bench_rc" \
  META_QUERY_EXIT_CODE="$query_rc" \
  META_QUERY_POS_MATCH="$pos_match" \
  META_QUERY_NEG_MATCH="$neg_match" \
  META_INDEX_BYTES="$index_bytes" \
  META_WALL="$wall" \
  META_PATH="$meta_path" \
  META_CMD_FILE="$cmd_file" \
  bench_py meta

  [[ "$build_rc" -eq 0 ]] || die "sshash build failed for $rid (exit $build_rc); see $build_log"
  printf 'finished %s in %.1fs\n' "$rid" "$wall"
}

run_sweep() {
  local n="${#BENCH_SEEDS[@]}"
  [[ "$N_SEEDS" -ge 1 && "$N_SEEDS" -le "$n" ]] || die "-n must be between 1 and $n"

  IFS=',' read -r -a ds_arr <<<"$DATASETS"
  local orders_list
  orders_list="$(resolve_orders)"
  IFS=',' read -r -a order_arr <<<"$orders_list"

  # A plain greedymini order extends its gm_k table to the m-mer by the centre
  # extension, which needs m - gm_k even: fail up front with a clear message
  # instead of after the first (failing) sshash build.
  local o
  for o in "${order_arr[@]}"; do
    o="$(trim "$o")"
    [[ "$o" == greedymini-* ]] || continue
    bench_py orders check-gm "$o" "$M" || die "order $o is not valid for m=$M"
  done

  local ds order i seed gm_file
  for ds in "${ds_arr[@]}"; do
    ds="$(trim "$ds")"
    [[ -n "$ds" ]] || continue
    for ((i = 0; i < N_SEEDS; i++)); do
      seed="${BENCH_SEEDS[$i]}"
      for order in "${order_arr[@]}"; do
        order="$(trim "$order")"
        [[ -n "$order" ]] || continue
        gm_file="$(order_gm_file "$order")"
        run_one "$ds" "$K" "$M" "$order" "$gm_file" "$seed" "$i" "${MODALITY:-canonical}"
      done
    done
  done
}

export_results() {
  bench_py export --results-dir "$RESULTS_DIR"
}
