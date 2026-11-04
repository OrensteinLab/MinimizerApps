# sshash binary resolution, benchmark runs (build/check/bench/query +
# equivalence gate), and export dispatch.

# wait_low_load: block until the (shared) machine is quiet before each timed step
# (BENCH_LOAD_MAX; BENCH_WAIT_LOW_LOAD=0 disables; benchmark/common/wait_low_load.sh).
source "${BASH_SOURCE[0]%/*}/../../common/wait_low_load.sh"
# bench_pin_cpus: pin every timed sshash step to least-busy CPUs of one node, memory bound to it (benchmark/common/pin_cpus.sh;
# pinning is not optional): build $THREADS threads, bench and query one CPU; step_*.json record the memory state
source "${BASH_SOURCE[0]%/*}/../../common/pin_cpus.sh"

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

# Idempotently ensures $DATA_DIR/<ds>/queries/{pos,neg}.txt (+cached
# .fa) exist for dataset `$1`, generating them from a throwaway native-scheme
# reference build the first time (scripts/equivalence.py). Sets
# QUERY_POS_FASTA / QUERY_NEG_FASTA on return.
ensure_query_sets() {
  local ds="$1" k="$2" m="$3" fasta out_dir json
  fasta="$(dataset_fasta_path "$ds")"
  [[ -f "$fasta" ]] || die "missing $fasta; run: $0 download -d $ds"
  out_dir="$DATA_DIR/$ds/queries"
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

# runs/{dataset}/{modality}/{scheme}[/k<gm_k>_<projection>_t<mod_t>]/{seed_idx} -- k/m
# aren't in the path (fixed for the whole benchmark, recorded in meta.json instead).
# args: dataset scheme seed_idx modality gm_k projection mod_t
run_rel_path() {
  local ds="$1" scheme="$2" seed_idx="$3" modality="${4:-canonical}"
  echo "${ds}/${modality}/$(scheme_dir_name "$scheme" "${5:-}" "${6:-}" "${7:-}")/${seed_idx}"
}

# args: dataset k m scheme gm_file seed seed_idx modality projection gm_k gm_table mod_t
# (the last four are the variant's identity fields, empty when not applicable)
run_one() {
  local ds="$1" k="$2" m="$3" scheme="$4" gm_file="$5" seed="${6:-}" seed_idx="${7:-0}" modality="${8:-canonical}"
  local projection="${9:-}" gm_k="${10:-}" gm_table="${11:-}" mod_t="${12:-}"

  local rid run_dir index_path meta_path
  rid="$(run_rel_path "$ds" "$scheme" "$seed_idx" "$modality" "$gm_k" "$projection" "$mod_t")"
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

  build_scheme_flags "$scheme" "$seed" "$gm_file" "$modality" "$projection" "$mod_t"

  local build_log="$run_dir/build.log" build_time="$run_dir/build_time.txt"
  local build_cmd=("$SSHASH_BIN" build -i "$fasta" -k "$k" -m "$m" -t "$THREADS" --verbose --bucket-stats "${SCHEME_FLAGS[@]}" -d "$tmp_dir" -o "$index_path")
  local wrapped

  local cmd_file
  cmd_file="$(mktemp)"
  {
    printf '%s\n' "${build_cmd[*]}"
  } >"$cmd_file"

  local t0 wall build_rc
  wait_low_load 32 0 $((THREADS + 2))
  bench_pin_cpus "$THREADS" "$run_dir/step_build.json"
  printf '%s\n' "$BENCH_PIN_LIST" >"$run_dir/cpus.txt"
  wrapped=(${BENCH_PIN[@]+"${BENCH_PIN[@]}"} "${build_cmd[@]}")
  [[ -x "$TIME_BIN" ]] && wrapped=("$TIME_BIN" -v -o "$build_time" "${wrapped[@]}")
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
    # bench and query are single-threaded random lookups: one CPU (+1 spare core), and the node must
    # have room for the index in huge pages (their speed depends on the page size, see step records)
    wait_low_load 8 6 2
    BENCH_PIN_SPARE=1 bench_pin_cpus 1 "$run_dir/step_bench.json"
    printf '%s\n' "$BENCH_PIN_LIST" >"$run_dir/cpus_bench.txt"
    set +e
    ${BENCH_PIN[@]+"${BENCH_PIN[@]}"} "$SSHASH_BIN" bench -i "$index_path" >"$bench_log" 2>&1
    bench_rc=$?
    set -e
    printf '%s\n' "$SSHASH_BIN bench -i $index_path" >>"$cmd_file"

    local qpath
    qpath="$(dataset_query_path "$ds")"
    if [[ -n "$qpath" && -f "$qpath" ]]; then
      local query_log="$run_dir/query.log" query_time="$run_dir/query_time.txt"
      local query_cmd=("$SSHASH_BIN" query -i "$index_path" -q "$qpath" --verbose)
      local wrapped_q
      wait_low_load 8 6 2
      BENCH_PIN_SPARE=1 bench_pin_cpus 1 "$run_dir/step_query.json"
      printf '%s\n' "$BENCH_PIN_LIST" >"$run_dir/cpus_query.txt"
      wrapped_q=(${BENCH_PIN[@]+"${BENCH_PIN[@]}"} "${query_cmd[@]}")
      [[ -x "$TIME_BIN" ]] && wrapped_q=("$TIME_BIN" -v -o "$query_time" "${wrapped_q[@]}")
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
  META_SCHEME="$scheme" \
  META_MODALITY="$modality" \
  META_SEED="${seed:-}" \
  META_SEED_IDX="$seed_idx" \
  META_GM_FILE="${gm_file:-}" \
  META_GM_K="$gm_k" \
  META_PROJECTION="$projection" \
  META_GM_TABLE="$gm_table" \
  META_MOD_T="$mod_t" \
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
  local schemes_list
  schemes_list="$(resolve_schemes)"
  expand_scheme_variants "$schemes_list" "$K" "$M"

  # A plain greedymini variant extends its gm_k table to the m-mer by the center
  # extension, which needs m - gm_k even: fail up front with a clear message
  # instead of after the first (failing) sshash build.
  local v
  for v in "${SCHEME_VARIANTS[@]}"; do
    unpack_variant "$v"
    [[ "$V_SCHEME" == greedymini ]] || continue
    bench_py schemes check-gm "$V_SCHEME" "$M" "$V_GM_FILE" || die "scheme $V_SCHEME ($V_GM_TABLE) is not valid for m=$M"
  done

  local ds i seed
  for ds in "${ds_arr[@]}"; do
    ds="$(trim "$ds")"
    [[ -n "$ds" ]] || continue
    for ((i = 0; i < N_SEEDS; i++)); do
      seed="${BENCH_SEEDS[$i]}"
      for v in "${SCHEME_VARIANTS[@]}"; do
        unpack_variant "$v"
        run_one "$ds" "$K" "$M" "$V_SCHEME" "$V_GM_FILE" "$seed" "$i" "${MODALITY:-canonical}" \
          "$V_PROJECTION" "$V_GM_K" "$V_GM_TABLE" "$V_MOD_T"
      done
    done
  done
}

export_results() {
  bench_py export --results-dir "$RESULTS_DIR"
}
