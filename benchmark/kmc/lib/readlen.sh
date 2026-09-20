# Read-length confound check: build size/length-matched derived datasets from an
# already-downloaded source, then benchmark both through the normal harness.
# See README: "Read-length confound check".

readlen_summary_field() {
  # Print a lex-row field for dataset/k from an existing kmc3_summary.csv (nothing if absent).
  local dataset="$1" k="$2" field="$3"
  local csv="$RESULTS_DIR/kmc3_summary.csv"
  [[ -f "$csv" ]] || return 1
  awk -F, -v ds="$dataset" -v k="$k" -v field="$field" '
    NR == 1 {
      for (i = 1; i <= NF; i++) h[$i] = i
      next
    }
    $h["dataset"] == ds && $h["k"] == k && $h["order_label"] == "default" {
      print $h[field]
      exit
    }
  ' "$csv"
}

readlen_target_kmers() {
  # Explicit --rl-target-kmers wins; else look up --rl-against's default (lex) n_total_kmers.
  # "default" specifically (not "lexicographic", which is the narrower true-canonical-lex
  # comparison scheme, not run everywhere) -- it's the one scheme every benchmark run has.
  if [[ -n "$RL_TARGET_KMERS" ]]; then
    echo "$RL_TARGET_KMERS"
    return 0
  fi
  [[ -n "$RL_AGAINST" ]] || die "pass --rl-target-kmers or --rl-against"
  local k="${KS%%,*}" v
  v="$(readlen_summary_field "$RL_AGAINST" "$k" n_total_kmers)"
  [[ -n "$v" ]] || die "can't find a default-scheme n_total_kmers for --rl-against $RL_AGAINST k=$k in $RESULTS_DIR/kmc3_summary.csv (run it first, e.g. 'run -d $RL_AGAINST -k $k'), or pass --rl-target-kmers explicitly"
  echo "${v%.*}"
}

readlen_build_dataset() {
  local out_id="$1" join="$2" target_kmers="$3" k="${KS%%,*}"
  local -a cmd=(
    make-readlen-dataset
    --source "$RL_SOURCE" --out-id "$out_id" --join "$join"
    --target-kmers "$target_kmers" --k "$k"
    --probe-reads "${RL_PROBE_READS:-200000}"
  )
  [[ "$FORCE" -eq 1 ]] && cmd+=(--force)
  bench_py "${cmd[@]}"
}

cmd_readlen_dataset() {
  [[ -n "$RL_SOURCE" ]] || die "--rl-source is required"
  [[ -n "$RL_OUT_ID" ]] || die "--rl-out-id is required"
  [[ -n "$RL_JOIN" ]] || die "--rl-join is required"
  local target
  target="$(readlen_target_kmers)"
  readlen_build_dataset "$RL_OUT_ID" "$RL_JOIN" "$target"
}

cmd_readlen_check() {
  [[ -n "$RL_SOURCE" ]] || die "--rl-source is required"
  [[ -n "$MEM_OVERRIDE" ]] || die "--mem is required (derived datasets aren't in data/datasets.csv)"
  local join="${RL_JOIN:-4}"
  local target
  target="$(readlen_target_kmers)"

  local short_id="${RL_SOURCE}_shortsub"
  local long_id="${RL_SOURCE}_longsub"

  echo "=== read-length confound check: target ~$target k-mers (matching ${RL_AGAINST:-explicit --rl-target-kmers}) ==="
  [[ "$N_SEEDS" -gt 1 ]] || echo "warning: -n $N_SEEDS repeat(s); pass -n 5 (or more) for a real wall_time_z, not just a single noisy sample" >&2

  echo "--- building $short_id: join=1, native read length, size-matched (isolates SIZE) ---"
  readlen_build_dataset "$short_id" 1 "$target"

  echo "--- building $long_id: join=$join, ~${join}x read length, size-matched (isolates LENGTH) ---"
  readlen_build_dataset "$long_id" "$join" "$target"

  local ds
  for ds in "$short_id" "$long_id"; do
    echo "=== benchmarking $ds ==="
    DATASETS="$ds"
    run_schemes
  done

  export_results
  echo "=== done: compare greedymini_dna's wall_time_z / wall_time_ratio for $short_id vs $long_id in the summary sheet ==="
}
