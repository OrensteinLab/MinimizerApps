#!/usr/bin/env bash
# Single-config runner: build index -> map reads -> score accuracy -> meta.json.
set -uo pipefail

# wait_low_load: block until the (shared) machine is quiet before each timed step
# (BENCH_LOAD_MAX, BENCH_WAIT_LOW_LOAD=0 disables; see benchmark/common/wait_low_load.sh).
source "$(dirname "${BASH_SOURCE[0]}")/../../common/wait_low_load.sh"
# bench_pin_cpus: every timed minimap2 step is pinned to one NUMA node, one thread per core, memory bound to it (benchmark/common/pin_cpus.sh)
source "$(dirname "${BASH_SOURCE[0]}")/../../common/pin_cpus.sh"

run_one() {
  local scheme="$1" seed="${2:-0}" seed_idx="${3:-0}"
  local preset k t
  if ! scheme_is_known "$scheme"; then
    echo "error: unknown scheme '$scheme' (schemes: ${ALL_SCHEMES//,/ }; greedymini's projection is PROJECTION=at0-cg1|at1-cg0)" >&2
    return 1
  fi
  preset=$(dataset_preset)
  k=${K:-$(preset_default_k "$preset")}
  t=${MOD_T:-$(( 4 + (k - 4) % W ))}
  if scheme_is_mod "$scheme" && { [ "$t" -le 0 ] || [ "$t" -ge "$k" ]; }; then
    echo "error: mod-minimizer needs 0 < t < k (k=$k w=$W gives t=$t); set MOD_T" >&2
    return 1
  fi
  # A GreedyMini table ranks length gm_k + 2i: the k-mer (greedymini) or the t-mer
  # (mod-minimizer-greedymini; w20_k9: k=19, w=10, t=9). Where mod-minimizer-greedymini's
  # table does not fit t (map-ont's k=15, t=5) the scheme does not apply and is skipped;
  # a greedymini table that does not fit k, or a table that cannot be read, is an error.
  local fits=0
  scheme_table_fits "$scheme" "$k" "$t" || fits=$?
  if [ "$fits" = 1 ] && [ "$scheme" = mod-minimizer-greedymini ]; then
    echo "== [$scheme] skip: its table $(scheme_gm_table "$scheme") does not fit t=$t (k=$k w=$W) ==" >&2
    return 0
  elif [ "$fits" != 0 ]; then
    echo "error: [$scheme] table $(scheme_gm_table "$scheme") does not fit k=$k t=$t, or cannot be read" >&2
    return 1
  fi
  local cond="k${k}_w${W}" map_extra=()
  if [ -n "${MID_OCC:-}" ]; then
    cond="${cond}_f${MID_OCC}"
    map_extra=(-f "$MID_OCC")
  fi
  # <scheme>[/k<gm_k>_<projection>_t<mod_t>] (scheme_names.scheme_dir_name): the
  # greedymini table k / projection and the mod-minimizer t are part of the
  # directory, so runs that differ in them never collide.
  local sdir gm_k="" proj="" gm_table=""
  sdir=$(scheme_dir_name "$scheme" "$t") || return 1
  if scheme_needs_gm "$scheme"; then
    gm_k=$(scheme_gm_k "$scheme") || return 1
    proj=$(scheme_projection) || return 1
    gm_table=$(scheme_gm_table "$scheme")
  fi
  local run_dir="$MINIMAP2_RESULTS_DIR/runs/$DATASET/$cond/$sdir/$seed_idx"
  local meta_json="$run_dir/meta.json"

  # meta.json is written last, only by a run that built, mapped and scored
  # cleanly, so its presence marks a complete run.
  if [ -f "$meta_json" ]; then
    if [ "${FORCE:-0}" != 1 ]; then
      echo "== [$sdir] skip: complete run in $run_dir (FORCE=1 reruns it) ==" >&2
      return 0
    fi
    # rerun: drop the old success marker first, so a failure below can't
    # leave the previous run's meta.json claiming success
    rm -f "$meta_json"
  fi
  mkdir -p "$run_dir"

  local data_dir="$MINIMAP2_DATA_DIR/$DATASET"
  local genome="$data_dir/genome.fa"
  local reads="$data_dir/reads.fa"
  local truth="$data_dir/truth.tsv"
  if [ ! -f "$genome" ] || [ ! -f "$reads" ]; then
    echo "error: dataset '$DATASET' not found under $data_dir (run: $0 make-dataset)" >&2
    return 1
  fi

  build_scheme_flags "$scheme" "$seed" "$t"

  local mmi="$run_dir/index.mmi"
  local build_time="$run_dir/build_time.txt"
  local build_log="$run_dir/build.log"
  local paf="$run_dir/map.paf"
  local map_time="$run_dir/map_time.txt"
  local map_log="$run_dir/map.log"
  local acc_json="$run_dir/accuracy.json"
  rm -f "$acc_json"

  # -x $preset (dataset_preset, lib/constants.sh): map-pb for the PacBio CLR
  # datasets -- the preset Winnowmap's paper (Jain et al. 2020) and
  # minimap2's docs recommend for them, which also enables HPC
  # (homopolymer-compressed) k-mers (`-H`) -- and map-ont for ONT reads (no
  # HPC). -k/-w passed *after* still win (confirmed via `-d` build log's own
  # "kmer size: $k; skip: $W" report) -- only k/w/H need to be at build
  # time (baked into the index); the other preset effects (scoring,
  # chaining bandwidth, etc.) are map-time parameters and must be repeated
  # on the map invocation too, where the loaded index's own k/w/H silently
  # win over the preset's (a harmless, expected stderr warning).
  wait_low_load 16 0 $((THREADS + 2))
  bench_pin_cpus "$THREADS" "$run_dir/step_build.json"
  printf '%s\n' "$BENCH_PIN_LIST" > "$run_dir/cpus_build.txt"
  echo "== [$sdir] building index (-x $preset k=$k w=$W ${SCHEME_FLAGS[*]:-}) ==" >&2
  if ! /usr/bin/time -v -o "$build_time" ${BENCH_PIN[@]+"${BENCH_PIN[@]}"} \
      "$MINIMAP2_BIN" -d "$mmi" -x "$preset" -k "$k" -w "$W" "${SCHEME_FLAGS[@]}" -t "$THREADS" "$genome" 2> "$build_log"; then
    echo "  build FAILED, see $build_log" >&2
    tail -5 "$build_log" >&2
    return 1
  fi

  wait_low_load 16 0 $((THREADS + 2))
  bench_pin_cpus "$THREADS" "$run_dir/step_map.json"
  printf '%s\n' "$BENCH_PIN_LIST" > "$run_dir/cpus_map.txt"
  echo "== [$sdir] mapping ${reads##*/} ==" >&2
  if ! /usr/bin/time -v -o "$map_time" ${BENCH_PIN[@]+"${BENCH_PIN[@]}"} \
      "$MINIMAP2_BIN" -x "$preset" -c "${map_extra[@]}" -t "$THREADS" "$mmi" "$reads" > "$paf" 2> "$map_log"; then
    echo "  map FAILED, see $map_log" >&2
    tail -5 "$map_log" >&2
    return 1
  fi

  # "distinct minimizers" is NOT the density_factor numerator -- it collapses
  # every repeated occurrence of the same minimizer key down to one count,
  # so it always undercounts total selected positions on a genome with any
  # repeat content (see mm_idx_stat in index.c: it separately tracks `sum`,
  # the true total-occurrences count, and prints it only via "average
  # occurrences"/"average spacing", not directly). Use "average spacing"
  # (genome_len / total_occurrences) for density_factor instead --
  # write_meta.py's own comment has the full reasoning.
  local density
  density=$(grep -o 'distinct minimizers: [0-9]*' "$build_log" | grep -o '[0-9]*' | head -1)
  local avg_occ
  avg_occ=$(grep -o 'average occurrences: [0-9.]*' "$build_log" | grep -o '[0-9.]*$' | head -1)
  local avg_spacing
  avg_spacing=$(grep -o 'average spacing: [0-9.]*' "$build_log" | grep -o '[0-9.]*$' | head -1)

  # Effective cap minimap2 actually used (its own default, or -f's value).
  local mid_occ
  mid_occ=$(grep -o 'mid_occ = [0-9]*' "$map_log" | grep -o '[0-9]*$' | head -1)

  # Real reads (no truth.tsv, e.g. the paper's D3): mapped fraction only.
  # A scoring failure fails the run: no meta.json, so export never sees it.
  if ! "$PYTHON" "$HERE_DIR/lib/score_accuracy.py" --paf "$paf" --truth "$truth" --reads "$reads" --json-out "$acc_json" >&2; then
    echo "  scoring FAILED for $run_dir (no meta.json written)" >&2
    return 1
  fi

  # Written to a temp name and renamed, so a half-written meta.json never
  # looks like a complete run.
  if ! "$PYTHON" "$HERE_DIR/lib/write_meta.py" \
      --scheme "$scheme" --gm-k "$gm_k" --projection "$proj" --gm-table "$gm_table" --seed "$seed" --seed-idx "$seed_idx" --k "$k" --w "$W" --mod-t "$t" --dataset "$DATASET" \
      --mid-occ-cap "${MID_OCC:-}" --mid-occ "${mid_occ:-}" --preset "$preset" --threads "$THREADS" \
      --build-time "$build_time" --map-time "$map_time" \
      --density "${density:-0}" --avg-occurrences "${avg_occ:-0}" --avg-spacing "${avg_spacing:-0}" \
      --accuracy "$acc_json" --out "$meta_json.tmp" >/dev/null; then
    echo "  write_meta FAILED for $run_dir" >&2
    rm -f "$meta_json.tmp"
    return 1
  fi
  # Indexes are 5-8 GB each; the scheme analysis reads only seed 0's of every variant, so the others
  # are deleted once mapped (their size is kept in index_bytes.txt, written before meta.json marks
  # the run complete). KEEP_ALL_INDEXES=1 keeps them.
  stat -c %s "$mmi" > "$run_dir/index_bytes.txt"
  mv -f "$meta_json.tmp" "$meta_json"
  if [ "$seed_idx" != 0 ] && [ "${KEEP_ALL_INDEXES:-0}" != 1 ]; then
    rm -f "$mmi"
  fi

  echo "  -> $meta_json"
}
