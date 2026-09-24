#!/usr/bin/env bash
# sshash benchmark: download datasets, run sshash build/check/bench/query
# across minimizer orders x seeds, export via `python -m scripts`.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export BENCHMARK_DIR="$HERE"

source "$HERE/lib/constants.sh"
source "$HERE/lib/env.sh"
load_bench_config

source "$HERE/lib/orders.sh"
source "$HERE/lib/datasets.sh"
source "$HERE/lib/run.sh"
source "$HERE/lib/experiment.sh"

# DATASETS_CSV and DEFAULT_GM_FILE are set by load_bench_config from config.yaml.

usage() {
  cat <<EOF
Usage: $(basename "$0") [command] [options]

Commands (default: run):
  download     Download eulertigs FASTA + query FASTQ from data/datasets.csv
  run          Run sshash build/check/bench/query across orders x seeds, export
  export       Export CSV/XLSX only (\`python -m scripts export\`)
  all          download + run + export
  experiment   Run a multi-config sweep from experiments/*.yaml
               Use 'experiment list' to show available experiments

Run options:
  -d, --datasets LIST     Comma-separated dataset ids (default: cod,kestrel,human)
  -k K                    k-mer length (default: 31)
  -m M                    Minimizer length, fixed across ALL orders by design
                          (default: 21, giving w=k-m+1=11 -- the mod-minimizer
                          paper's own tested SSHash regime; see README)
  -O, --orders LIST       native,xor-hash,xor-pcg-hash,greedymini,
                          mod-minimizer-native-hash,mod-minimizer-xor-hash,
                          mod-minimizer-pcg-xor-hash,mod-minimizer-greedymini,
                          or "all" (default: all). Every order (including
                          native) is run once per seed -- sshash has no
                          seed-less baseline the way kmc's lex does. Each
                          mod-minimizer-* order (except -greedymini, fixed at
                          t=11) additionally runs once per t in {10, 11}.
  -M, --modality VALUE   canonical (default, single-probe) or regular
                          (two-probe, sshash build --regular). Separate axis
                          from -O; a mod-minimizer-* order is always regular
                          regardless. E.g. -O native -M regular is the
                          apples-to-apples comparison against mod-minimizer's
                          own (mandatory) two-probe path.
  --greedy-mini-file PATH GreedyMini .gm table for -O greedymini
                          (default: orders/GreedyMini/w11_k15_canonical.gm;
                          mod-minimizer-greedymini's own O_t table is fixed to
                          orders/GreedyMini/w21_k11_canonical.gm, config.yaml
                          mod_t_gm_file, not overridable from here)
  -n, --n-seeds N         First N predefined seeds from config.yaml (default: 5, max: ${#BENCH_SEEDS[@]})
  -t, --threads N         sshash build/query -t (default: 8)
  --sshash PATH           sshash binary (default: $REPO_ROOT/build/sshash)
  --data-dir PATH         Dataset root (default: SSHASH_DATA_DIR or SSHASH_BENCH_ROOT)
  --results-dir PATH      Results root (default: SSHASH_RESULTS_DIR or .../results)
  --dataset-file PATH     Ad hoc local FASTA to benchmark under dataset id "adhoc"
                          (bypasses datasets.csv entirely -- for local smoke testing)
  --keep-index            Keep the built .sshash index after a run (default: removed
                          after check/bench/query/equivalence to save disk)
  --force                 Re-run even if a run's meta.json already exists

Download options:
  --timeout SEC           curl connect timeout (default: 600)
  --skip-missing          Skip failed URLs instead of failing the whole download

Examples:
  $0 download -d cod
  $0 run -d cod -n 1
  $0 run -d cod,kestrel,human -O all -n 5 -t 16
  $0 run -d cod,kestrel -O native -M regular -n 5 -t 16
  $0 run --dataset-file sshash/data/unitigs_stitched/salmonella_enterica_k31_ust.fa.gz -n 1
  $0 experiment sweep_orders
  $0 experiment sweep_m
  $0 experiment list
  $0 export
EOF
}

parse_options() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      -d|--datasets) DATASETS="$2"; shift 2 ;;
      -k) K="$2"; shift 2 ;;
      -m) M="$2"; shift 2 ;;
      -O|--orders) ORDERS="$2"; shift 2 ;;
      -M|--modality)
        MODALITY="$2"
        [[ "$MODALITY" == "canonical" || "$MODALITY" == "regular" ]] || \
          die "--modality must be 'canonical' or 'regular', got: $MODALITY"
        shift 2 ;;
      --greedy-mini-file) GM_FILE="$2"; shift 2 ;;
      -n|--n-seeds) N_SEEDS="$2"; shift 2 ;;
      -t|--threads) THREADS="$2"; shift 2 ;;
      --sshash) SSHASH="$2"; shift 2 ;;
      --data-dir) DATA_DIR="$2"; shift 2 ;;
      --results-dir) RESULTS_DIR="$2"; shift 2 ;;
      --dataset-file) DATASET_FILE="$2"; shift 2 ;;
      --run-tag) RUN_TAG="$2"; shift 2 ;;
      --timeout) DOWNLOAD_TIMEOUT="$2"; shift 2 ;;
      --skip-missing) SKIP_MISSING=1; shift ;;
      --keep-index) KEEP_INDEX=1; shift ;;
      --force) FORCE=1; shift ;;
      -h|--help) usage; exit 0 ;;
      *) die "unknown option: $1 (try -h)" ;;
    esac
  done
}

if [[ $# -gt 0 ]]; then
  case "$1" in
    download|run|export|all|experiment)
      CMD="$1"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
  esac
fi

if [[ "$CMD" == "experiment" ]]; then
  EXPERIMENT="${1:-list}"
  [[ $# -gt 0 ]] && shift
fi

parse_options "$@"

load_dotenv
[[ -n "${SSHASH_BENCH_ROOT:-}" ]] || die "SSHASH_BENCH_ROOT is not set (see benchmark/sshash/.env)"

if [[ -z "$DATA_DIR" ]]; then
  DATA_DIR="${SSHASH_DATA_DIR:-$SSHASH_BENCH_ROOT}"
fi
if [[ -z "$RESULTS_DIR" ]]; then
  RESULTS_DIR="${SSHASH_RESULTS_DIR:-$SSHASH_BENCH_ROOT/results}"
fi

# --dataset-file registers an "adhoc" dataset id pointing at an arbitrary
# local FASTA (no download, no query FASTQ) -- used for local smoke testing
# without touching data/datasets.csv (see README "Smoke testing").
if [[ -n "$DATASET_FILE" ]]; then
  [[ -f "$DATASET_FILE" ]] || die "--dataset-file not found: $DATASET_FILE"
  ADHOC_ABS="$(cd "$(dirname "$DATASET_FILE")" && pwd)/$(basename "$DATASET_FILE")"
  ADHOC_CSV="$(mktemp)"
  {
    printf 'dataset,species,eulertigs_filename,eulertigs_size,eulertigs_md5,eulertigs_url,query_filename,query_size,query_md5,query_url,local_path,local_query_path\n'
    cat "$DATASETS_CSV" | tail -n +2
    printf 'adhoc,ad hoc local FASTA,,,,,,,,,%s,\n' "$ADHOC_ABS"
  } >"$ADHOC_CSV"
  DATASETS_CSV="$ADHOC_CSV"
  DATASETS="adhoc"
fi

case "$CMD" in
  download)
    cmd_download
    ;;
  export)
    export_results
    ;;
  run)
    SSHASH_BIN="$(resolve_sshash)"
    SSHASH_VERSION="$(sshash_version "$SSHASH_BIN")"
    run_schemes
    export_results
    ;;
  all)
    cmd_download
    SSHASH_BIN="$(resolve_sshash)"
    SSHASH_VERSION="$(sshash_version "$SSHASH_BIN")"
    run_schemes
    export_results
    ;;
  experiment)
    if [[ "$EXPERIMENT" == "list" ]]; then
      list_experiments
    else
      run_experiment
    fi
    ;;
  *)
    die "unknown command: $CMD"
    ;;
esac
