#!/usr/bin/env bash
# KMC3 benchmark: download/prepare FASTQ, run kmc, export via `python -m scripts`.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export BENCHMARK_DIR="$HERE"

source "$HERE/lib/constants.sh"
source "$HERE/lib/env.sh"
load_bench_config

source "$HERE/lib/gm.sh"
source "$HERE/lib/datasets.sh"
source "$HERE/lib/run.sh"
source "$HERE/lib/experiment.sh"
source "$HERE/lib/readlen.sh"

# DATASETS_CSV and PAPER_CSV are set by load_bench_config from config.yaml.
PAPER_PLATFORM="xeon_hdd"

usage() {
  cat <<EOF
Usage: $(basename "$0") [command] [options]

Commands (default: run):
  download     Download gzip FASTQ from data/datasets.csv
  prepare      Decompress gzip FASTQ for the plain input arm
  run          Run kmc counts and export results
  export       Export CSV/XLSX only (\`python -m scripts export\`)
  all          download [+prepare if plain] + run + export
  experiment   Run a multi-config sweep from experiments/*.yaml
               Use 'experiment list' to show available experiments
  readlen-dataset  Build one size-matched derived dataset (see Read-length options)
  readlen-check    Build both size-matched variants, benchmark, export (see Read-length options)

Run options:
  -d, --datasets LIST     Comma-separated dataset ids (default: all five paper sets)
  -k, --k LIST            k values (default: 28,55)
  --input LIST            gzip,plain (default: gzip)
  -t, --threads N         kmc -t (default: 12)
  --mem N                 Override -m for all datasets (default: from datasets.csv)
  --bins N                kmc -n number of bins (64..2000; default: 512)
  --strict-mem            kmc -sm (hard RAM cap; stresses bin/memory path)
  --run-tag TAG           Replace auto config dir name (see Results layout)
  --kmc PATH              kmc binary (default: $REPO_ROOT/bin/kmc)
  -O, --signature-order S default|lex|hash|pcg|gm|all (default: all)
                          lex = true canonical lex, no is_allowed denylist
                          (standalone comparison scheme; not part of -O all)
  --custom-order PATH     GreedyMini .gm or _dna.gm (overrides per-k default w{k-p+1}_k{p})
  --order-label LABEL     Label for hash|pcg|gm (required unless -O all)
  -n, --n-seeds N         First N predefined seeds for custom schemes; also repeats the
                          default (lex) baseline N times (default: 1, max: ${#BENCH_SEEDS[@]})
  --signature-len N       Pass -pN to kmc
  --data-dir PATH         FASTQ root (default: KMC3_DATA_DIR or KMC3_BENCH_ROOT)
  --results-dir PATH      Results root (default: KMC3_RESULTS_DIR or .../results)
  --paper-csv PATH        Paper reference CSV
  --paper-platform P      xeon_hdd|xeon_ssd (default: xeon_hdd)
  --force                 Re-run / re-decompress even if outputs exist
  --keep-tmp              Keep kmc temp directories
  --keep-db               Keep .kmc_pre/.kmc_suf sidecars

Download options:
  --timeout SEC           curl connect timeout (default: 600)
  --skip-missing          Skip failed URLs; skip dataset if none succeed

Read-length options (readlen-dataset / readlen-check; see README "Read-length confound check"):
  --rl-source ID          Already-downloaded dataset id to build from (required)
  --rl-out-id ID          New dataset id to create (readlen-dataset only)
  --rl-join N             Consecutive source reads joined into 1 pseudo-read
                          (readlen-dataset: required; readlen-check: long variant, default 4;
                          the short variant is always join=1)
  --rl-target-kmers N     Target total k-mers (overrides --rl-against lookup)
  --rl-against ID         Size-match target: read its lex n_total_kmers from an existing
                          kmc3_summary.csv (default: f_vesca; run it first, or pass
                          --rl-target-kmers directly)
  --rl-probe-reads N      Source records sampled to estimate avg read length (default: 200000)
  (readlen-check also takes the Run options above: -k -t --mem --bins --strict-mem -n
  -O --force; --mem is required since derived datasets aren't in data/datasets.csv)

Examples:
  $0 download -d f_vesca
  $0 run -d f_vesca -k 28 -n 1
  $0 run -d f_vesca -k 28 -O all -n 5 --bins 1536 -t 24 --mem 8 --strict-mem --force
  $0 run -d f_vesca -k 28 -O lex -n 5 --bins 512 -t 12
  $0 experiment sweep_config
  $0 experiment sweep_gm_lut_fixed_dataset
  $0 experiment sweep_gm_lut_fixed_config
  $0 experiment list
  $0 export
  $0 readlen-dataset --rl-source h_sapiens_2 --rl-out-id h_sapiens_2_shortsub --rl-join 1 --rl-against f_vesca
  $0 readlen-check --rl-source h_sapiens_2 --rl-against f_vesca -k 28 -n 5 --mem 12 --bins 512 -t 12
EOF
}

parse_options() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      -d|--datasets) DATASETS="$2"; shift 2 ;;
      -k|--k) KS="$2"; shift 2 ;;
      --input) INPUTS="$2"; shift 2 ;;
      -t|--threads) THREADS="$2"; shift 2 ;;
      --mem) MEM_OVERRIDE="$2"; shift 2 ;;
      --bins) BINS="$2"; shift 2 ;;
      --strict-mem) STRICT_MEM=1; shift ;;
      --run-tag) RUN_TAG="$2"; shift 2 ;;
      --kmc) KMC="$2"; shift 2 ;;
      -O|--signature-order) SIGNATURE_ORDER="$2"; shift 2 ;;
      --custom-order) CUSTOM_ORDER="$2"; shift 2 ;;
      --order-label) ORDER_LABEL="$2"; shift 2 ;;
      -n|--n-seeds) N_SEEDS="$2"; shift 2 ;;
      --signature-len) SIGNATURE_LEN="$2"; shift 2 ;;
      --data-dir) DATA_DIR="$2"; shift 2 ;;
      --results-dir) RESULTS_DIR="$2"; shift 2 ;;
      --paper-csv) PAPER_CSV="$2"; shift 2 ;;
      --paper-platform) PAPER_PLATFORM="$2"; shift 2 ;;
      --timeout) DOWNLOAD_TIMEOUT="$2"; shift 2 ;;
      --skip-missing) SKIP_MISSING=1; shift ;;
      --rl-source) RL_SOURCE="$2"; shift 2 ;;
      --rl-out-id) RL_OUT_ID="$2"; shift 2 ;;
      --rl-join) RL_JOIN="$2"; shift 2 ;;
      --rl-target-kmers) RL_TARGET_KMERS="$2"; shift 2 ;;
      --rl-against) RL_AGAINST="$2"; shift 2 ;;
      --rl-probe-reads) RL_PROBE_READS="$2"; shift 2 ;;
      --force) FORCE=1; shift ;;
      --keep-tmp) KEEP_TMP=1; shift ;;
      --keep-db) KEEP_DB=1; shift ;;
      -h|--help) usage; exit 0 ;;
      *) die "unknown option: $1 (try -h)" ;;
    esac
  done
}

if [[ $# -gt 0 ]]; then
  case "$1" in
    download|prepare|run|export|all|experiment|readlen-dataset|readlen-check)
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

if [[ "$CMD" != "experiment" ]]; then
  normalize_signature_order
fi

load_dotenv
[[ -n "${KMC3_BENCH_ROOT:-}" ]] || die "KMC3_BENCH_ROOT is not set (see benchmark/.env)"

if [[ -z "$DATA_DIR" ]]; then
  DATA_DIR="${KMC3_DATA_DIR:-$KMC3_BENCH_ROOT}"
fi
if [[ -z "$RESULTS_DIR" ]]; then
  RESULTS_DIR="${KMC3_RESULTS_DIR:-$KMC3_BENCH_ROOT/results}"
fi

case "$CMD" in
  download)
    cmd_download
    ;;
  prepare)
    cmd_prepare
    ;;
  export)
    export_results
    ;;
  run)
    KMC_BIN="$(resolve_kmc)"
    KMC_VERSION="$(kmc_version "$KMC_BIN")"
    run_schemes
    export_results
    ;;
  all)
    cmd_download
    if inputs_include_plain; then
      cmd_prepare
    fi
    KMC_BIN="$(resolve_kmc)"
    KMC_VERSION="$(kmc_version "$KMC_BIN")"
    run_schemes
    export_results
    ;;
  readlen-dataset)
    cmd_readlen_dataset
    ;;
  readlen-check)
    KMC_BIN="$(resolve_kmc)"
    KMC_VERSION="$(kmc_version "$KMC_BIN")"
    cmd_readlen_check
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
