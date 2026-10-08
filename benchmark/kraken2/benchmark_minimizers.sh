#!/usr/bin/env bash
# benchmark_minimizers — time MinimizerScanner::NextMinimizer for selected
# schemes (same -S / -m / -z / -g / -M interface as run_benchmark.sh).
# Reports dens_theory/expected/particular as CSV.
# Build via: make benchmark_minimizers (from repo root).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=benchmark_setup.sh
source "${SCRIPT_DIR}/benchmark_setup.sh"
# shellcheck source=benchmark_interface.sh
source "${SCRIPT_DIR}/benchmark_interface.sh"

CODEBASE="$(cd "${SCRIPT_DIR}/../../kraken2" && pwd)"
SRC_DIR="${CODEBASE}/src"
BIN="${SRC_DIR}/benchmark_minimizers"
bench_iface_init
# shellcheck source=bench_seeds.sh
source "${SCRIPT_DIR}/bench_seeds.sh"
# Density microbench uses NCBI library FASTAs, not the paper exclusion DB.
USE_FULL_NCBI=1

L="${L:-31}"
K="${K:-}"
N="${N:-100000000}"
BATCH="${BATCH:-1024}"
PARALLEL_JOBS="${PARALLEL_JOBS:-0}"
TOGGLE="${TOGGLE:-0xe37e28c4271b5a2d}"
TOGGLE_SET=0
N_SEEDS="${N_SEEDS:-1}"
FORCE=0
PYTHON="${PYTHON:-python3}"
SPACES="${SPACES:-7}"
TIMING_FILE="${TIMING_FILE:-}"
LOOKUP_NS="${LOOKUP_NS:-}"
RESOLVE_NS="${RESOLVE_NS:-}"
OUT_FILE="${OUT_FILE:-}"

show_help() {
    cat <<EOF
Usage: $0 [options]

Time NextMinimizer for selected schemes and report density + density factor
(= density × (w+1)) and super-mer mean/max lengths:
  theory, expected (random DNA), particular (DB library FASTAs).

Prefer: make benchmark_minimizers   (compiles with bench macros, then runs)

Per-seed mode (default): for seed i = 1..N (-n) the toggle mask is the same
one run_benchmark.sh builds the seed-i databases with (bench_const.py
DEFAULT_SEEDS[i-1], binary string -> 0x hex), every selected scheme (-S) is run
with --seed-idx i at -l/-k/-s, and its rows are MERGED into OUT_FILE
(rows with the same scheme + seed_idx are replaced, others kept). Schemes whose
row for (seed, toggle, l, k, s) already exists are skipped unless -f.
With an explicit -T HEX a single run is made (seed_idx 0, merged the same way).

Options (shared with run_benchmark.sh):
$(bench_common_options_help)
  -j N          Parallel density jobs (schemes × library files; default: OMP auto)

Microbench-only:
  -l INT        minimizer length ℓ (default: 31). greedymini needs ℓ = gm_k + 2i
                (w5_k15_canonical: 15, 17, ..., 31); the centre extension then
                keeps the table's density factor (~1.467) at every such ℓ.
                See kraken2/minimizers/README.md for the sweep commands.
  -k INT        k-mer length (default: l+4 ⇒ w=5)
  -N INT        approx DNA length for expected density (default: ${N})
  -b INT        batch size for latency samples (default: ${BATCH})
  -n N          number of seeds i = 1..N (default: ${N_SEEDS}; seed 1 = the
                default toggle). Same list as run_benchmark.sh -n
  -T HEX        single run with this toggle mask (seed_idx 0); excludes -n
  -f            recompute rows that are already in OUT_FILE
  -s INT        minimizer-spaces s (default: ${SPACES}; kraken2-build NT is 7)
  -F FILE       classify timing file for lookup/resolve ns (or TIMING_FILE env)
  -L FLOAT      ns per TOK_LOOKUP override (or LOOKUP_NS env)
  -h, --help    Show help

Environment:
  KRAKEN2_BENCHMARK_DIR   Parent of db/ / data/ / runs/
  SHARED_DB               Override DB selected by -m/-z/-N/-E
  OUT_FILE                Density CSV (default: runs/<lib>/tables/minimizer_density.csv)
  PARALLEL_JOBS           Same as -j (0 = auto)
  GM_FILE / MOD_GM_FILE
  TIMING_FILE / LOOKUP_NS / RESOLVE_NS
EOF
    exit 0
}

while getopts ":o:d:S:g:M:j:l:k:N:b:T:s:F:L:n:fmzEh-" opt; do
    case "${opt}" in
        S) bench_handle_S_opt "$@" ;;
        j) PARALLEL_JOBS="${OPTARG}" ;;
        l) L="${OPTARG}" ;;
        k) K="${OPTARG}" ;;
        N) N="${OPTARG}" ;;
        b) BATCH="${OPTARG}" ;;
        T) TOGGLE="${OPTARG}"; TOGGLE_SET=1 ;;
        n) N_SEEDS="${OPTARG}" ;;
        f) FORCE=1 ;;
        s) SPACES="${OPTARG}" ;;
        F) TIMING_FILE="${OPTARG}" ;;
        L) LOOKUP_NS="${OPTARG}" ;;
        h) show_help ;;
        -)
            case "${OPTARG}" in
                help) show_help ;;
                *) echo "Unknown option --${OPTARG}" >&2; exit 1 ;;
            esac
            ;;
        \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
        :)  echo "Option -${OPTARG} requires an argument." >&2; exit 1 ;;
        *)
            if ! bench_handle_common_opt "${opt}"; then
                echo "Invalid option: -${opt}" >&2
                exit 1
            fi
            ;;
    esac
done

shift $((OPTIND - 1)) || true

bench_validate_scheme_filter
[[ "${1:-}" == "-h" || "${1:-}" == "--help" ]] && show_help

if ! [[ "${PARALLEL_JOBS}" =~ ^[0-9]+$ ]]; then
    echo "Error: -j / PARALLEL_JOBS must be a non-negative integer." >&2
    exit 1
fi

bench_resolve_shared_db
DB_DIR="${SHARED_DB}"
if [[ -z "${OUT_FILE}" ]]; then
    OUT_FILE="$(bench_run_root)/tables/minimizer_density.csv"
fi

bench_require_gm_files

if [[ ! -x "${BIN}" ]]; then
    echo "Error: missing ${BIN}. Run: make benchmark_minimizers" >&2
    exit 1
fi
if [[ ! -d "${DB_DIR}" ]]; then
    echo "Error: database directory not found: ${DB_DIR}" >&2
    echo "  Set KRAKEN2_BENCHMARK_DIR and run download_library.sh (use -m/-z for sample DBs)." >&2
    exit 1
fi
if [[ ! -d "${DB_DIR}/library" ]]; then
    echo "Error: missing ${DB_DIR}/library" >&2
    echo "  Run ./benchmark/kraken2/download_library.sh -o \"\${KRAKEN2_BENCHMARK_DIR}\"." >&2
    exit 1
fi
if [[ -z "$(find "${DB_DIR}/library" -type f \( -name '*.fna' -o -name '*.fa' -o -name '*.fasta' -o -name '*.faa' \) -print -quit)" ]]; then
    echo "Error: no library FASTA under ${DB_DIR}/library" >&2
    echo "  Run ./benchmark/kraken2/download_library.sh -o \"\${KRAKEN2_BENCHMARK_DIR}\"." >&2
    exit 1
fi

mkdir -p "$(dirname "${OUT_FILE}")"

if [[ -n "${SEEDS_STR:-}" ]]; then
    read -r -a SEEDS <<< "${SEEDS_STR}"
fi
if ! [[ "${N_SEEDS}" =~ ^[0-9]+$ ]] || [[ "${N_SEEDS}" -lt 1 ]]; then
    echo "Error: -n / N_SEEDS must be a positive integer." >&2
    exit 1
fi
if [[ "${TOGGLE_SET}" -eq 0 && "${#SEEDS[@]}" -lt "${N_SEEDS}" ]]; then
    echo "Error: need at least ${N_SEEDS} seeds; have ${#SEEDS[@]}." >&2
    exit 1
fi

# K empty => l+4 in the tool; the CSV check needs the number.
K_EFF="${K:-$((L + 4))}"

echo "========================================"
echo "[*] DB:      ${DB_DIR} (${RUN_LAYOUT})"
echo "[*] Schemes: ${DEFAULT_SCHEMES[*]}"
echo "[*] l=${L} k=${K_EFF} s=${SPACES}  expected-DNA n=${N}"
echo "[*] Jobs:    ${PARALLEL_JOBS} (0 = OMP auto)"
echo "[*] Output:  ${OUT_FILE}"
echo "========================================"

# One (seed_idx, toggle) run for the schemes in $@ (names); merges into OUT_FILE.
run_one() {
    local seed_idx="$1" toggle="$2"
    shift 2
    local -a schemes=("$@")
    local tmp status=0
    tmp="$(mktemp "${OUT_FILE}.XXXXXX.csv")"
    local -a args=(-d "${DB_DIR}" -o "${tmp}" -l "${L}" -k "${K_EFF}" -n "${N}"
                   -b "${BATCH}" -T "${toggle}" -s "${SPACES}"
                   -j "${PARALLEL_JOBS}" --seed-idx "${seed_idx}"
                   -S "${schemes[@]}")
    local sch
    for sch in "${schemes[@]}"; do
        case "${sch}" in
            greedymini|greedymini-at1) args+=(-g "${GM_FILE}") ;;
            mod-minimizer-greedymini|mod-minimizer-greedymini-at1) args+=(-M "${MOD_GM_FILE}") ;;
        esac
    done
    [[ -n "${TIMING_FILE}" ]] && args+=(-F "${TIMING_FILE}")
    [[ -n "${LOOKUP_NS}" ]] && args+=(-L "${LOOKUP_NS}")
    [[ -n "${RESOLVE_NS}" ]] && args+=(-R "${RESOLVE_NS}")
    echo "[*] seed ${seed_idx} toggle ${toggle}: ${schemes[*]}"
    echo "[*] ${BIN} ${args[*]}"
    "${BIN}" "${args[@]}" || status=$?
    if [[ ${status} -eq 0 ]]; then
        "${PYTHON}" "${SCRIPT_DIR}/density_csv.py" merge "${OUT_FILE}" "${tmp}" "${SPACES}" || status=$?
    fi
    rm -f "${tmp}"
    return "${status}"
}

fail=0
if [[ "${TOGGLE_SET}" -eq 1 ]]; then
    seed_list=(0)
    toggles=("${TOGGLE}")
else
    seed_list=()
    toggles=()
    for ((i = 1; i <= N_SEEDS; i++)); do
        seed_list+=("${i}")
        # Bit-string (as k2 build --toggle-mask / build_db -T, strtoull base 2) -> hex.
        toggles+=("$(printf '0x%x' $((2#${SEEDS[$((i - 1))]})))")
    done
fi
for idx in "${!seed_list[@]}"; do
    sid="${seed_list[$idx]}"
    tg="${toggles[$idx]}"
    if [[ "${FORCE}" -eq 1 ]]; then
        todo=("${DEFAULT_SCHEMES[@]}")
    else
        mapfile -t todo < <("${PYTHON}" "${SCRIPT_DIR}/density_csv.py" missing \
            "${OUT_FILE}" "${sid}" "${tg}" "${L}" "${K_EFF}" "${SPACES}" "${DEFAULT_SCHEMES[@]}")
    fi
    if [[ ${#todo[@]} -eq 0 ]]; then
        echo "[*] seed ${sid}: density rows already in ${OUT_FILE}; skipping (-f to redo)."
        continue
    fi
    run_one "${sid}" "${tg}" "${todo[@]}" || { echo "Error: density run failed (seed ${sid})." >&2; fail=1; }
done

if [[ -f "${OUT_FILE}" ]]; then
    echo ""
    echo "========================================"
    echo "[*] Density results (${OUT_FILE})"
    echo "========================================"
    cat "${OUT_FILE}"
fi
exit "${fail}"
