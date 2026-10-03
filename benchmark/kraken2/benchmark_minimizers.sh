#!/usr/bin/env bash
# benchmark_minimizers — time MinimizerScanner::NextMinimizer for selected
# schemes (same -S / -m / -z / -g / -G / -M interface as run_benchmark.sh).
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
# Density microbench uses NCBI library FASTAs, not the paper exclusion DB.
USE_FULL_NCBI=1

L="${L:-31}"
K="${K:-}"
N="${N:-100000000}"
BATCH="${BATCH:-1024}"
PARALLEL_JOBS="${PARALLEL_JOBS:-0}"
TOGGLE="${TOGGLE:-0xe37e28c4271b5a2d}"
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

Options (shared with run_benchmark.sh):
$(bench_common_options_help)
  -j N          Parallel density jobs (schemes × library files; default: OMP auto)

Microbench-only:
  -l INT        minimizer length ℓ (default: 31). GreedyMini needs ℓ = gm_k + 2i
                (w5_k15_canonical: 15, 17, ..., 31); the centre extension then
                keeps the table's density factor (~1.467) at every such ℓ.
                See kraken2/minimizers/README.md for the sweep commands.
  -k INT        k-mer length (default: l+4 ⇒ w=5)
  -N INT        approx DNA length for expected density (default: ${N})
  -b INT        batch size for latency samples (default: ${BATCH})
  -T HEX        toggle mask (default: ${TOGGLE})
  -s INT        minimizer-spaces s (default: ${SPACES}; kraken2-build NT is 7)
  -F FILE       classify timing file for lookup/resolve ns (or TIMING_FILE env)
  -L FLOAT      ns per TOK_LOOKUP override (or LOOKUP_NS env)
  -h, --help    Show help

Environment:
  KRAKEN2_BENCHMARK_DIR   Parent of db/ / data/ / runs/
  SHARED_DB               Override DB selected by -m/-z/-N/-E
  OUT_FILE                Density CSV (default: runs/<lib>/tables/minimizer_density.csv)
  PARALLEL_JOBS           Same as -j (0 = auto)
  GM_FILE / GM_DNA_FILE / MOD_GM_FILE
  TIMING_FILE / LOOKUP_NS / RESOLVE_NS
EOF
    exit 0
}

while getopts ":o:d:S:g:G:M:j:l:k:N:b:T:s:F:L:mzEh-" opt; do
    case "${opt}" in
        S) bench_handle_S_opt "$@" ;;
        j) PARALLEL_JOBS="${OPTARG}" ;;
        l) L="${OPTARG}" ;;
        k) K="${OPTARG}" ;;
        N) N="${OPTARG}" ;;
        b) BATCH="${OPTARG}" ;;
        T) TOGGLE="${OPTARG}" ;;
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

args=(-d "${DB_DIR}" -o "${OUT_FILE}" -l "${L}" -n "${N}" -b "${BATCH}"
      -T "${TOGGLE}" -s "${SPACES}" -j "${PARALLEL_JOBS}" -S "${DEFAULT_SCHEMES[@]}")
if [[ -n "${K}" ]]; then
    args+=(-k "${K}")
fi
if [[ "${NEED_GM}" -eq 1 ]]; then
    args+=(-g "${GM_FILE}")
fi
if [[ "${NEED_GM_DNA}" -eq 1 ]]; then
    args+=(-G "${GM_DNA_FILE}")
fi
if [[ "${NEED_MOD_GM}" -eq 1 ]]; then
    args+=(-M "${MOD_GM_FILE}")
fi
if [[ -n "${TIMING_FILE}" ]]; then
    args+=(-F "${TIMING_FILE}")
fi
if [[ -n "${LOOKUP_NS}" ]]; then
    args+=(-L "${LOOKUP_NS}")
fi
if [[ -n "${RESOLVE_NS}" ]]; then
    args+=(-R "${RESOLVE_NS}")
fi

echo "========================================"
echo "[*] DB:      ${DB_DIR} (${RUN_LAYOUT})"
echo "[*] Schemes: ${DEFAULT_SCHEMES[*]}"
echo "[*] Jobs:    ${PARALLEL_JOBS} (0 = OMP auto)"
echo "[*] ${BIN} ${args[*]}"
echo "========================================"
echo ""

"${BIN}" "${args[@]}"
status=$?

if [[ ${status} -eq 0 && -f "${OUT_FILE}" ]]; then
    echo ""
    echo "========================================"
    echo "[*] Density results (${OUT_FILE})"
    echo "========================================"
    cat "${OUT_FILE}"
fi
exit "${status}"
