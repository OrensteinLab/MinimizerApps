#!/usr/bin/env bash
# download_library — NCBI taxonomy + standard libraries, then db/standard_p{10,0.5}.
#
# Same library set as this tree already has (archaea, bacteria, viral, human,
# UniVec_Core). Plasmid is omitted. Does not run `k2 build`; indexes are
# built per minimizer scheme by run_benchmark.sh.
#
# Usage: ./benchmark/kraken2/download_library.sh [-o BASE_DIR] [-f FRAC] [-j THREADS]
#                                         [-S] [-Z] [-F]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=benchmark_setup.sh
source "${SCRIPT_DIR}/benchmark_setup.sh"

CODEBASE="$(cd "${SCRIPT_DIR}/../../kraken2" && pwd)"
BASE_DIR="${KRAKEN2_BENCHMARK_DIR:-${SCRIPT_DIR}/benchmark_data}"
SAMPLE_FRAC="${SAMPLE_FRAC:-0.1}"
SAMPLE_SMOKE_FRAC="${SAMPLE_SMOKE_FRAC:-0.005}"
THREADS="${THREADS:-16}"
SKIP_SAMPLE=0
SMOKE_ONLY=0
FORCE_TAXONOMY="${FORCE_TAXONOMY:-0}"

show_help() {
    cat <<EOF
Usage: $0 [options]

Download NCBI taxonomy and the standard Kraken 2 libraries (see
scripts/standard_installation.sh). Then subsample into db/standard_p10 (10%) and
db/standard_p0.5 (0.5%) unless -S.

Options:
  -o BASE_DIR   Base directory (default: \$KRAKEN2_BENCHMARK_DIR or ./benchmark_data)
  -f FRAC       Fraction for standard_p10 (default: ${SAMPLE_FRAC})
  -j N          Threads for library download and subsample (default: ${THREADS})
  -S            Skip standard_p10 / standard_p0.5
  -Z            Only (re)create db/standard_p0.5 from existing db/standard
  -F            Force taxonomy refresh
  -h            Help
EOF
    exit 0
}

while getopts ":o:f:j:SZFh" opt; do
    case "${opt}" in
        o) BASE_DIR="${OPTARG}" ;;
        f) SAMPLE_FRAC="${OPTARG}" ;;
        j) THREADS="${OPTARG}" ;;
        S) SKIP_SAMPLE=1 ;;
        Z) SMOKE_ONLY=1 ;;
        F) FORCE_TAXONOMY=1 ;;
        h) show_help ;;
        \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
        :)  echo "Option -${OPTARG} requires an argument." >&2; exit 1 ;;
    esac
done

# shellcheck source=benchmark_interface.sh
source "${SCRIPT_DIR}/benchmark_interface.sh"

SHARED_DB="${SHARED_DB:-$(bench_shared_db_for standard)}"
SAMPLE_DB="${SAMPLE_DB:-$(bench_canonical_shared_db standard_p10)}"
SAMPLE_SMOKE_DB="${SAMPLE_SMOKE_DB:-$(bench_canonical_shared_db standard_p0.5)}"

mkdir -p "${SHARED_DB}" "${BASE_DIR}"

create_sample_db_at() {
    if (($# < 2)) || (($# % 2 != 0)); then
        echo "create_sample_db_at: expected FRAC DST pairs" >&2
        return 1
    fi
    local args=()
    while (($#)); do
        local frac="$1" dst="$2"
        shift 2
        echo "[*] Creating sample at ${dst} (keep ~${frac} of sequences)..."
        rm -rf "${dst}"
        mkdir -p "${dst}/library"
        ln -sfn "${SHARED_DB}/taxonomy" "${dst}/taxonomy"
        args+=("${frac}" "${dst}")
    done
    python3 "${SCRIPT_DIR}/create_sample_db.py" -j "${THREADS}" \
        "${SHARED_DB}" "${args[@]}"
    local i
    for ((i = 0; i < ${#args[@]}; i += 2)); do
        echo "[*] Sample ready: ${args[i+1]}"
    done
}

if [[ "${SMOKE_ONLY}" -eq 1 ]]; then
    if [[ ! -d "${SHARED_DB}/library" ]]; then
        echo "Error: missing ${SHARED_DB}/library (run $0 without -Z first)." >&2
        exit 1
    fi
    create_sample_db_at "${SAMPLE_SMOKE_FRAC}" "${SAMPLE_SMOKE_DB}"
    echo "standard_p0.5 ready: ${SAMPLE_SMOKE_DB}"
    exit 0
fi

taxonomy_needs_refresh() {
    local nodes="${SHARED_DB}/taxonomy/nodes.dmp"
    [[ ! -s "${nodes}" ]] && return 0
    local newest
    newest=$(find "${SHARED_DB}/library" -type f \
        \( -name 'assembly_summary.txt' -o -name 'library.fna' -o -name 'manifest.txt' \) \
        -printf '%T@\n' 2>/dev/null | sort -n | tail -1)
    [[ -z "${newest}" ]] && return 1
    local nodes_mtime
    nodes_mtime=$(stat -c '%Y' "${nodes}")
    newest=${newest%.*}
    (( newest > nodes_mtime ))
}

echo "========================================================="
echo " download_library  (standard_installation libraries)"
echo "========================================================="
echo "  Base:      ${BASE_DIR}"
echo "  Full DB:   ${SHARED_DB}"
echo "========================================================="

if [[ "${FORCE_TAXONOMY}" == "1" ]] || taxonomy_needs_refresh; then
    echo "[*] Downloading / refreshing NCBI taxonomy..."
    rm -rf "${SHARED_DB}/taxonomy"
    k2 download-taxonomy --db "${SHARED_DB}"
else
    echo "[*] Taxonomy up to date relative to libraries."
fi

mkdir -p "${SHARED_DB}/library"
# Match scripts/standard_installation.sh (nucleotide).
for lib in archaea bacteria viral plasmid human UniVec_Core; do
    if [[ ! -d "${SHARED_DB}/library/${lib}" ]]; then
        echo "[*] Downloading library: ${lib}..."
        extra_args=()
        if [[ "${lib}" == "human" ]]; then
            extra_args+=(--no-masking)
        fi
        k2 download-library --db "${SHARED_DB}" --library "${lib}" \
            --threads "${THREADS}" --resume "${extra_args[@]}"
    else
        echo "[*] Library '${lib}' already present. Skipping."
    fi
done

if taxonomy_needs_refresh; then
    echo "[*] Libraries newer than taxonomy; refreshing taxonomy again..."
    rm -rf "${SHARED_DB}/taxonomy"
    k2 download-taxonomy --db "${SHARED_DB}"
fi

echo "[*] Repairing taxonomy for taxids missing from NCBI taxdump..."
k2 repair-taxonomy --db "${SHARED_DB}"

if [[ "${SKIP_SAMPLE}" -eq 0 ]]; then
    create_sample_db_at \
        "${SAMPLE_FRAC}" "${SAMPLE_DB}" \
        "${SAMPLE_SMOKE_FRAC}" "${SAMPLE_SMOKE_DB}"
else
    echo "[*] Skipping db/standard_p10 and db/standard_p0.5 (-S)."
fi

echo ""
echo "download_library complete."
echo "  Full DB:   ${SHARED_DB}"
echo "  p10:       ${SAMPLE_DB}"
echo "  p0.5:      ${SAMPLE_SMOKE_DB}"
echo "Reads: ./benchmark/kraken2/download_data.sh -o ${BASE_DIR}"
