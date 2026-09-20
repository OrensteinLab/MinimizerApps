#!/usr/bin/env bash
# prepare_strex_data — Wood et al. 2019 strain-exclusion *library* (Zenodo FASTA).
# "strex" = strain-exclusion library, not a standard_p* subsample.
#
# Classify reads (strain-exclusion + Table S5 ARGOS): ./benchmark/kraken2/download_data.sh
#
# Usage: ./benchmark/kraken2/prepare_strex_data.sh [-o BASE_DIR] [-j THREADS] [-F]

set -euo pipefail

EVAL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODEBASE="$(cd "${EVAL_DIR}/../../kraken2" && pwd)"
# shellcheck source=benchmark_setup.sh
source "${EVAL_DIR}/benchmark_setup.sh"

BASE_DIR="${KRAKEN2_BENCHMARK_DIR:-${EVAL_DIR}/benchmark_data}"
THREADS="${THREADS:-16}"
FORCE=0

show_help() {
    cat <<EOF
Usage: $0 [options]

Download Zenodo strain-exclusion FASTA and add it to db/strain_exclusion
(taxonomy + library.fna). Does not download classify reads.

Options:
  -o BASE_DIR   Base directory (default: \$KRAKEN2_BENCHMARK_DIR or ./benchmark_data)
  -j N          Threads for k2 add-to-library (default: ${THREADS})
  -F            Force re-download / re-extract / re-add
  -h            Help
EOF
    exit 0
}

while getopts ":o:j:Fh" opt; do
    case "${opt}" in
        o) BASE_DIR="${OPTARG}" ;;
        j) THREADS="${OPTARG}" ;;
        F) FORCE=1 ;;
        h) show_help ;;
        \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
        :)  echo "Option -${OPTARG} requires an argument." >&2; exit 1 ;;
    esac
done

# shellcheck source=benchmark_interface.sh
source "${EVAL_DIR}/benchmark_interface.sh"

bench_iface_init
SHARED="$(bench_shared_db_for strain_exclusion)"
# shellcheck source=zenodo_ref.sh
source "${EVAL_DIR}/zenodo_ref.sh"

mkdir -p "${SHARED}/taxonomy" "${SHARED}/library"

echo "========================================================="
echo " prepare_strex_data  (strain-exclusion library)"
echo "========================================================="
echo "  Base:   ${BASE_DIR}"
echo "  DB:     ${SHARED}"
echo "========================================================="

strex_ensure_zenodo 0

if [[ -z "${STREX_FASTA:-}" || -z "${NODES_DMP:-}" || -z "${NAMES_DMP:-}" ]]; then
    echo "Error: could not find strain_excluded.fna + nodes.dmp + names.dmp in ${ZENODO_EXTRACT}" >&2
    exit 1
fi

cp -f "${NODES_DMP}" "${SHARED}/taxonomy/nodes.dmp"
cp -f "${NAMES_DMP}" "${SHARED}/taxonomy/names.dmp"

if [[ ! -d "${SHARED}/library/added" || "${FORCE}" -eq 1 ]]; then
    echo "[*] k2 add-to-library ${STREX_FASTA}"
    rm -rf "${SHARED}/library/added"
    k2 add-to-library --db "${SHARED}" --threads "${THREADS}" --file "${STREX_FASTA}"
else
    echo "[*] Library already added under ${SHARED}/library/added"
fi

echo ""
echo "prepare_strex_data complete."
echo "  Exclusion DB: ${SHARED}"
echo "Build: ./benchmark/kraken2/run_benchmark.sh -E -o ${BASE_DIR}"
echo "NCBI libraries: ./benchmark/kraken2/download_library.sh -o ${BASE_DIR}"
