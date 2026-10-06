#!/usr/bin/env bash
# prepare_mock_data — real mock-community shotgun run -> data/mock/ (reads + composition).
#
#   data/mock/R1.fastq R2.fastq           ~250k pairs, every k-th pair of the run
#   data/mock/composition.csv             expected composition (copy of
#                                         tables/mock_zymo_d6300_expected.csv)
#   data/mock/SOURCE.txt                  provenance (run, ENA URLs, stride)
#
# Mock: ZymoBIOMICS Microbial Community Standard (D6300; 8 bacteria + 2 yeasts),
# SRA SRR12324251 / BioProject PRJNA648136, Illumina NovaSeq 6000, 2x150 bp,
# ~29.4M pairs. Expected composition: Zymo D6300 product sheet, Table 1.
# Streams ~3.5 GB from ENA; only the sample is stored.
#
# Classify it with ./benchmark/kraken2/run_benchmark.sh -E (see README.md).
#
# Usage: ./benchmark/kraken2/prepare_mock_data.sh [-o BASE_DIR] [-r RUN] [-p PAIRS] [-F]

set -euo pipefail

EVAL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=benchmark_setup.sh
source "${EVAL_DIR}/benchmark_setup.sh"

BASE_DIR="${KRAKEN2_BENCHMARK_DIR:-${EVAL_DIR}/benchmark_data}"
MOCK_RUN="${MOCK_RUN:-SRR12324251}"
PAIRS="${PAIRS:-250000}"
FORCE=0
TABLES="${EVAL_DIR}/tables"
COMPOSITION="${TABLES}/mock_zymo_d6300_expected.csv"

show_help() {
    cat <<EOF
Usage: $0 [options]

Stream the ZymoBIOMICS D6300 shotgun run from ENA, keep every k-th read pair
(deterministic, no RNG) and copy the expected composition next to the reads.

Options:
  -o BASE_DIR   Base directory (default: \$KRAKEN2_BENCHMARK_DIR or ./benchmark_data)
  -r RUN        SRA/ENA run accession (default: ${MOCK_RUN}); the composition
                table is for D6300, so only use another D6300 run
  -p N          Read pairs to keep (default: ${PAIRS})
  -F            Force re-download
  -h            Help
EOF
    exit 0
}

while getopts ":o:r:p:Fh" opt; do
    case "${opt}" in
        o) BASE_DIR="${OPTARG}" ;;
        r) MOCK_RUN="${OPTARG}" ;;
        p) PAIRS="${OPTARG}" ;;
        F) FORCE=1 ;;
        h) show_help ;;
        \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
        :)  echo "Option -${OPTARG} requires an argument." >&2; exit 1 ;;
    esac
done

if ! [[ "${PAIRS}" =~ ^[0-9]+$ ]]; then
    echo "Error: -p must be a non-negative integer." >&2
    exit 1
fi

DATA_DIR="${DATA_DIR:-${BASE_DIR}/data}"
OUT_DIR="${DATA_DIR}/mock"

echo "========================================================="
echo " prepare_mock_data"
echo "========================================================="
echo "  Base:     ${BASE_DIR}"
echo "  Run:      ${MOCK_RUN} (ZymoBIOMICS D6300, PRJNA648136)"
echo "  Output:   ${OUT_DIR}  (${PAIRS} pairs)"
echo "========================================================="

mkdir -p "${OUT_DIR}"
if [[ "${FORCE}" -eq 0 && -s "${OUT_DIR}/R1.fastq" && -s "${OUT_DIR}/R2.fastq" ]]; then
    echo "[*] Reads already present in ${OUT_DIR} (-F to redo)."
else
    nice -n 19 python3 "${EVAL_DIR}/prepare_mock_reads.py" \
        --run "${MOCK_RUN}" --out-dir "${OUT_DIR}" --pairs "${PAIRS}"
fi
cp -f "${COMPOSITION}" "${OUT_DIR}/composition.csv"

echo ""
echo "prepare_mock_data complete."
echo "  Reads/composition: ${OUT_DIR}"
echo "Run: ./benchmark/kraken2/run_benchmark.sh -E -o ${BASE_DIR}"
