#!/usr/bin/env bash
# prepare_cami_data — CAMI II toy short-read sample -> data/cami/ (reads + truth).
#
#   data/cami/R1.fastq R2.fastq      ~250k pairs, every k-th pair of the sample
#   data/cami/truth.tsv              #seqid taxid domain   (CAMI read -> taxid)
#   data/cami/truth_extra.tsv        #seqid taxid genome_id novelty_category
#   data/cami/taxonomic_profile.txt  CAMI gold-standard profile of the full sample
#   data/cami/SOURCE.txt             provenance (URLs, md5, sample, stride)
#
# Source: CAMI II (Meyer et al. 2022, Nat Methods 19:429) toy datasets,
# PUBLISSO repository frl:6425521 (https://frl.publisso.de/data/frl:6425521/).
# Default sample: marine, sample 0 (2x150 bp, ~16.6M pairs, 977 genomes).
# Download ~6 GB (reads 5.2 GB + setup 0.9 GB); removed after sampling unless -k.
#
# Classify it with ./benchmark/kraken2/run_benchmark.sh -E (see README.md).
#
# Usage: ./benchmark/kraken2/prepare_cami_data.sh [-o BASE_DIR] [-s marine|strain|plant]
#                                                 [-n SAMPLE] [-p PAIRS] [-k] [-F]

set -euo pipefail

EVAL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=benchmark_setup.sh
source "${EVAL_DIR}/benchmark_setup.sh"

BASE_DIR="${KRAKEN2_BENCHMARK_DIR:-${EVAL_DIR}/benchmark_data}"
CAMI_SET="${CAMI_SET:-marine}"
CAMI_SAMPLE="${CAMI_SAMPLE:-0}"
PAIRS="${PAIRS:-250000}"
KEEP_RAW=0
FORCE=0
CAMI_URL="${CAMI_URL:-https://frl.publisso.de/data/frl:6425521}"

show_help() {
    cat <<EOF
Usage: $0 [options]

Download one CAMI II short-read sample, subsample it deterministically (every
k-th pair, original read IDs kept) and write truth.tsv from the CAMI
gold-standard read mapping.

Options:
  -o BASE_DIR   Base directory (default: \$KRAKEN2_BENCHMARK_DIR or ./benchmark_data)
  -s SET        marine (default; the only one validated), strain, plant
  -n N          Sample number within the dataset (default: ${CAMI_SAMPLE})
  -p N          Read pairs to keep (default: ${PAIRS})
  -k            Keep the downloaded tarballs (default: delete after sampling)
  -F            Force re-download / re-sample
  -h            Help
EOF
    exit 0
}

while getopts ":o:s:n:p:kFh" opt; do
    case "${opt}" in
        o) BASE_DIR="${OPTARG}" ;;
        s) CAMI_SET="${OPTARG}" ;;
        n) CAMI_SAMPLE="${OPTARG}" ;;
        p) PAIRS="${OPTARG}" ;;
        k) KEEP_RAW=1 ;;
        F) FORCE=1 ;;
        h) show_help ;;
        \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
        :)  echo "Option -${OPTARG} requires an argument." >&2; exit 1 ;;
    esac
done

case "${CAMI_SET}" in
    marine) SUBDIR="marine/short_read";          PFX="marmgCAMI2" ;;
    strain) SUBDIR="strain/short_read";          PFX="strmgCAMI2" ;;
    plant)  SUBDIR="plant_associated/short_read"; PFX="rhimgCAMI2" ;;
    *) echo "Error: -s must be marine, strain or plant." >&2; exit 1 ;;
esac
if ! [[ "${CAMI_SAMPLE}" =~ ^[0-9]+$ ]] || ! [[ "${PAIRS}" =~ ^[0-9]+$ ]]; then
    echo "Error: -n and -p must be non-negative integers." >&2
    exit 1
fi

DATA_DIR="${DATA_DIR:-${BASE_DIR}/data}"
OUT_DIR="${DATA_DIR}/cami"
RAW_DIR="${BASE_DIR}/downloads/cami_${CAMI_SET}"
READS_TGZ="${PFX}_sample_${CAMI_SAMPLE}_reads.tar.gz"
SETUP_TGZ="${PFX}_setup.tar.gz"

pick_nodes() {
    local cand
    for cand in \
        "${BASE_DIR}/db/strain_exclusion/taxonomy" \
        "${BASE_DIR}/db/standard/taxonomy" \
        "${BASE_DIR}/db/taxonomy"
    do
        if [[ -s "${cand}/nodes.dmp" ]]; then
            echo "${cand}/nodes.dmp"
            return 0
        fi
    done
    return 1
}

fetch() {  # name  -> RAW_DIR/name, resumable, md5-verified
    local name="$1" want got
    want="$(grep -E "  ${SUBDIR}/${name}\$" "${RAW_DIR}/md5sums.txt" | cut -d' ' -f1)"
    if [[ -z "${want}" ]]; then
        echo "Error: no md5 for ${SUBDIR}/${name} in md5sums.txt" >&2
        exit 1
    fi
    if [[ -f "${RAW_DIR}/${name}" ]]; then
        got="$(md5sum "${RAW_DIR}/${name}" | cut -d' ' -f1)"
        if [[ "${got}" == "${want}" ]]; then
            echo "[*] Already present (md5 ok): ${name}"
            return 0
        fi
    fi
    echo "[*] Downloading ${name} ..."
    curl -fSL -C - --retry 5 -o "${RAW_DIR}/${name}" "${CAMI_URL}/${SUBDIR}/${name}"
    got="$(md5sum "${RAW_DIR}/${name}" | cut -d' ' -f1)"
    if [[ "${got}" != "${want}" ]]; then
        echo "Error: md5 mismatch for ${name} (got ${got}, want ${want})" >&2
        exit 1
    fi
}

echo "========================================================="
echo " prepare_cami_data"
echo "========================================================="
echo "  Base:     ${BASE_DIR}"
echo "  CAMI II:  ${CAMI_SET} sample ${CAMI_SAMPLE}  (${CAMI_URL}/${SUBDIR})"
echo "  Output:   ${OUT_DIR}  (${PAIRS} pairs)"
echo "========================================================="

if [[ "${FORCE}" -eq 0 && -s "${OUT_DIR}/R1.fastq" && -s "${OUT_DIR}/R2.fastq" && -s "${OUT_DIR}/truth.tsv" ]]; then
    echo "[*] ${OUT_DIR} already prepared (-F to redo)."
    exit 0
fi

mkdir -p "${RAW_DIR}" "${OUT_DIR}"
EXTRACT="${RAW_DIR}/extract"
mkdir -p "${EXTRACT}"

curl -fsSL -o "${RAW_DIR}/md5sums.txt" "${CAMI_URL}/md5sums.txt"
fetch "${READS_TGZ}"
fetch "${SETUP_TGZ}"

echo "[*] Extracting reads, mapping and metadata ..."
tar -xzf "${RAW_DIR}/${READS_TGZ}" -C "${EXTRACT}"
tar -xzf "${RAW_DIR}/${SETUP_TGZ}" -C "${EXTRACT}" \
    "simulation_short_read/metadata.tsv" \
    "simulation_short_read/genome_to_id.tsv" \
    "simulation_short_read/taxonomic_profile_${CAMI_SAMPLE}.txt" \
    "simulation_short_read/config.ini"

SIM="${EXTRACT}/simulation_short_read"
READS_DIR="$(echo "${SIM}"/*_sample_"${CAMI_SAMPLE}"/reads)"
if [[ ! -s "${READS_DIR}/anonymous_reads.fq.gz" || ! -s "${READS_DIR}/reads_mapping.tsv.gz" ]]; then
    echo "Error: anonymous_reads.fq.gz / reads_mapping.tsv.gz not found under ${SIM}" >&2
    exit 1
fi

nodes_args=()
if NODES="$(pick_nodes)"; then
    nodes_args=(--nodes "${NODES}")
else
    echo "  Note: no nodes.dmp yet; truth.tsv domain column will be 'all' (scoring derives it)." >&2
fi

echo "[*] Sampling ${PAIRS} pairs + truth ..."
nice -n 19 python3 "${EVAL_DIR}/prepare_cami_reads.py" \
    --reads "${READS_DIR}/anonymous_reads.fq.gz" \
    --mapping "${READS_DIR}/reads_mapping.tsv.gz" \
    --metadata "${SIM}/metadata.tsv" \
    --out-dir "${OUT_DIR}" \
    --pairs "${PAIRS}" \
    "${nodes_args[@]}"

cp -f "${SIM}/taxonomic_profile_${CAMI_SAMPLE}.txt" "${OUT_DIR}/taxonomic_profile.txt"
{
    echo "dataset	CAMI II ${CAMI_SET} short read, sample ${CAMI_SAMPLE}"
    echo "source	${CAMI_URL}/${SUBDIR}/"
    echo "reads_tarball	${READS_TGZ}	md5 $(grep -E "  ${SUBDIR}/${READS_TGZ}\$" "${RAW_DIR}/md5sums.txt" | cut -d' ' -f1)"
    echo "setup_tarball	${SETUP_TGZ}	md5 $(grep -E "  ${SUBDIR}/${SETUP_TGZ}\$" "${RAW_DIR}/md5sums.txt" | cut -d' ' -f1)"
    echo "paper	Meyer et al. 2022, Nat Methods 19:429-440, doi:10.1038/s41592-022-01431-4"
    echo "subsample	every k-th read pair (see prepare_cami_reads.py), ${PAIRS} pairs, original IDs"
    echo "truth	reads_mapping.tsv.gz tax_id (CAMI NCBI taxdump 2017-02)"
} > "${OUT_DIR}/SOURCE.txt"

rm -rf "${EXTRACT}"
if [[ "${KEEP_RAW}" -eq 0 ]]; then
    rm -f "${RAW_DIR}/${READS_TGZ}" "${RAW_DIR}/${SETUP_TGZ}" "${RAW_DIR}/md5sums.txt"
    rmdir "${RAW_DIR}" 2>/dev/null || true
fi

echo ""
echo "prepare_cami_data complete."
echo "  Reads/truth: ${OUT_DIR}"
echo "Library coverage of the truth (excluded reads per library):"
echo "  python3 ${EVAL_DIR}/score_cami.py --nodes <db>/taxonomy/nodes.dmp \\"
echo "      --truth ${OUT_DIR}/truth.tsv --library-dir <db>"
echo "Run: ./benchmark/kraken2/run_benchmark.sh -E -o ${BASE_DIR}"
