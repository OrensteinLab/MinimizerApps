#!/usr/bin/env bash
# download_data — fetch both classify datasets under data/.
#
#   data/strain-exclusion/   Mason / Zenodo accuracy FASTQs + truth.tsv
#   data/fda-argos/          Table S5 SRR runs (benchmark/kraken2/tables/table_s5.csv)
#
# NCBI taxonomy/libraries: ./benchmark/kraken2/download_library.sh
# Strain-exclusion FASTA:  ./benchmark/kraken2/prepare_strex_data.sh
#
# Usage: ./benchmark/kraken2/download_data.sh [-o BASE_DIR] [-p PAIRS] [-F] [-s] [-a]

set -euo pipefail

EVAL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=benchmark_setup.sh
source "${EVAL_DIR}/benchmark_setup.sh"

CODEBASE="$(cd "${EVAL_DIR}/../../kraken2" && pwd)"
BASE_DIR="${KRAKEN2_BENCH_ROOT:-${EVAL_DIR}/benchmark_data}"
TABLES="${EVAL_DIR}/tables"
PAIRS_PER_RUN="${PAIRS_PER_RUN:-10000}"
FORCE=0
SKIP_STREX=0
SKIP_ARGOS=0
ESM_URL="https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-019-1891-0/MediaObjects/13059_2019_1891_MOESM1_ESM.xlsx"

show_help() {
    cat <<EOF
Usage: $0 [options]

Download classify reads for both datasets:
  data/strain-exclusion/   Zenodo Mason-style accuracy FASTQs
  data/fda-argos/          the 25 Table S5 FDA-ARGOS SRR accessions

Library/taxonomy is not downloaded here (see download_library.sh /
scripts/standard_installation.sh).

Options:
  -o BASE_DIR   Base directory (default: \$KRAKEN2_BENCH_ROOT or ./benchmark_data)
  -p N          Pairs per ARGOS run (default: ${PAIRS_PER_RUN})
  -s            Skip strain-exclusion
  -a            Skip FDA-ARGOS
  -F            Force re-download / re-extract
  -h            Help
EOF
    exit 0
}

while getopts ":o:p:saFh" opt; do
    case "${opt}" in
        o) BASE_DIR="${OPTARG}" ;;
        p) PAIRS_PER_RUN="${OPTARG}" ;;
        s) SKIP_STREX=1 ;;
        a) SKIP_ARGOS=1 ;;
        F) FORCE=1 ;;
        h) show_help ;;
        \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
        :)  echo "Option -${OPTARG} requires an argument." >&2; exit 1 ;;
    esac
done

DATA_DIR="${KRAKEN2_DATA_DIR:-${BASE_DIR}/data}"
SE_DIR="${DATA_DIR}/strain-exclusion"
ARGOS_DIR="${DATA_DIR}/fda-argos"
XLSX="${EVAL_DIR}/kraken2_sup.xlsx"

# shellcheck source=zenodo_ref.sh
source "${EVAL_DIR}/zenodo_ref.sh"

mkdir -p "${SE_DIR}" "${ARGOS_DIR}" "${TABLES}"

ensure_sup_tables() {
    if [[ -s "${TABLES}/table_s3.csv" && -s "${TABLES}/table_s5.csv" && "${FORCE}" -eq 0 ]]; then
        echo "[*] Tables already present in ${TABLES}"
        return 0
    fi
    if [[ ! -s "${XLSX}" || "${FORCE}" -eq 1 ]]; then
        echo "[*] Downloading Additional file 1 XLSX..."
        wget -O "${XLSX}" "${ESM_URL}"
    fi
    python3 "${EVAL_DIR}/export_sup_tables.py" --xlsx "${XLSX}" --out-dir "${TABLES}"
}

pick_taxonomy() {
    local cand
    for cand in \
        "${BASE_DIR}/db/standard/taxonomy" \
        "${BASE_DIR}/db/strain_exclusion/taxonomy" \
        "${BASE_DIR}/db/taxonomy" \
        "${BASE_DIR}/db_strain_exclusion/taxonomy"
    do
        if [[ -s "${cand}/names.dmp" && -s "${cand}/nodes.dmp" ]]; then
            echo "${cand}"
            return 0
        fi
    done
    return 1
}

echo "========================================================="
echo " download_data"
echo "========================================================="
echo "  Base:       ${BASE_DIR}"
echo "  Strain-excl: ${SE_DIR}"
echo "  FDA-ARGOS:  ${ARGOS_DIR}"
echo "========================================================="

ensure_sup_tables

if [[ "${SKIP_STREX}" -eq 0 ]]; then
    echo "[*] Strain-exclusion (Zenodo accuracy FASTQs)..."
    strex_ensure_zenodo "${FORCE}"
    if [[ -z "${STREX_FASTA:-}" && -z "${ACC_R1:-}" && -z "${GENOME_FQ_DIRS:-}" ]]; then
        echo "Error: locate_strex.py did not find FASTQs or FASTA in ${ZENODO_EXTRACT}" >&2
        exit 1
    fi
    ACC_OUT1="${SE_DIR}/R1.fastq"
    ACC_OUT2="${SE_DIR}/R2.fastq"
    TRUTH="${SE_DIR}/truth.tsv"
    if [[ -n "${ACC_R1:-}" && -n "${ACC_R2:-}" ]]; then
        echo "[*] Copying packaged accuracy FASTQs → ${SE_DIR}"
        if [[ "${FORCE}" -eq 1 || ! -f "${ACC_OUT1}" || -L "${ACC_OUT1}" || ! -f "${ACC_OUT2}" || -L "${ACC_OUT2}" ]]; then
            rm -f "${ACC_OUT1}" "${ACC_OUT2}"
            cp -f "${ACC_R1}" "${ACC_OUT1}"
            cp -f "${ACC_R2}" "${ACC_OUT2}"
        else
            echo "[*] Already present: ${ACC_OUT1}"
        fi
    elif [[ -n "${GENOME_FQ_DIRS:-}" ]]; then
        echo "[*] Building accuracy subset (first 1000 fragments / genome)..."
        rm -f "${ACC_OUT1}" "${ACC_OUT2}"
        IFS=',' read -r -a _dirs <<< "${GENOME_FQ_DIRS}"
        for d in "${_dirs[@]}"; do
            shopt -s nullglob
            files1=("${d}"/*_1.fq "${d}"/*_1.fastq)
            shopt -u nullglob
            for f1 in "${files1[@]}"; do
                [[ -f "${f1}" ]] || continue
                f2="${f1/_1.fq/_2.fq}"
                [[ -f "${f2}" ]] || f2="${f1/_1.fastq/_2.fastq}"
                if [[ ! -f "${f2}" ]]; then
                    echo "Warning: missing mate for ${f1}" >&2
                    continue
                fi
                head -qn 4000 "${f1}" >> "${ACC_OUT1}"
                head -qn 4000 "${f2}" >> "${ACC_OUT2}"
            done
        done
    else
        echo "Error: no accuracy FASTQs in the Zenodo extract." >&2
        exit 1
    fi
    python3 "${EVAL_DIR}/make_strex_truth.py" \
        --fastq "${ACC_OUT1}" \
        --table-s3 "${TABLES}/table_s3.csv" \
        --out "${TRUTH}"
else
    echo "[*] Skipping strain-exclusion (-s)"
fi

if [[ "${SKIP_ARGOS}" -eq 0 ]]; then
    echo "[*] FDA-ARGOS from ${TABLES}/table_s5.csv..."
    argos_extra=()
    [[ "${FORCE}" -eq 1 ]] && argos_extra+=(--force)
    python3 "${EVAL_DIR}/prepare_fda_argos_reads.py" \
        --table-s5 "${TABLES}/table_s5.csv" \
        --out-dir "${ARGOS_DIR}" \
        --pairs "${PAIRS_PER_RUN}" \
        "${argos_extra[@]}"
    TAX_DIR=""
    TAX_DIR="$(pick_taxonomy || true)"
    if [[ -n "${TAX_DIR}" ]]; then
        python3 "${EVAL_DIR}/make_strex_truth.py" \
            --fastq "${ARGOS_DIR}/R1.fastq" \
            --table-s5 "${TABLES}/table_s5.csv" \
            --names-dmp "${TAX_DIR}/names.dmp" \
            --nodes-dmp "${TAX_DIR}/nodes.dmp" \
            --out "${ARGOS_DIR}/truth.tsv"
    else
        echo "  Note: no names.dmp yet; skip ARGOS truth.tsv." >&2
        echo "  Run ./benchmark/kraken2/download_library.sh then re-run:" >&2
        echo "    python3 ${EVAL_DIR}/make_strex_truth.py --fastq ${ARGOS_DIR}/R1.fastq \\" >&2
        echo "      --table-s5 ${TABLES}/table_s5.csv --names-dmp ${BASE_DIR}/db/standard/taxonomy/names.dmp \\" >&2
        echo "      --nodes-dmp ${BASE_DIR}/db/standard/taxonomy/nodes.dmp --out ${ARGOS_DIR}/truth.tsv" >&2
    fi
else
    echo "[*] Skipping FDA-ARGOS (-a)"
fi

echo ""
echo "download_data complete."
echo "  Strain-excl: ${SE_DIR}"
echo "  FDA-ARGOS: ${ARGOS_DIR}"
echo "Library: ./benchmark/kraken2/download_library.sh -o ${BASE_DIR}"
