# Shared Zenodo 3365797 paths + fetch (Wood et al. 2019 reference data).
# Source after BASE_DIR and EVAL_DIR are set.
# shellcheck shell=bash

ZENODO_URL="${ZENODO_URL:-https://zenodo.org/records/3365797/files/kraken2-reference-data.tar.gz?download=1}"
ZENODO_DIR="${BASE_DIR}/zenodo"
ZENODO_TAR="${ZENODO_DIR}/kraken2-reference-data.tar.gz"
ZENODO_EXTRACT="${ZENODO_DIR}/extract"
ZENODO_LOCATE="${ZENODO_DIR}/locate.env"

strex_ensure_zenodo() {
    local force="${1:-0}"
    mkdir -p "${ZENODO_DIR}"
    if [[ ! -s "${ZENODO_TAR}" || "${force}" -eq 1 ]]; then
        echo "[*] Downloading Zenodo reference tarball (~9 GB)..."
        wget -O "${ZENODO_TAR}.partial" "${ZENODO_URL}"
        mv "${ZENODO_TAR}.partial" "${ZENODO_TAR}"
    else
        echo "[*] Zenodo tarball present: ${ZENODO_TAR}"
    fi
    if [[ ! -d "${ZENODO_EXTRACT}" || "${force}" -eq 1 ]]; then
        echo "[*] Extracting ${ZENODO_TAR} → ${ZENODO_EXTRACT}..."
        rm -rf "${ZENODO_EXTRACT}"
        mkdir -p "${ZENODO_EXTRACT}"
        tar -C "${ZENODO_EXTRACT}" -xzf "${ZENODO_TAR}"
    fi
    python3 "${EVAL_DIR}/locate_strex.py" --root "${ZENODO_EXTRACT}" --out "${ZENODO_LOCATE}"
    # shellcheck disable=SC1090
    source "${ZENODO_LOCATE}"
}
