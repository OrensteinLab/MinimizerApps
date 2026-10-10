#!/usr/bin/env bash
# setup.sh — install host build deps + Miniconda if needed, then
# create/update the benchmark conda env with bioconda kraken2.
#
# Usage: ./benchmark/kraken2/setup.sh [ENV_NAME]
#   ENV_NAME defaults to kraken_env (matches environment.yml).
#
# Installs (when missing):
#   - apt: build-essential (g++/make/OpenMP), zlib1g-dev, time (GNU /usr/bin/time)
#   - Miniconda under $HOME/miniconda3

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODEBASE="$(cd "${SCRIPT_DIR}/../../kraken2" && pwd)"
ENV_NAME="${1:-kraken_env}"
MINICONDA_PREFIX="${MINICONDA_PREFIX:-${HOME}/miniconda3}"

have_cmd() { command -v "$1" >/dev/null 2>&1; }

run_sudo() {
    if [[ "$(id -u)" -eq 0 ]]; then
        "$@"
    elif have_cmd sudo; then
        sudo "$@"
    else
        echo "Error: need root or sudo to install: $*" >&2
        exit 1
    fi
}

install_apt_deps() {
    local -a missing=()
    # Headers needed to compile k2mask and other zlib-linked binaries.
    [[ -f /usr/include/zlib.h ]] || missing+=(zlib1g-dev)
    # g++ with OpenMP (benchmark/kraken2/README.md); build-essential pulls g++/make.
    if ! have_cmd g++; then
        missing+=(build-essential)
    else
        _omp_tmp="$(mktemp -d)"
        if ! echo 'int main(){return 0;}' \
            | g++ -x c++ -fopenmp - -o "${_omp_tmp}/t" 2>/dev/null; then
            missing+=(build-essential)
        fi
        rm -rf "${_omp_tmp}"
        unset _omp_tmp
    fi
    have_cmd make || missing+=(make)
    # GNU time (/usr/bin/time -v); bash builtin time has no -v.
    [[ -x /usr/bin/time ]] || missing+=(time)

    if [[ ${#missing[@]} -eq 0 ]]; then
        echo "[*] Host build deps already present (g++, zlib, make, time)"
        return 0
    fi

    if ! have_cmd apt-get; then
        echo "Error: missing packages (${missing[*]}) and apt-get not available." >&2
        echo "  Install a C++ toolchain with OpenMP and zlib development headers, then re-run." >&2
        exit 1
    fi

    echo "[*] Installing apt packages: ${missing[*]}"
    run_sudo apt-get update -qq
    DEBIAN_FRONTEND=noninteractive run_sudo apt-get install -y "${missing[@]}"
}

install_miniconda() {
    local arch installer url tmp
    arch="$(uname -m)"
    case "${arch}" in
        x86_64)  installer="Miniconda3-latest-Linux-x86_64.sh" ;;
        aarch64) installer="Miniconda3-latest-Linux-aarch64.sh" ;;
        *)
            echo "Error: unsupported architecture '${arch}' for Miniconda auto-install." >&2
            exit 1
            ;;
    esac
    url="https://repo.anaconda.com/miniconda/${installer}"
    tmp="$(mktemp -d)"
    trap 'rm -rf "${tmp}"' RETURN

    echo "[*] Installing Miniconda to ${MINICONDA_PREFIX}"
    if have_cmd curl; then
        curl -fsSL "${url}" -o "${tmp}/${installer}"
    elif have_cmd wget; then
        wget -q -O "${tmp}/${installer}" "${url}"
    else
        echo "Error: need curl or wget to download Miniconda." >&2
        exit 1
    fi
    bash "${tmp}/${installer}" -b -p "${MINICONDA_PREFIX}"
    "${MINICONDA_PREFIX}/bin/conda" init bash >/dev/null
}

ensure_conda() {
    # Already on PATH (e.g. activated shell).
    if have_cmd conda; then
        return 0
    fi

    # Common install locations / fresh MINICONDA_PREFIX.
    local candidate
    for candidate in \
        "${MINICONDA_PREFIX}" \
        "${HOME}/miniconda3" \
        "${HOME}/anaconda3" \
        "${HOME}/mambaforge" \
        "${HOME}/miniforge3" \
        "${HOME}/.conda"
    do
        if [[ -x "${candidate}/bin/conda" ]]; then
            # shellcheck disable=SC1091
            source "${candidate}/etc/profile.d/conda.sh"
            return 0
        fi
        if [[ -f "${candidate}/etc/profile.d/conda.sh" ]]; then
            # shellcheck disable=SC1091
            source "${candidate}/etc/profile.d/conda.sh"
            have_cmd conda && return 0
        fi
    done

    install_miniconda
    # shellcheck disable=SC1091
    source "${MINICONDA_PREFIX}/etc/profile.d/conda.sh"
    if ! have_cmd conda; then
        echo "Error: conda still not available after Miniconda install." >&2
        exit 1
    fi
}

install_apt_deps
ensure_conda

echo "[*] Creating/updating conda env '${ENV_NAME}' from ${CODEBASE}/environment.yml"
conda env create -f "${CODEBASE}/environment.yml" -n "${ENV_NAME}" 2>/dev/null \
    || conda env update -f "${CODEBASE}/environment.yml" -n "${ENV_NAME}" --prune

echo ""
echo "Done. Next steps:"
echo "  source ~/.bashrc   # if this was a fresh Miniconda install"
echo "  conda activate ${ENV_NAME}"
echo "  cp benchmark/kraken2/.env.example benchmark/kraken2/.env"
echo "  # edit benchmark.env — set KRAKEN2_BENCH_ROOT, CONDA_ENV"
echo "  make -j\"\$(nproc)\" all"
