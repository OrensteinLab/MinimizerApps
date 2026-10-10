# Shared PATH / conda setup for benchmark scripts. Source from other scripts:
#   source "$(dirname "${BASH_SOURCE[0]}")/benchmark_setup.sh"
#
# This file lives in benchmark/kraken2/ (not upstream Kraken 2 scripts/),
# outside the kraken2/ source tree so upstream diffs stay harness-free.
# Kraken 2 wrappers stay in repo scripts/ (k2, kraken2-build, …).

_eval_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_bench_codebase="$(cd "${_eval_dir}/../../kraken2" && pwd)"
_k2_scripts="${_bench_codebase}/scripts"

if [[ -f "${_eval_dir}/.env" ]]; then
    # shellcheck source=/dev/null
    source "${_eval_dir}/.env"
fi

_bench_conda_bin() {
    local env="$1"
    [[ -z "${env}" ]] && return 0
    if [[ -d "${env}/bin" ]]; then
        echo "${env}/bin"
        return 0
    fi
    local home="${HOME:-}"
    if [[ -n "${home}" ]]; then
        local root
        for root in "${home}/.conda" "${home}/miniconda3" "${home}/anaconda3" "${home}/mambaforge" "${home}/miniforge3"; do
            if [[ -d "${root}/envs/${env}/bin" ]]; then
                echo "${root}/envs/${env}/bin"
                return 0
            fi
        done
    fi
    if command -v conda >/dev/null 2>&1; then
        local base
        base="$(conda info --base 2>/dev/null || true)"
        if [[ -n "${base}" && -d "${base}/envs/${env}/bin" ]]; then
            echo "${base}/envs/${env}/bin"
        fi
    fi
}

_bench_prepend_path() {
    local dir="$1"
    [[ -z "${dir}" || ! -d "${dir}" ]] && return 0
    local rest=":${PATH}:"
    rest="${rest//:${dir}:/:}"
    rest="${rest#:}"
    rest="${rest%:}"
    export PATH="${dir}:${rest}"
}

# Repo wrappers + binaries must beat conda-installed kraken2 on PATH.
# Prepend in reverse priority order (last prepend wins the front).
if [[ -n "${CONDA_PREFIX:-}" ]]; then
    _bench_prepend_path "${CONDA_PREFIX}/bin"
elif [[ -n "${CONDA_ENV:-}" ]]; then
    _conda_bin="$(_bench_conda_bin "${CONDA_ENV}")"
    _bench_prepend_path "${_conda_bin}"
fi
_bench_prepend_path "${_bench_codebase}/src"
_bench_prepend_path "${_k2_scripts}"

if ! command -v python3 >/dev/null 2>&1; then
    echo "Error: python3 not found after PATH setup." >&2
    echo "  Set CONDA_ENV in benchmark/kraken2/.env (e.g. kraken_env)." >&2
    exit 1
fi
# Prefer conda python when CONDA_ENV / CONDA_PREFIX is set.
if [[ -n "${CONDA_PREFIX:-}" ]]; then
    _py="${CONDA_PREFIX}/bin/python3"
elif [[ -n "${CONDA_ENV:-}" ]]; then
    _conda_bin="$(_bench_conda_bin "${CONDA_ENV}")"
    _py="${_conda_bin}/python3"
else
    _py="$(command -v python3)"
fi
if [[ ! -x "${_py}" ]]; then
    echo "Error: python3 not found (${_py:-unset})." >&2
    echo "  Set CONDA_ENV in benchmark/kraken2/.env or conda activate your env." >&2
    exit 1
fi

# bioconda kraken2 in the conda env (perl, wget, python runtime).
_conda_root=""
if [[ -n "${CONDA_PREFIX:-}" ]]; then
    _conda_root="${CONDA_PREFIX}"
elif [[ -n "${CONDA_ENV:-}" ]]; then
    _conda_bin="$(_bench_conda_bin "${CONDA_ENV}")"
    _conda_root="$(cd "${_conda_bin}/.." && pwd)"
fi
if [[ -n "${_conda_root}" ]]; then
    if ! compgen -G "${_conda_root}/conda-meta/kraken2-"*.json >/dev/null 2>&1; then
        echo "Error: bioconda kraken2 is not installed in conda env '${CONDA_ENV:-${CONDA_PREFIX}}'." >&2
        echo "  conda install -c bioconda -c conda-forge kraken2" >&2
        echo "  or: ./benchmark/kraken2/setup.sh" >&2
        echo "  See: https://anaconda.org/channels/bioconda/packages/kraken2/overview" >&2
        exit 1
    fi
fi

export PYTHON="${_py}"
export K2="${_k2_scripts}/k2"

unset _eval_dir _bench_codebase _k2_scripts _bench_conda_bin _bench_prepend_path _conda_bin _py _conda_root
