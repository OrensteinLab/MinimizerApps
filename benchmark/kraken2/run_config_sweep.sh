#!/usr/bin/env bash
# Run a YAML-defined multi-config experiment sweep (see experiments/*.yaml).
# Mirrors benchmark/kmc/run_config_sweep.sh's pattern.
#
# Usage:
#   ./benchmark/kraken2/run_config_sweep.sh list
#   ./benchmark/kraken2/run_config_sweep.sh sweep_schemes

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODEBASE="$(cd "${SCRIPT_DIR}/../../kraken2" && pwd)"
# shellcheck source=lib/experiment.sh
source "${SCRIPT_DIR}/lib/experiment.sh"

if [[ "${1:-}" == "list" || $# -eq 0 ]]; then
    list_experiments
    exit 0
fi

EXPERIMENT="$1"
run_experiment
