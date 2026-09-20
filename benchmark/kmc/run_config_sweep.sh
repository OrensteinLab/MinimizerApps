#!/usr/bin/env bash
# Wrapper for experiments/sweep_config.yaml (see benchmark/README.md).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "$HERE/kmc3_bench.sh" experiment sweep_config "$@"
