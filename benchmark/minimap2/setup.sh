#!/usr/bin/env bash
# One-time conda env setup, matching kraken2's own CONDA_ENV pattern.
# `pbsim2` provides a real ONT/PacBio error-model long-read simulator to
# eventually replace lib/make_dataset.py's lightweight built-in one (see that
# file's own docstring); `samtools` is for anyone who wants to convert this
# harness's PAF output to sorted BAM for external tools.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
conda env create -f "$HERE/environment.yml" || conda env update -f "$HERE/environment.yml"
echo "Activate with: conda activate minimap2_bench"
echo "Then point PYTHON at its python3 in lib/constants.sh (or export PYTHON=...)."
