# Run the lab-authored benchmark harness's test suite (see tests/).
#
# `make venv` creates ../../kmc/.venv (same venv lib/constants.sh's bench_py()
# auto-detects and prefers) and installs pytest into it. `make test` then
# picks it up automatically; override with PYTHON= to use something else.

VENV := $(abspath ../../kmc/.venv)
VENV_PYTHON ?= python3.11

ifneq ($(wildcard $(VENV)/bin/python),)
PYTHON ?= $(VENV)/bin/python
else
PYTHON ?= python3
endif

.PHONY: test venv

venv:
	$(VENV_PYTHON) -m venv $(VENV)
	$(VENV)/bin/pip install --upgrade pip
	$(VENV)/bin/pip install pytest

test:
	$(PYTHON) -m pytest tests/
