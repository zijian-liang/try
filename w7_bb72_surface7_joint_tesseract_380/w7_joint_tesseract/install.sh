#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
PYTHON_BIN="${PYTHON_BIN:-python3}"
"$PYTHON_BIN" -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -c 'import stim; from tesseract_decoder import tesseract; print("Dependencies ready. Stim", stim.__version__)'
printf '%s\n' 'Run: .venv/bin/python main.py --workers 380'
