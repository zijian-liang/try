#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
QEC_PYTHON="${PYTHON_BIN:-python3.12}"
QEC_ENV="${VENV_DIR:-.venv-qec312}"
"$QEC_PYTHON" -c 'import sys,platform; assert sys.version_info[:2] == (3,12), "Use Python 3.12 for the tested decoder wheel"; print(sys.version); print(platform.machine(), platform.libc_ver())'
if [[ ! -x "$QEC_ENV/bin/python" ]]; then
    "$QEC_PYTHON" -m venv "$QEC_ENV"
fi
"$QEC_ENV/bin/python" -m ensurepip --upgrade
"$QEC_ENV/bin/python" -m pip install --upgrade pip
"$QEC_ENV/bin/python" -m pip install -r requirements.txt
"$QEC_ENV/bin/python" -c 'import stim; from tesseract_decoder import tesseract; print("Ready; Stim",stim.__version__)'
printf '%s\n' 'Run: bash run_linux.sh --workers 380'
