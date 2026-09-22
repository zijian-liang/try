#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
QEC_ENV="${VENV_DIR:-.venv}"
if [[ -z "${QEC_PYTHON:-}" ]]; then
    if [[ -n "${PYTHON_BIN:-}" ]]; then
        QEC_PYTHON="$PYTHON_BIN"
    elif [[ -x "$QEC_ENV/bin/python" ]]; then
        QEC_PYTHON="$QEC_ENV/bin/python"
    else
        QEC_PYTHON=python3
    fi
fi
if [[ ! -x "$QEC_ENV/bin/python" ]]; then
    "$QEC_PYTHON" -m venv "$QEC_ENV"
fi
"$QEC_ENV/bin/python" -m pip install --upgrade pip
"$QEC_ENV/bin/python" -m pip install -r requirements.txt
"$QEC_ENV/bin/python" -c 'import stim; from tesseract_decoder import tesseract; print("Dependencies ready. Stim", stim.__version__)'
printf '%s\n' 'Run: bash run_linux.sh --workers 8'
