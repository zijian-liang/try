#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [[ -z "${QEC_PYTHON:-}" ]]; then
    if [[ -n "${PYTHON_BIN:-}" ]]; then
        QEC_PYTHON="$PYTHON_BIN"
    elif [[ -x "${VENV_DIR:-.venv}/bin/python" ]]; then
        QEC_PYTHON="${VENV_DIR:-.venv}/bin/python"
    else
        QEC_PYTHON=python3
    fi
fi
exec "$QEC_PYTHON" -u main.py "$@"
