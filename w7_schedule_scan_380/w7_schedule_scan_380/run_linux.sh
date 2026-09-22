#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
QEC_PYTHON="${PYTHON_BIN:-.venv-qec312/bin/python}"
exec "$QEC_PYTHON" -u main.py "$@"
