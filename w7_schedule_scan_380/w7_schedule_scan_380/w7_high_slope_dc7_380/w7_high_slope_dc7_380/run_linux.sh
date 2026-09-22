#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
QEC_PYTHON="${PYTHON_BIN:-python3}"
exec "$QEC_PYTHON" -u main.py "$@"
