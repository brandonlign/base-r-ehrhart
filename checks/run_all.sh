#!/bin/sh
# Runs every standard-library check and writes the output to checks/output.txt.
# codim4_symbolic.py needs sympy and is run separately.
set -eu
cd "$(dirname "$0")"
PY="${PYTHON:-python3}"
{
  for f in ehrhart.py rows.py bases_2_3.py bases_4_5.py disks.py; do
    echo "== $f"
    "$PY" "$f"
  done
} 2>&1 | tee output.txt
