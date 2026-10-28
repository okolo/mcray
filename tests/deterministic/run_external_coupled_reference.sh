#!/usr/bin/env bash
# Fetch/build/test/archive finite-rate PD+decay deterministic refinements.
set -euo pipefail
if [[ $# -lt 4 || $# -gt 5 ]]; then
    echo "Usage: $0 NUCLEI_DATA PREPARED_REFERENCES TEST_DATA_DIR REMOTE_REF [initial|full]" >&2
    exit 2
fi
script_dir="$(cd "$(dirname "$0")" && pwd)"
data="$(cd "$1" && pwd)"
references="$(cd "$2" && pwd)"
python3 "$script_dir/../server/run_tests.py" --server --ref "$4" --data-dir "$3" \
    --input "$data/manifest.json" --input "$data/decay_manifest.json" --input "$data/photon_manifest.json" \
    --input "$references/cmb/shared_manifest.json" --input "$references/cmb-inoue/shared_manifest.json" -- \
    python3 tests/deterministic/run_coupled_matrix.py --external-server --data "$data" --references "$references" --profile "${5:-initial}"
