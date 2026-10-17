#!/usr/bin/env bash
# Fetch, build, test and archive deterministic PD diagnostics, including failure.
set -euo pipefail
if [[ $# != 4 ]]; then
    echo "Usage: $0 NUCLEI_DATA PREPARED_REFERENCES TEST_DATA_DIR REMOTE_REF" >&2
    exit 2
fi
script_dir="$(cd "$(dirname "$0")" && pwd)"
data="$(cd "$1" && pwd)"
references="$(cd "$2" && pwd)"
python3 "$script_dir/../server/run_tests.py" --server --ref "$4" --data-dir "$3" \
    --input "$data/manifest.json" --input "$data/photon_manifest.json" \
    --input "$references/cmb/shared_manifest.json" --input "$references/cmb-inoue/shared_manifest.json" -- \
    python3 tests/deterministic/run_pd_matrix.py --external-server --data "$data" --references "$references"
