#!/usr/bin/env bash
# Fetch/build/test/archive v2 preparation only. No acceptance histories run.
set -euo pipefail
if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
    echo "Requires an external native Linux amd64 server" >&2
    exit 2
fi
if [[ $# -lt 4 || $# -gt 5 ]]; then
    echo "Usage: $0 PREVIOUS_RESULTS PREVIOUS_SUMMARY TEST_DATA_DIR REMOTE_REF [WORKERS]" >&2
    exit 2
fi
script_dir="$(cd "$(dirname "$0")" && pwd)"
previous="$(cd "$1" && pwd)"
summary="$(python3 -c 'import pathlib,sys; print(pathlib.Path(sys.argv[1]).resolve())' "$2")"
workers="${5:-16}"
python3 "$script_dir/../server/run_tests.py" --server --ref "$4" --data-dir "$3" --jobs "$workers" \
    --input "$previous/manifest.json" --input "$summary" -- \
    python3 tests/crpropa/prepare_transport_photon_v2.py \
    --previous-results "$previous" --previous-summary "$summary" --workers "$workers" --hours 20
