#!/usr/bin/env bash
# Independent method pilots: resolve source, build, test, run, archive. No acceptance.
set -euo pipefail
if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
    echo "Requires an external native Linux amd64 server" >&2
    exit 2
fi
if [[ $# -lt 4 || $# -gt 6 ]]; then
    echo "Usage: $0 INPUT_BUNDLE PILOT_REFERENCES TEST_DATA_DIR REMOTE_REF [WORKERS] [smoke|pilot]" >&2
    exit 2
fi
script_dir="$(cd "$(dirname "$0")" && pwd)"
bundle="$(cd "$1" && pwd)"
references="$(cd "$2" && pwd)"
workers="${5:-16}"
mkdir -p "$3"
results="$(cd "$3" && pwd)"
profile="${6:-pilot}"
python3 "$script_dir/../server/run_tests.py" --server --ref "$4" --data-dir "$results" --jobs "$workers" \
    --input "$bundle/inputs-manifest.json" --input "$references/reference-build.json" -- \
    python3 tests/crpropa/run_rare_transport_pilots.py --external-server \
    --bundle "$bundle" --references "$references" --profile "$profile" --workers "$workers" \
    --campaign-ledger "$results/rare-transport-v2-${profile}-ledger.json"
