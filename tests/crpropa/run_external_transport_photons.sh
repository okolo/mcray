#!/usr/bin/env bash
# Fetch/build/test/archive pre-acceptance nuclear transport and photon pilots.
set -euo pipefail
if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
    echo "Requires an external native Linux amd64 server" >&2
    exit 2
fi
if [[ $# -lt 3 || $# -gt 5 ]]; then
    echo "Usage: $0 INPUT_BUNDLE TEST_DATA_DIR REMOTE_REF [smoke|pilot] [initial|full]" >&2
    exit 2
fi
script_dir="$(cd "$(dirname "$0")" && pwd)"
bundle="$(cd "$1" && pwd)"
# The workload path is relative to the fresh checkout, never to this launcher.
python3 "$script_dir/../server/run_tests.py" --server --ref "$3" --data-dir "$2" \
    --input "$bundle/inputs-manifest.json" -- \
    python3 tests/crpropa/run_transport_photon_pilots.py --external-server \
    --bundle "$bundle" --profile "${4:-pilot}" --matrix "${5:-initial}"
