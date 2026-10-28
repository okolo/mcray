#!/usr/bin/env bash
# Fetch/build/test/archive the frozen initial-point acceptance on Linux amd64.
set -euo pipefail
if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
    echo "Requires an external native Linux amd64 server" >&2
    exit 2
fi
if [[ $# -lt 5 || $# -gt 6 ]]; then
    echo "Usage: $0 INPUT_BUNDLE PILOT_REFERENCES DETERMINISTIC_CACHE TEST_DATA_DIR REMOTE_REF [WORKERS]" >&2
    exit 2
fi
script_dir="$(cd "$(dirname "$0")" && pwd)"
bundle="$(cd "$1" && pwd)"
references="$(cd "$2" && pwd)"
cache="$(cd "$3" && pwd)"
workers="${6:-8}"
mkdir -p "$4"
results_root="$(cd "$4" && pwd)"
protocol_sha="$(python3 -c 'import hashlib,pathlib,sys; print(hashlib.sha256(pathlib.Path(sys.argv[1]).read_bytes()).hexdigest())' "$script_dir/protocols/transport_photons_v1.json")"
python3 "$script_dir/../server/run_tests.py" --server --ref "$5" --data-dir "$4" --jobs "$workers" \
    --input "$bundle/inputs-manifest.json" --input "$references/reference-build.json" \
    --input "$script_dir/protocols/transport_photons_v1.json" -- \
    python3 tests/crpropa/run_transport_photon_acceptance.py --external-server \
    --bundle "$bundle" --references "$references" --deterministic-cache "$cache" --workers "$workers" \
    --campaign-ledger "$results_root/transport-photon-acceptance-v1-ledger.json" \
    --expected-protocol-sha "$protocol_sha"
