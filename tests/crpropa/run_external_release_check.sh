#!/usr/bin/env bash
# One-shot eight-hour ordinary release check; retain partial evidence on timeout.
set -euo pipefail
if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
    echo "Requires an external native Linux amd64 server" >&2
    exit 2
fi
if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: $0 BASE_DIR PUBLISHED_COMMIT [WORKERS]" >&2
    exit 2
fi
base="$(cd "$1" && pwd)"
revision="$2"
workers="${3:-16}"
script_dir="$(cd "$(dirname "$0")" && pwd)"
results="$base/mcray-test-data"
mkdir -p "$results"
exec 9>"$results/release-check-8h.lock"
if ! flock -n 9; then
    echo "Another release-check launcher is active; no duplicate run started." >&2
    exit 2
fi
if [[ -e "$results/release-check-8h-ledger.json" ]]; then
    echo "Release check already reserved. Preserve its ledger and results; do not rerun." >&2
    exit 2
fi
export MCRAY_RELEASE_DEADLINE
MCRAY_RELEASE_DEADLINE="$(python3 -c 'import time; print(time.time()+8*3600)')"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

# Reuse retained server results. Restore originals only when their payload
# directory is absent; the workload verifies the frozen hashes before use.
restore_results() {
    local run="$1" marker="$2" retained
    retained="$results/work/$run/results"
    if [[ -f "$retained/$marker" ]]; then
        printf '%s\n' "$retained"
        return
    fi
    retained="$base/release-check-inputs/$run"
    if [[ ! -f "$retained/$marker" ]]; then
        test -f "$results/$run.tar.gz"
        mkdir -p "$base/release-check-inputs"
        tar -xzf "$results/$run.tar.gz" -C "$base/release-check-inputs"
    fi
    test -f "$retained/$marker"
    printf '%s\n' "$retained"
}
bundle="$base/rate-inputs"
references="$base/acceptance-references/20261016T102021Z-03eb6c57/transport-photons"
if [[ ! -f "$references/reference-build.json" ]]; then
    pilot="$(restore_results 20261016T102021Z-03eb6c57 transport-photons/reference-build.json)"
    references="$pilot/transport-photons"
fi
coupled="$(restore_results 20261019T205338Z-943099d8 deterministic-coupled/matrix.json)"
previous="$(restore_results 20261020T153324Z-d619dbcb manifest.json)"
summary="$results/20261020T153324Z-d619dbcb.json"
test -f "$summary"
remaining="$(python3 -c 'import os,time; print(max(1, int(float(os.environ["MCRAY_RELEASE_DEADLINE"])-time.time())))')"
python3 "$script_dir/../server/run_tests.py" --server --ref "$revision" \
    --data-dir "$results" --jobs "$workers" --max-runtime-seconds "$remaining" \
    --input "$bundle/inputs-manifest.json" --input "$references/reference-build.json" \
    --input "$previous/manifest.json" --input "$summary" -- \
    python3 tests/crpropa/run_release_check.py --external-server \
    --bundle "$bundle" --references "$references" \
    --deterministic-cache "$coupled/deterministic-coupled" \
    --previous-results "$previous" --previous-summary "$summary" \
    --workers "$workers" --campaign-ledger "$results/release-check-8h-ledger.json"
