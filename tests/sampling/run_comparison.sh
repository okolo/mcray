#!/usr/bin/env bash
# Workload for tests/server/run_tests.py on an external Linux amd64 server.
set -euo pipefail
: "${MCRAY_RESULTS_DIR:?Use the external server runner}"
: "${MCRAY_NR_REFERENCE_DIR:?Provide licensed NR headers on the server}"
[[ $(uname -s) == Linux && $(uname -m) == x86_64 ]] || { echo 'Requires external Linux x86_64' >&2; exit 2; }
cmake -S src/app/crbeam -B build/sampling -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_STANDARD=17 -DMCRAY_NR_REFERENCE_DIR="$MCRAY_NR_REFERENCE_DIR"
cmake --build build/sampling --target em_sampling_probe -j "${MCRAY_JOBS:-4}"
cp build/sampling/tests/sampling_reference/manifest.json "$MCRAY_RESULTS_DIR/sampling-sources.json"
probe="$PWD/build/sampling/tests/em_sampling_probe"
"$probe" kernels 200000 16384 1e-3 > "$MCRAY_RESULTS_DIR/kernels.tsv"
# Per-primary moments allow valid errors for correlated cascade descendants.
# These defaults are exploratory; no physical acceptance criterion is implied.
for seed in 20261012 20261013 20261014; do
    for sampler in old current; do
        run_seed="$seed"
        if [[ "$sampler" == current ]]; then run_seed=$((seed + 100)); fi
        "$probe" cascade "$sampler" "${MCRAY_CASCADE_N:-10000}" "$run_seed" \
            "$MCRAY_RESULTS_DIR/${sampler}-${seed}.tsv" 0.1 "${MCRAY_RATE_VARIATION:-0}"
    done
    python3 tests/sampling/summarize.py "$MCRAY_RESULTS_DIR/old-${seed}.tsv" \
        "$MCRAY_RESULTS_DIR/current-${seed}.tsv" > "$MCRAY_RESULTS_DIR/comparison-${seed}.tsv"
done
