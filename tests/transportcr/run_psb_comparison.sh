#!/usr/bin/env bash
set -euo pipefail
reference_dir="$1"
results_dir="$2"
project_dir="$(cd "$(dirname "$0")/../.." && pwd)"
scenario="$project_dir/build/Release/tests/nuclear_scenario"
comparison="$project_dir/build/Release/tests/nuclear_comparison"
mkdir -p "$results_dir"
cd "$project_dir/bin"
"$scenario" --process psb --scan --output "$results_dir/rates.tsv"
distance="$(python3 - "$results_dir/rates.tsv" <<'PY'
import re, sys
with open(sys.argv[1]) as f:
    print(re.search(r'L_mpc=([^ ]+)', f.readline()).group(1))
PY
)"
energy=2.238721138568347e18
python3 "$project_dir/tests/transportcr/prepare_case.py" \
    --transportcr "$reference_dir" --output "$results_dir/reference" --process psb \
    --energy-per-nucleon "$energy" --distance-mpc "$distance" \
    --resolution 16 --step-mpc 0.25 --grid-span-decades 0.1
(
    cd "$results_dir/reference"
    "$reference_dir/src/Release/propagation" case.xsw > run.log 2>&1
)
"$scenario" --process psb --batches 20 --primaries 10000 --seed 20261011 \
    --energy-per-nucleon "$energy" --output "$results_dir/mc.tsv" > "$results_dir/mc.log" 2>&1
"$comparison" --mc "$results_dir/mc.tsv" --batches 20 --primaries 10000 --process psb \
    --energy-per-nucleon "$energy" \
    --initial "$results_dir/reference/results/case/step0/all" \
    --final "$results_dir/reference/results/case/fin/all" \
    --gtest_output="xml:$results_dir/comparison.xml" | tee "$results_dir/comparison.log"
