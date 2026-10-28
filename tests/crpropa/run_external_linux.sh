#!/usr/bin/env bash
# Run long decay comparison ensembles on an external native Linux amd64 server.
set -euo pipefail
if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
    echo "This script requires an external native Linux amd64 server" >&2
    exit 2
fi
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 /path/to/prepared-data /path/to/new-results [--profile pilot|extended] [--resume]" >&2
    exit 2
fi
script_dir="$(cd "$(dirname "$0")" && pwd)"
data_dir="$1"
output_dir="$2"
shift 2
python3 "$script_dir/run_decay_comparisons.py" --data "$data_dir" --output "$output_dir" \
    --profile extended --external-server "$@"
