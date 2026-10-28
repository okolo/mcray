#!/usr/bin/env bash
# Fetch, build, test and archive long network rate diagnostics on external Linux.
set -euo pipefail
if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
    echo "Requires an external native Linux amd64 server" >&2
    exit 2
fi
if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: $0 INPUT_BUNDLE TEST_DATA_DIR [REMOTE_REF]" >&2
    exit 2
fi
script_dir="$(cd "$(dirname "$0")" && pwd)"
bundle="$(cd "$1" && pwd)"
python3 - "$bundle" <<'PY'
import hashlib,json,pathlib,sys
root=pathlib.Path(sys.argv[1]);manifest=json.loads((root/'inputs-manifest.json').read_text())
for name,expected in manifest['payload_sha256'].items():
    path=root/name
    if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:raise SystemExit('Checksum mismatch: '+name)
print('All transferred input/tool checksums match')
PY
python3 "$script_dir/../server/run_tests.py" --server --ref "${3:-talys}" --data-dir "$2" \
    --input "$bundle/inputs-manifest.json" -- \
    python3 "$script_dir/run_network_rate_scan.py" --external-server --profile extended \
    --data "$bundle/nuclei" --cmb-reference "$bundle/cmb" --inoue-reference "$bundle/cmb-inoue"
