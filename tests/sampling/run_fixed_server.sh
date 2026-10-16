#!/usr/bin/env bash
# Run the frozen, unpublished correction from a companion Git bundle.
set -euo pipefail
sampling_dir="$(cd "$(dirname "$0")" && pwd)"
package_dir="$(cd "$sampling_dir/.." && pwd)"
cd "$package_dir"
sha256sum --check SHA256SUMS
python_command="${PYTHON:-python3}"
requirements="$sampling_dir/secondary-requirements.txt"
if "$python_command" -c 'import sys; sys.exit(sys.version_info >= (3,12))'; then
    requirements="$sampling_dir/secondary-requirements-py310.txt"
fi
"$python_command" -c 'import sys; assert sys.version_info >= (3,10), "Python 3.10 or newer is required"'
"$python_command" -m venv "$package_dir/.venv"
"$package_dir/.venv/bin/pip" install -r "$requirements"
branch="$("$package_dir/.venv/bin/python" -c 'import json; print(json.load(open("package-manifest.json"))["branch"])')"
exec "$package_dir/.venv/bin/python" "$sampling_dir/run_secondary_server.py" \
    --server --repo "$package_dir/mcray.bundle" --ref "$branch" "$@"
