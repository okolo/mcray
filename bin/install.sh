#!/usr/bin/env bash
# Build CRbeam and install the pinned full Fe56 nuclear table package.
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash bin/install.sh [--install-deps] [--jobs N]

Build Release, run CTest, and prepare/verify bin/tables/nuclear.
Supported platform: Linux amd64; --install-deps requires Ubuntu 22.04/24.04.
--install-deps  Install Ubuntu build/Python packages using apt-get and sudo.
--jobs N        Build/test parallelism (default: MCRAY_BUILD_JOBS or 2).
--help          Show this help without changing files.

Generated sources, Python environment, temporary files and installation logs
stay in .crbeam-install/ inside this checkout. Existing nuclear data are
verified and preserved; incomplete or incompatible packages are never replaced.
EOF
}

fail() { echo "error: $*" >&2; exit 1; }
jobs="${MCRAY_BUILD_JOBS:-2}"
install_deps=false
while [[ $# -gt 0 ]]; do
    case "$1" in
        --install-deps) install_deps=true; shift ;;
        --jobs)
            [[ $# -ge 2 ]] || fail "--jobs requires a positive integer"
            jobs="$2"; shift 2 ;;
        --help|-h) usage; exit 0 ;;
        *) fail "Unknown argument: $1 (use --help)" ;;
    esac
done
[[ "$jobs" =~ ^[1-9][0-9]*$ ]] || fail "--jobs must be a positive integer"
[[ "$(uname -s)" == Linux && "$(uname -m)" == x86_64 ]] || fail "Linux amd64 is required"

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$root"
if $install_deps; then
    [[ -f /etc/os-release ]] || fail "Cannot identify Ubuntu"
    # Read the operating system's identity, never a downloaded shell fragment.
    . /etc/os-release
    [[ "$ID" == ubuntu && ( "$VERSION_ID" == 22.04 || "$VERSION_ID" == 24.04 ) ]] ||
        fail "--install-deps supports Ubuntu 22.04/24.04; install dependencies manually otherwise"
    privilege=()
    if [[ "$(id -u)" != 0 ]]; then
        command -v sudo >/dev/null || fail "sudo is required for --install-deps"
        privilege=(sudo)
    fi
    "${privilege[@]}" apt-get update
    "${privilege[@]}" apt-get install -y ca-certificates git build-essential cmake \
        gfortran libgsl-dev libboost-filesystem-dev libboost-system-dev \
        libyaml-cpp-dev libssl-dev python3 python3-venv python3-pip util-linux
fi
for tool in git cmake ctest make gcc g++ gfortran python3 flock; do
    command -v "$tool" >/dev/null || fail "Missing $tool; use --install-deps or install dependencies manually"
done

state="$root/.crbeam-install"
tables="$root/bin/tables/nuclear"
revision=34b10e27375c0e9894d7c984ac15fb96e12975ef
mkdir -p "$state/tmp"
exec 9>"$state/install.lock"
flock -n 9 || fail "Another installer is using this checkout"
run="$(mktemp -d "$state/run.XXXXXXXX")"
export TMPDIR="$state/tmp"
export PIP_CACHE_DIR="$state/pip-cache"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "Installation log/work: $run"
exec 3>&1
exec >"$run/install.log" 2>&1
status() { echo "$*"; echo "$*" >&3; }
trap 'code=$?; if (( code != 0 )); then status "Installation failed ($code). Retained log/work: $run"; fi' EXIT
application_revision="$(git rev-parse HEAD)"
echo "Application source: $application_revision"
echo "Nuclear source: $revision"

verify_package() {
    local directory="$1"
    python3 "$root/tools/prepare_nuclear_photon_data.py" --verify --output "$directory"
    # Hash verification alone also accepts smaller packages. Require the full
    # GZK/decay/photon extension and current preparers before reporting success.
    python3 - "$directory" "$root" <<'PY'
import hashlib
import json
from pathlib import Path
import sys

directory, root = map(Path, sys.argv[1:])
for name in ('fallback_decays.tsv', 'gzk_inventory.json', 'decays.tsv',
             'decay_stable.tsv', 'pd_photon_exceptions.tsv', 'pd_photons.tsv',
             'pd_photon_transitions.tsv', 'decay_photon_channels.tsv', 'decay_photons.tsv'):
    if not (directory / name).is_file():
        raise SystemExit('Incomplete full-process nuclear package: ' + name)
for manifest, tool in (('manifest.json', 'prepare_gzk_data.py'),
                       ('photon_manifest.json', 'prepare_nuclear_photon_data.py')):
    data = json.loads((directory / manifest).read_text())
    if not data.get('payload_sha256'):
        raise SystemExit('Missing payload hash inventory: ' + manifest)
    expected = hashlib.sha256((root / 'tools' / tool).read_bytes()).hexdigest()
    if data.get('preparation_sha256') != expected:
        raise SystemExit('Existing package uses a different preparer: ' + manifest +
                         '; preserve it and select a fresh checkout for installation')
PY
}

existing=false
if [[ -e "$tables" || -L "$tables" ]]; then
    status "Checking existing nuclear package; it will not be replaced"
    verify_package "$tables"
    existing=true
fi

status "Building CRbeam Release ($jobs parallel jobs)"
MCRAY_BUILD_JOBS="$jobs" bash "$root/bin/build.sh" Release
status "Running CTest"
ctest --test-dir "$root/build/Release" --output-on-failure --parallel "$jobs"

if ! $existing; then
    # Every attempt gets new directories, so an interrupted preparation is
    # retained for diagnosis but is never reused as a valid installed package.
    status "Fetching the pinned nuclear reference and preparing tables"
    repo="$run/CRPropa3-data"
    git clone --filter=blob:none --sparse https://github.com/CRPropa/CRPropa3-data.git "$repo"
    git -C "$repo" fetch origin "$revision"
    git -C "$repo" checkout --detach "$revision"
    [[ "$(git -C "$repo" rev-parse HEAD)" == "$revision" ]] || fail "Wrong nuclear source revision"

    python3 -m venv "$run/venv"
    python="$run/venv/bin/python"
    "$python" -m pip install --disable-pip-version-check \
        'numpy==2.2.6' 'scipy==1.15.3' 'GitPython==3.1.45'
    "$python" -m pip freeze > "$run/python-packages.txt"
    "$python" "$root/tools/prepare_talys_data.py" --repo "$repo" \
        --output "$run/pd" --grid-policy crpropa-padding
    "$python" "$root/tools/prepare_nuclear_decay_data.py" --repo "$repo" \
        --talys-data "$run/pd" --output "$run/decay"
    "$python" "$root/tools/prepare_gzk_data.py" --repo "$repo" \
        --decay-data "$run/decay" --output "$run/gzk"
    "$python" "$root/tools/prepare_nuclear_photon_data.py" --repo "$repo" \
        --data "$run/gzk" --decay-reference "$run/decay" --output "$run/nuclear"
    verify_package "$run/nuclear"
    # GNU mv's no-clobber operation preserves even a concurrently created
    # destination. Refuse success if the staged package was not moved.
    mv -T --no-clobber "$run/nuclear" "$tables"
    [[ ! -e "$run/nuclear" ]] || fail "Nuclear destination appeared during preparation; preserved both packages"
fi
verify_package "$tables"
status "Installation completed. CRbeam: $root/build/Release/CRbeam"
status "Nuclear tables: $tables"
status "Run from $root/bin, for example: ../build/Release/CRbeam /absolute/path/to/simulation.yaml"
status "See examples/config.full.yaml for opt-in nuclear process settings."
