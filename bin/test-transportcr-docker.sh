#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "$0")/.." && pwd)"
reference_dir="${1:-$project_dir/../transportcr}"
if [[ $# -gt 0 ]]; then shift; fi
reference_revision="${TRANSPORTCR_REVISION:-94eca0c886df80b1d0b1dd09dae95b0483c2273c}"
context_dir="$(mktemp -d "${TMPDIR:-/tmp}/mcray-reference.XXXXXX")"
trap 'rm -rf "$context_dir"' EXIT
# Copy only the pinned committed reference, never user results or local changes.
git clone --no-local --no-checkout "$reference_dir" "$context_dir/transportcr"
git -C "$context_dir/transportcr" checkout --detach "$reference_revision"
cp "$project_dir/docker/transportcr-validation.Dockerfile" "$context_dir/Dockerfile"
docker build --platform linux/amd64 --file "$project_dir/docker/linux-amd64.Dockerfile" \
    --tag mcray-test-linux-amd64:local "$project_dir"
docker build --platform linux/amd64 --tag mcray-transportcr-linux-amd64:local "$context_dir"
mkdir -p "$project_dir/build"
results_dir="$(mktemp -d "$project_dir/build/transportcr-validation.XXXXXX")"
printf '%s\n' "$reference_revision" > "$results_dir/reference-revision.txt"
docker run --rm --platform linux/amd64 --mount "type=bind,source=$results_dir,target=/results" \
    mcray-transportcr-linux-amd64:local "$@"
echo "Validation outputs: $results_dir"
