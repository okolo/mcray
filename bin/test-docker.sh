#!/usr/bin/env bash
set -euo pipefail
script_dir="$(cd "$(dirname "$0")" && pwd)"
project_dir="$(cd "$script_dir/.." && pwd)"
image="mcray-test-linux-amd64:local"
docker build --platform linux/amd64 \
    --file "$project_dir/docker/linux-amd64.Dockerfile" \
    --tag "$image" "$project_dir"
docker run --rm --platform linux/amd64 "$image" "$@"
