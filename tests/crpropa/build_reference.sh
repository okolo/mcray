#!/usr/bin/env bash
# Freeze the approved CRPropa sources and build a reproducible Linux amd64 adapter.
set -euo pipefail
if [[ $# != 2 ]]; then
    echo "Usage: $0 /path/to/CRPropa3 /path/to/new-build-output" >&2
    exit 2
fi
script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_dir="$(cd "$script_dir/../.." && pwd)"
reference_dir="$(cd "$1" && pwd)"
output_dir="$(python3 -c 'import pathlib,sys; print(pathlib.Path(sys.argv[1]).resolve())' "$2")"
revision=73fa865b3304113da2059755a83f5939c1039d68
[[ "$(git -C "$reference_dir" rev-parse HEAD)" == "$revision" ]] || { echo "Wrong CRPropa revision" >&2; exit 1; }
python3 - "$reference_dir" "$output_dir" "$repo_dir" "$revision" <<'PY'
import hashlib,json,pathlib,shutil,subprocess,sys,tarfile
reference,output,repo=map(pathlib.Path,sys.argv[1:4]);revision=sys.argv[4]
output.mkdir(parents=True,exist_ok=False)
subprocess.run(['git','-C',str(reference),'archive','--format=tar','-o',str(output/'reference.tar'),revision],check=True)
context=output/'context';(context/'reference').mkdir(parents=True)
with tarfile.open(output/'reference.tar') as archive: archive.extractall(context/'reference',filter='data')
shutil.copytree(repo/'tests/crpropa',context/'adapters')
shutil.copyfile(repo/'tests/nuclear_decay_probe.cpp',context/'adapters/mcray_decay_probe.cpp')
shutil.copyfile(repo/'tests/shared_photon_table.h',context/'adapters/shared_photon_table.h')
shutil.copyfile(repo/'docker/crpropa-validation.Dockerfile',context/'Dockerfile')
hashes={str(p.relative_to(context)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(context.rglob('*')) if p.is_file()}
(output/'source-manifest.json').write_text(json.dumps({'crpropa_revision':revision,'source_sha256':hashes,
 'mcray_revision':subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip(),
 'mcray_status':subprocess.check_output(['git','-C',str(repo),'status','--porcelain'],text=True)},indent=2)+'\n')
PY
docker image inspect mcray-test-linux-amd64:local > "$output_dir/base-image.json"
docker build --platform linux/amd64 --tag mcray-crpropa-validation:local \
    "$output_dir/context" > "$output_dir/build.log" 2>&1
docker image inspect mcray-crpropa-validation:local > "$output_dir/reference-image.json"
