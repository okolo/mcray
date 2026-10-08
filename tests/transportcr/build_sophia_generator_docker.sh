#!/usr/bin/env bash
# Freeze and build the independent TransportCR SOPHIA generator in a new directory.
set -euo pipefail
if [[ $# != 2 ]]; then
    echo "Usage: $0 /path/to/propagation/sources /path/to/new-output" >&2
    exit 2
fi
source_dir="$1"
output_dir="$(python3 -c 'from pathlib import Path; import sys; print(Path(sys.argv[1]).expanduser().resolve())' "$2")"
image="mcray-test-linux-amd64:local"
python3 - "$source_dir" "$output_dir" <<'PY'
import hashlib
import json
from pathlib import Path
import shutil
import sys
source = Path(sys.argv[1]).expanduser().resolve()
output = Path(sys.argv[2]).expanduser().resolve()
files = ['TParticle.h', 'sophia/Sophia.cpp', 'sophia/Sophia.h',
         'sophia/GitVersion.cpp', 'sophia/sophia2cpp.f']
files += ['sophia/SOPHIA/' + name + '.f' for name in
          ('eventgen', 'sampling', 'inpoutput', 'jetset74dp', 'rndm')]
for name in files:
    if not (source / name).is_file():
        raise SystemExit('Missing generator source: ' + str(source / name))
output.mkdir(parents=True, exist_ok=False)
hashes = {}
for name in files:
    target = output / 'source' / name
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source / name, target)
    hashes[name] = hashlib.sha256(target.read_bytes()).hexdigest()
(output / 'build').mkdir()
(output / 'source-manifest.json').write_text(json.dumps(
    {'source': str(source), 'sha256': hashes}, indent=2) + '\n')
PY
docker image inspect --format '{{.Id}}' "$image" > "$output_dir/image-id.txt"
docker run --rm --platform linux/amd64 --entrypoint bash \
    --mount "type=bind,source=$output_dir,target=/work" "$image" -c '
    set -euo pipefail
    cd /work/build
    gfortran -O2 -std=legacy -ffixed-line-length-200 -c \
        /work/source/sophia/sophia2cpp.f \
        /work/source/sophia/SOPHIA/eventgen.f \
        /work/source/sophia/SOPHIA/sampling.f \
        /work/source/sophia/SOPHIA/inpoutput.f \
        /work/source/sophia/SOPHIA/jetset74dp.f \
        /work/source/sophia/SOPHIA/rndm.f
    g++ -std=c++17 -O2 -D_DEBUG \
        /work/source/sophia/Sophia.cpp /work/source/sophia/GitVersion.cpp \
        sophia2cpp.o eventgen.o sampling.o inpoutput.o jetset74dp.o rndm.o \
        -lgfortran -o sophiatab
    sha256sum sophiatab > binary-sha256.txt
    g++ --version > compiler-versions.txt
    gfortran --version >> compiler-versions.txt
' > "$output_dir/build.log" 2>&1
