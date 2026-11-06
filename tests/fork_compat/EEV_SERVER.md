# One-command EeV comparison on Ubuntu 22.04

Use an external **Ubuntu 22.04 x86-64** server. The entry point is
**run_eev_field_server.py**, a standalone Python-standard-library script.
It prepares a fresh checkout, builds, tests and runs the
[frozen 162-million-history comparison](EEV_FIELD_ENSEMBLE.md).
It defaults to sixteen worker processes and enforces **eight hours total**
from invocation, including source/reference retrieval, compilation, checks,
collection, archive verification and metadata packaging. Package installation
and transferring files happen beforehand. No root access is needed for the run.

Requirements: Internet access for the repository and pinned reference/GoogleTest
archives (or the transferred offline inputs below), at least 64 GiB free in
the test-data filesystem for final collection, and preferably sixteen CPU
cores. Workers are isolated processes, with OMP/BLAS threads fixed to one.
Fewer cores are supported with --workers 1..16; counts/seeds do not change.
Do not run the final campaign on the workstation.

## Install dependencies once

```sh
sudo apt-get update
sudo apt-get install -y --no-install-recommends \
  build-essential cmake gfortran git curl ca-certificates python3 \
  libgsl-dev libboost-filesystem-dev libboost-system-dev \
  libyaml-cpp-dev libssl-dev tmux
mkdir -p "$HOME/eev-runner" "$HOME/mcray-test-data"
```

Ubuntu 22.04's system CMake/GCC/Python are sufficient. The script downloads
GoogleTest 1.18.0 and verifies its pinned SHA-256 before configuring CMake.
The historical snapshot is pinned to
e1fc71f613efe7fecf3210e15339c34814245cd4; the adapter verifies its seventeen
source-file hashes and shared unit/time assumptions before compilation.
No external nuclear or interaction tables are required by this test.

## Transfer the current unpublished candidate

A Git bundle makes this workflow usable before a push. On the development
machine, from the public mcray repository root, after committing intended
changes and checking the tree is clean:

```sh
mkdir -p build/eev-transfer
revision=$(git rev-parse HEAD)
git bundle create build/eev-transfer/candidate.bundle codex/fork-field-observer
cp tests/fork_compat/run_eev_field_server.py build/eev-transfer/
printf '%s\n' "$revision" > build/eev-transfer/REVISION
scp build/eev-transfer/candidate.bundle build/eev-transfer/REVISION \
  build/eev-transfer/run_eev_field_server.py \
  USER@SERVER:eev-runner/
```

The bundle contains the branch's reachable committed history. It does not
include ignored build/results or uncommitted files and does not push anything.
The script clones the bundle into a unique run directory and verifies the
exact supplied commit. Keep the bundle and REVISION paired. If the candidate
is already published, --repository may instead be its Git URL; still supply
an exact forty-character commit, not a moving branch name.

## Start the complete eight-hour run

On the server:

```sh
tmux new -s mcray-eev
cd "$HOME/eev-runner"
python3 run_eev_field_server.py --server \
  --repository "$PWD/candidate.bundle" \
  --revision "$(cat REVISION)" \
  --data-dir "$HOME/mcray-test-data" --workers 16
```

Detach with Ctrl-B then D; return with tmux attach -t mcray-eev. The server's
network/session/scheduler must permit the job to run. The script does not
change the system or upload artifacts. Do not start a second final run to
retry a disappointing scientific outcome or a partial attempt.

To check the entire setup without consuming final seeds, add **--fixture-only**.
That still clones, builds, runs CTest/Python guards and the disjoint fixture;
it never starts the 162-million-history collection. The complete command
already performs this prerequisite, so a separate fixture-only invocation is
optional. Fixture outputs cannot be interpreted as acceptance evidence.

If Internet access is unavailable, transfer these two archives too:

- The immutable reference archive from
  https://codeload.github.com/okolo/mcray/tar.gz/e1fc71f613efe7fecf3210e15339c34814245cd4.
- GoogleTest from
  https://codeload.github.com/google/googletest/tar.gz/refs/tags/v1.18.0,
  SHA-256 6e3191c1455468b3fc35a417fb565c1c5071aee1b7e7f85e30cf48a98d37d8b5.

Add --reference-archive /absolute/reference.tar.gz and
--gtest-archive /absolute/googletest-1.18.0.tar.gz to the same command.
All reference extraction rejects links/path traversal; reference inputs remain
read-only after extraction. Ordinary CMake tests do not need the old snapshot.

## Results and interruption

The script prints a unique directory under --data-dir:

- server-status.json: exact commit, host, elapsed time, outer status/exit code.
- build.log, ctest.log, ctest.xml, python-guards.log and fixture.log: prerequisites.
- fixture/: separately seeded software-check records and verified summaries.
- ensemble/report.json: the 144 simultaneous intervals and scientific result.
- ensemble/*-*.jsonl.gz: eighteen raw per-primary streams; *.verified.json
  records their raw/compressed hashes and counts after full decompression.
- ensemble/manifest.json and metadata.tar.gz: source/reference/input/core
  hashes, compiler command/version, reservation state, executable and logs.
- server-metadata.tar.gz: outer logs/status, without duplicating raw streams.

A completed collection may report compatible_at_reported_resolution=false:
that is a detected discrepancy, not a reason to rerun or change the protocol.
Compatibility means only zero-containing difference intervals for the frozen
observables at the stated resolution. Numerical stopping differences remain
separate diagnostics; whole-program physical equivalence is not asserted.

The stable eev-field-ensemble-v1.reservation.json in --data-dir is created
exclusively before the first final source process. It blocks a repeated final
reservation even with a different output directory. Keep it and the data;
do not bypass it by choosing a new storage directory. Interrupted/failed final
attempts remain spent and need separately reviewed recovery.

Timeout/error returns nonzero and retains partial data/logs. The supervisor
stops descendants, including workers in separate process sessions. A hard
eight-hour watchdog covers final verification/packaging too. A machine failure
or SIGKILL can leave incomplete gzip/metadata; such output is never a successful
comparison. Ordinary script interruption is recorded as failure. There is no
automatic resume, history replacement or implicit increase of sample count.
