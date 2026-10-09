# External Linux amd64 test runner

Calculations expected to exceed 30 minutes must run on an external server,
never on the workstation. This includes detached Docker calculations.

Prerequisites: native Linux x86-64, Python 3, Git, CMake 3.24+, make, GCC/G++,
gfortran, GSL and Boost filesystem/system development packages. See
../../docker/linux-amd64.Dockerfile for the current Ubuntu dependency list.
GitHub access is required for source and pinned GoogleTest downloads.

Copy this runner to the server, then run:

```sh
python3 run_tests.py --server --ref talys --data-dir /srv/mcray-test-data
```

It clones the remote, fetches the requested current ref once, checks out its
exact commit, builds Release and runs CTest. It uses a fresh checkout for every
run. Unpushed source changes are not available. The exact commit, runner hash,
commands, platform, logs, CTest XML/diagnostics and status are archived. The
manifest records input-file hashes supplied through repeated --input options.
For scientific runs, include data manifests, configurations and seeds; pin
reference programs independently. Compiler/library/image identities must be
captured by the workload if the comparison requires them.

Append a workload command after --, for example a reviewed scenario script:

```sh
python3 run_tests.py --server --ref talys --data-dir /srv/mcray-test-data \
  --input /srv/data/manifest.json -- \
  bash /srv/workloads/approved-comparison.sh
```

The workload runs only after build and CTest succeed. Its working directory is
the frozen source checkout. Write all intended result files below the absolute
MCRAY_RESULTS_DIR environment variable; explicitly pass that location to tools
that need an output argument. External data are not copied automatically.
The runner does not launch SOPHIA by default or supply missing generator sources.
Use tmux or a server scheduler to keep an approved long run alive.

Outputs in --data-dir are a unique timestamp/UUID .tar.gz and a small .json
summary containing its SHA-256. The work/<run-id>/ directory retains the source,
build and uncompressed outputs for diagnosis. Failed commands produce a nonzero
exit and are archived as failed, not passed. Keyboard interruption is recorded;
machine failure or SIGKILL can leave only the work directory and an unfinished
archive. No automatic resumption is provided. No existing results are overwritten.

DVC is optional and is not installed/configured by the runner. After completion,
register selected finalized bundles in an approved DVC metadata repository and
upload them to its configured remote. No Git commit/push, artifact upload or
credential changes are performed automatically. Do not track live output files.

Portable runner workflow checks (no compilation or long calculation):

```sh
python3 tests/server/test_runner.py
```
