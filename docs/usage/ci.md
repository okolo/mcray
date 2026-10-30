# Automated build and short tests

The [C++ build and tests workflow](../../.github/workflows/ci.yml) runs on pull
requests, pushes to main, and manual dispatch. Its Linux amd64 Release job uses
an Ubuntu 22.04 hosted runner, installs the documented C++/Fortran, Boost, GSL,
yaml-cpp and OpenSSL development dependencies, then configures the supported
src/app/crbeam CMake project with tests enabled and builds all targets.

The job runs the complete short CTest suite, including the worker lifetime
regression, then both postprocessing Python suites and the review-publication
guard tests. MCRAY_CRBEAM points to the freshly built executable, so the loader's
executable round-trip case is exercised rather than skipped. NumPy 1.26.4 and
PyYAML 6.0.3 are installed in a project-local venv; NumPy's minimum supported
postprocessing version is declared in tools/postprocess/requirements.txt.

A separate Linux amd64 ASan+UBSan job instruments C, C++ and Fortran code
with GCC AddressSanitizer and UndefinedBehaviorSanitizer in RelWithDebInfo
(`-O1 -g`, frame pointers retained). Executables use non-PIE linking to avoid
ASan shadow-memory/ASLR collisions. It builds every target and runs the full
short CTest suite with leak detection enabled and both sanitizers configured
to stop on the first error. Reports flow into CTest logs and the retained
JUnit/diagnostic artifact; no suppressions or excluded tests are configured.
System libraries are not rebuilt with instrumentation. This job does not run
Python suites, which remain in the Release job, or ThreadSanitizer; ASan/UBSan
success is not evidence that arbitrary data races are absent.

The checked revision, architecture, compiler/CMake and Python package versions
are logged. PR runs check GitHub's proposed merge revision; the actual SHA is
in environment.log. GoogleTest remains pinned by its release archive hash in
tests/CMakeLists.txt. Bash steps use pipefail so piping through tee cannot hide
a failing build or test command. Checkout/artifact actions are pinned to full
commit SHAs. The job has read-only repository permissions, does not use provider
secrets, and does not publish GitHub reviews or push changes. Logs, the resolved
Python requirements and CTest JUnit results are retained as an artifact for
14 days, including on failure. Large build outputs and numerical tables are
not uploaded. Each job has a 25-minute timeout and two build/test workers.

## Independent review execution evidence

The independent AI-review workflow remains a separate check. It runs on Ubuntu
22.04, installs the same C/C++/Fortran, Boost, GSL, yaml-cpp and OpenSSL development
dependencies, and checks them by compiling/linking a small toolchain probe.
The probe uses workflow-owned source outside the checkout; it does not configure
or execute the candidate project's CMake code. System compiler paths and /usr
package prefix avoid accidentally selecting an unrelated Conda installation.
Control scripts are copied and tested before checking out the frozen head,
so a PR branch predating this workflow need not contain the new helpers.
The reviewer receives a scratch build path and must keep the frozen head clean.

Before analysis, a read-only collector waits up to ten minutes for the matching
pull-request CI run. It validates the repository, workflow, frozen PR head/base,
latest run and current attempt, both required jobs and available CTest step metadata, and
current-attempt nonexpired artifacts. It reads only bounded environment/JUnit
members without extracting arbitrary archive paths. The environment's actual
build SHA must be the reviewed head or a merge commit whose parents are exactly
that frozen base and head. The build tree, report hashes, counts/skips and job
URLs are recorded. Reruns or PR changes during collection invalidate evidence;
an older success cannot replace a newer failed/pending run. A changed CI
workflow definition is reported as unavailable for automatic evidence reuse
and needs independent assessment/local checks.

GitHub sometimes returns an empty step list for a completed job. This is
explicitly recorded as `step_metadata_available: false`; the successful job,
unchanged CI definition and revision-bound JUnit/environment reports remain
required. When a nonempty step list is provided, a missing or unsuccessful
CTest step invalidates the evidence.

Verified means provenance and report consistency, not code correctness or
physical acceptance. Logs/reports are untrusted data, never model instructions.
The reviewer must independently analyze the diff and required architecture
consultation, assess relevant coverage, and distinguish CI execution from
checks it personally ran. Missing/failed/stale evidence is never success;
the installed toolchain also allows targeted project checks. Full duplicate
sanitizer builds or unavailable historical-fork data are not mandatory when
verified CI provides the relevant execution evidence. Material unresolved
verification gaps still require COMMENT; confirmed blockers require
REQUEST_CHANGES. Green checks never automatically produce APPROVE.

The formal publisher's head/base, completed-review, blocking-finding and
architectural-consultation gates are unchanged, including the check after
approval publication. Branch protection is unchanged. The job adds only
`actions: read` for CI evidence access; its existing review permissions remain.
Toolchain logs, evidence JSON and the two selected reports are retained for
14 days. Availability errors remain explicit in the evidence file. The review
job's success means execution/publication completed, not necessarily APPROVE.

The collector and existing publication guards have offline standard-library
tests, included automatically in the Release CI Python checks:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s .github/scripts -p 'test_*.py' -v
```

These tests are implementation/regression checks, not full physical acceptance.
They do not start external-data research ensembles or table generation.
For local commands, dependencies and reference-data checks, see
[the test guide](../../tests/README.md).

Runner labels are documented in the
[GitHub-hosted runner reference](https://docs.github.com/en/actions/reference/runners/github-hosted-runners).
