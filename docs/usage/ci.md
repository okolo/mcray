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

The independent AI-review workflow remains a separate check. Its approval
is not evidence that a C++ build ran: that job previously lacked the required
development dependencies. Use the dedicated build/test check for execution
evidence, and the review for analysis. This change does not alter repository
branch-protection rules or make a new check required automatically.

These tests are implementation/regression checks, not full physical acceptance.
They do not start external-data research ensembles or table generation.
For local commands, dependencies and reference-data checks, see
[the test guide](../../tests/README.md).

Runner labels are documented in the
[GitHub-hosted runner reference](https://docs.github.com/en/actions/reference/runners/github-hosted-runners).
