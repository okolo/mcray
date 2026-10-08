# Linux amd64 checks in Docker

The primary mcray target is Linux x86-64 (Docker: `linux/amd64`).
Checks use Ubuntu 24.04, system GCC/G++/gfortran, GSL,
Boost.Filesystem, and Boost.System.
The application and dependencies use compatible libstdc++;
OpenMP uses GNU libgomp.

From the repository root:

```bash
bash bin/test-docker.sh
```

The script builds mcray-test-linux-amd64:local and then runs CTest in a
disposable container. Build and run explicitly select the platform;
the Dockerfile also checks dpkg --print-architecture = amd64.
On Apple Silicon, Docker Desktop emulates the tests; on Linux amd64 they
run natively. A running Docker daemon is required.
The first run needs access to Docker Hub, Ubuntu repositories, and GitHub
for the pinned GoogleTest archive.
Sources are copied into the image: building does not modify host files or results.
.dockerignore excludes local builds, IDE files, Git, and calculation results.

The container fully builds CRbeam, mcray_test, and z2t, then runs:

- 52 GoogleTest unit cases: interpolation for double and long double,
  including zeros, Clone, and legacy-call compatibility; PSB getRate
  boundaries, channel support, Clone behavior, weighted bookkeeping, and
  light-nucleus multiplicities; nucleon GZK event bookkeeping and proton pair
  weight/rate scaling, safe near-threshold rejection, and independent PPP-kernel quadrature at six photon energies;
  GZK spectral batch errors, Bonferroni quantiles, interval outcomes, unresolved-zero guards,
  and log-cell reconstruction with count conservation;
- two real-CMB PSB integration checks, a native Inoue12 extension-to-zero
  check, a static CMB+Inoue12 check, and a GZK
  rate-convergence regression;
- a two-step all-off propagation and output verification fixture;
- crbeam_help: CRbeam help startup, including static and OpenMP initialization;
- z2t_smoke: a short travel-time calculation for z=0.1.

C++ builds require C++17. GoogleTest 1.18.0 is pinned by archive SHA-256
and used only for tests. CTest labels are unit, integration, and e2e/smoke.
See ../tests/README.md for framework and assertion guidance.

The build uses two parallel jobs by default; locally override the count
through MCRAY_BUILD_JOBS. This limits memory use under amd64 emulation.
Each test has a 30-second timeout. Build/test failures produce a nonzero
script exit code.
These checks are not full physical PSB/TALYS validation; long historical
--test calculations are not run automatically.

Test selection and verbose output:

```bash
bash bin/test-docker.sh -L unit -V
```

Rerun without building:

```bash
docker run --rm --platform linux/amd64 mcray-test-linux-amd64:local
```

Check architecture and compiler versions:

```bash
docker image inspect mcray-test-linux-amd64:local --format '{{.Os}}/{{.Architecture}}'
docker run --rm --platform linux/amd64 --entrypoint bash \
    mcray-test-linux-amd64:local -c 'uname -m; gcc --version; gfortran --version'
```

Tests are registered with CTest and available outside Docker after a normal build:

```bash
ctest --test-dir build/Release --output-on-failure
```

Ubuntu is selected by its 24.04 tag; apt packages receive current fixes.
The image does not pin every package version and is not an immutable
environment archive.
For scientific comparisons, record compiler, image, and input-table versions
alongside the results.

## Numerical tolerances

On Linux amd64, long double has extended precision and differs from
long double on macOS arm64.
The first Docker run exposed an overly strict == 2 check in the ordinary
linear-interpolation test: long double returned 1.99999999999999999989158.
The check now uses 8 * epsilon * 2 tolerance; library computation is unchanged.
Checks of exact zero and coefficients returned directly from tables still
use exact equality.

Before the GoogleTest migration, verification dated 2026-10-01 showed that
the full Release build and all four CTest tests
passed in Docker Linux amd64 (uname: x86_64), emulated on macOS arm64.
Tools: GCC/G++/gfortran 13.3.0, Boost 1.83, GSL 2.7.1, Ubuntu 24.04.
All four tests also passed natively on macOS arm64 with Clang 22.1.5.
A native run on physical Linux amd64 hardware has not yet been performed.

## TransportCR comparison

Run the slower isolated Fe-56 PSB comparison separately:

```bash
bash bin/test-transportcr-docker.sh
```

This builds a temporary clone of the pinned TransportCR reference using its
unchanged bin/build.sh, and compares static-CMB propagation with 200,000 mcray
primary histories in 20 independent batches. No TransportCR changes are committed.
Results, input XML, logs, reference revision, and GoogleTest XML are saved in a
new build/transportcr-validation.* directory. See ../docs/validation/psb_validation_results.md
for the physical setup, convergence checks, and remaining limitations.
This long comparison is not part of the default short CTest suite.

After migration to C++17 and GoogleTest 1.18.0, all 26 CTest cases passed
in both environments on 2026-10-01: 24 unit cases and two command-line
smoke checks. The full Release builds also succeeded.
