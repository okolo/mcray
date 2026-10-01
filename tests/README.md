# Testing mcray

The application and library require C++17. C++ unit and integration tests
use GoogleTest; CTest is the common test runner for local and Docker builds.
GoogleTest 1.18.0 is fetched only when BUILD_TESTING is enabled (the default).
Its release archive and SHA-256 are pinned in CMakeLists.txt.
CMake 3.16 or newer is required for test configuration; bin/build.sh uses
--fresh and requires CMake 3.24 or newer.

## Run tests

From the repository root:

```bash
bash bin/build.sh Release
ctest --test-dir build/Release --output-on-failure
ctest --test-dir build/Release -L unit --output-on-failure
bash bin/test-docker.sh
bash bin/test-docker.sh -L unit
```

The current suite contains 24 discovered GoogleTest unit cases and two
command-line smoke tests. Unit tests exercise the production library:

- Log interpolation for double and long double: zeros, positive intervals,
  external limits, Clone, legacy behavior, and invalid scale combinations.
- PSB branching: invalid deltaA boundaries, unchanged valid coefficients,
  and single-/double-nucleon channels.

CTest labels are unit for GoogleTest cases and e2e;smoke for CRbeam help
and the short z2t calculation. Every current test has a 30-second timeout.
These smoke checks do not validate nuclear propagation.

## Adding checks

Use GoogleTest for new C++ unit and integration checks. Unit checks isolate
a contract or numerical operation; integration checks combine data loading,
convolution, and processes with controlled inputs. Label new integration
checks integration and keep their inputs and tolerances explicit.
Physical scenarios and acceptance criteria still require agreement.

Check exact zeros and unchanged table coefficients with exact assertions.
For numerical computations, state an absolute/relative tolerance and retain
the tested type's precision, especially for long double. For Monte Carlo
outputs, compare distributions using statistical uncertainties; a fixed
seed alone does not establish cross-platform or cross-code equivalence.

The framework for richer end-to-end output checks has not yet been selected.
Historical --test calculations remain available separately and are not
automatically added to the short suite.

## Dependency setup

The first test configuration requires GitHub access. CMake caches the
download in the build directory. To configure offline, provide an extracted,
verified GoogleTest 1.18.0 source tree:

```bash
cmake -S src/app/crbeam -B build/offline \
    -DFETCHCONTENT_SOURCE_DIR_GOOGLETEST=/absolute/path/to/googletest-1.18.0
```

Include the compiler and dependency options appropriate to the platform;
bin/build.sh contains the supported macOS and Linux settings.
When using an externally supplied source tree, verify it against the pinned
archive yourself. Production-only configurations may use BUILD_TESTING=OFF
to avoid downloading or building GoogleTest.

The test executable links GoogleTest's main before mcray because the legacy
library also contains a main object that delegates to user_main.
