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

The short suite contains 37 discovered GoogleTest unit cases, three
convolution integration cases, two command-line smoke tests, and a two-step
no-interaction E2E fixture (44 CTest entries). Unit tests exercise the production library:

- Log interpolation for double and long double: zeros, positive intervals,
  external limits, Clone, legacy behavior, and invalid scale combinations.
- PSB branching: invalid deltaA boundaries, unchanged valid coefficients,
  and single-/double-nucleon channels.
- Nucleon GZK: supported rates, Clone, and SOPHIA event energy, weight and
  nucleon-number bookkeeping.
- Proton pairs: sampling weight versus rate scaling, Clone, continuous loss,
  and explicit absence of nuclear/neutron support.

CTest labels distinguish unit, integration, and e2e. CRbeam help and the
short z2t calculation also have the smoke label. Every current test has a 30-second timeout.
These smoke checks do not validate nuclear propagation.

## Adding checks

TransportCR GZK cases can use externally supplied finer SOPHIA tables:
pass --sophia-tables /path/to/sophia2 to tests/transportcr/prepare_case.py.
The helper validates the selected density (10 times --resolution), copies
the required files into the isolated case, and records their SHA-256 values
in sophia-manifest.json. It does not regenerate or overwrite reference data.

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

The approved nuclear validation plan uses a dedicated C++ scenario and
GoogleTest comparison executable, coordinated by CTest or a small runner.
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

## PSB propagation against TransportCR

See ../nuclear_validation_plan.md and ../psb_validation_results.md.
The real production PropagationEngine is used by nuclear_scenario; no
transport algorithm is duplicated in the tests.

```bash
# Pin and build TransportCR separately in the Linux amd64 environment.
bash bin/test-transportcr-docker.sh /absolute/path/to/transportcr

# Run locally using an already built TransportCR; choose a new output directory.
bash tests/transportcr/run_psb_comparison.sh \
    /absolute/path/to/transportcr /absolute/path/to/new-results
```

The default slow comparison uses 20 independent batches of 10,000 Fe-56
primaries, CMB at 2.7255 K, and the lookback interval of z=0.01 in a static
background. Outputs include per-primary weighted species yields and energies,
XML inputs, initial/final TransportCR spectra, logs, and GoogleTest XML.
Primary-history counts and baryon accounting are checked so caught propagation
errors cannot silently turn missing histories into a successful comparison.

The Docker wrapper copies only the pinned committed TransportCR revision
94eca0c886df80b1d0b1dd09dae95b0483c2273c into a temporary build context.
It does not change or commit in the reference checkout. Large outputs stay in
a unique build/transportcr-validation.* directory. The primary propagation
CLI CRbeam still lacks nuclear input; this is a library-based E2E check.

The XML generator defaults to a broad grid. For PSB-only runs, its
--grid-span-decades option can raise Emin without changing the retained
energy-per-nucleon distribution. TransportCR's fixed-energy source overrides
BaseEmax with its total injection energy, so XML Emax does not fully narrow
this grid; check the actual rates/output range. Initial/final output uses
interpolated total-energy columns. Dividing j(E)*E^2 by E and normalizing to
the initial Fe column removes common unit/bin-width/volume factors; this
conversion is validated with the all-off run. Record grid convergence before
treating the comparison as a completed physical validation.
