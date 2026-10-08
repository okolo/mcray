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

The short suite contains 52 discovered GoogleTest unit cases, five
background and convolution integration cases, two command-line smoke tests, and a two-step
no-interaction E2E fixture (61 CTest entries). Unit tests exercise the production library:

- Log interpolation for double and long double: zeros, positive intervals,
  external limits, Clone, legacy behavior, and invalid scale combinations.
- PSB branching: invalid deltaA boundaries, unchanged valid coefficients,
  and single-/double-nucleon channels.
- Nucleon GZK: supported rates, Clone, and SOPHIA event energy, weight and
  nucleon-number bookkeeping.
- Proton pairs: sampling weight versus rate scaling, Clone, continuous loss,
  explicit absence of nuclear/neutron support, near-threshold null
  interactions, and independent rest-frame integration of the PPP
  differential kernel compared with lab-energy sampling at photon
  energies 2.2, 3, 5, 10, 20, and 50 electron masses (10,000 samples each).

The scenario runner accepts --background cmb-inoue for the agreed static
CMB+Inoue12 Baseline experiment. FrozenInoueBaseline is a test-only fixture;
production Inoue12 already extends the first table plane to z=0; the fixture
only freezes its density. Integration tests check the native extension to z=0,
the compound background, and rates along the path.
TransportCR prepare_case.py accepts the same --background option and records
the EBL table hash. The production background implementation is unchanged.

Use --particle-output in nuclear_scenario to retain individual weighted
particles within each complete primary history. The same TSV columns remain;
energy_eV is weighted energy, so energy_eV/yield is a particle's energy.
With this output, nuclear_comparison --spectral-report /path/report.tsv writes
diagnostic spectra in centered logarithmic bins (10 per decade, six decades
below injection, separate underflow/overflow). Its errors use independent
batch means. These bins do not establish a global spectral acceptance result.
For zero observed products, the reported SE is unavailable (NA); a 99%
upper bound applies to the probability that a primary history populates the
bin, not to weighted mean multiplicity. Do not interpret missing products
as a precise physical zero.
GZK secondary count/energy moments are checked separately with GoogleTest.
Pair counts and shapes use different models in the two codes; compare energy
balance and check the mcray sampler independently rather than requiring
identical pair spectra.

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

For completed pair-only TransportCR runs, compare grid or step refinement:

```bash
python3 tests/transportcr/check_pair_convergence.py \
    --coarse /path/to/coarse-case --fine /path/to/fine-case
```

The script verifies compatible physical settings and lower energy bounds,
then compares proton loss and emitted pair energy. It writes JSON and returns
1 when either relative change exceeds 0.002. Grid resolution, spatial step,
and the empty upper energy range may differ; derived source normalizations
written back by TransportCR are ignored. This is a successive-refinement
check, not a proof of absolute error or of spectral-shape convergence.


## Approved GZK spectral comparison

See ../gzk_spectral_validation_plan.md. Seven secondary species in four
fixed energy intervals form a family of 28 observables. The approved method
uses 99% simultaneous Bonferroni Student intervals and a 5% spectral budget;
existing integrated moment checks remain unchanged.

The shared test-only gzk_spectral_statistics.h computes batch errors,
quantiles, interval classification, and fixed-bin boundaries. Eight GoogleTest
cases check correlated history contributions, analytic Cauchy quantiles,
tail handling, and agreement/disagreement/inconclusive outcomes.

Use the comparison runner's pilot report after preparing compatible
TransportCR grid and spatial-step controls:

```bash
build/Release/tests/nuclear_comparison \
    --process gzk --a 1 --energy-per-nucleon 1e20 \
    --batches 20 --primaries 7500 --mc /path/to/pilot.tsv \
    --initial /path/to/reference/step0/all --final /path/to/reference/fin/all \
    --gzk-grid-initial /path/to/grid-control/step0/all \
    --gzk-grid-final /path/to/grid-control/fin/all \
    --gzk-step-initial /path/to/step-control/step0/all \
    --gzk-step-final /path/to/step-control/fin/all \
    --gzk-spectral-report /path/to/wide-bin-pilot.tsv
```

The report is exploratory, not a global spectral acceptance assertion. It
includes underflow/overflow, history occupancy, batch SEs, simultaneous
half-widths, refinement changes, and conditional interval outcomes. Missing
controls cannot imply numerical convergence. Physical XML settings must
be verified compatible before supplying the spectra. The pilot temporarily
uses the sum of absolute grid/step changes as a practical refinement envelope;
this is not a rigorous error bound. Conditional outcomes require further
assessment of Student-interval applicability. Filling all batches alone
is not evidence of normality.

Reference products use the approved constant dN/dlog(E) reconstruction within
centered logarithmic cells. Counts are apportioned by log-energy intersection
with physical intervals; both tails are retained. The parser validates the
printed log grid. Three additional unit cases verify crossing-cell fractions,
count conservation, and rejection of missing grid cells. This reconstruction
is limited to the seven secondary species. Do not interpret the pilot as a
completed physical agreement test or use it to select favorable seeds.


## Independent SOPHIA table generator

The original TransportCR predecessor generator was found at
~/propagation/sources/sophia. See ../sophia_table_generation.md for source
provenance, runtime estimates, and pending full-statistics decisions.
Build an immutable working copy outside Git using the existing Linux amd64
image (build it with bin/test-docker.sh first):

```bash
bash tests/transportcr/build_sophia_generator_docker.sh \
    /absolute/path/to/propagation/sources /absolute/path/to/new-output
```

The script refuses an existing output directory, hashes all copied source
files, records the image identity and compiler versions, and builds with
assertions enabled. It does not modify the original sources or tables.
The resulting build/sophiatab is a Linux amd64 executable; run it in the
same Docker image, in a separate output directory. Its CLI is:

```text
sophiatab p|n bins_per_decade kmax_GeV trials_per_cell seed
```

Small pilot tables validate format, reproducibility, and loader integration;
they are never numerical-convergence or physical-acceptance references.
Full table statistics and MC acceptance statistics are separate decisions.


After both generation jobs succeed, validate them before copying into a
TransportCR case. Each directory must contain exit-code.txt with value 0:

```bash
python3 tests/transportcr/validate_sophia_tables.py \
    --p-directory /path/to/full200-p --n-directory /path/to/full200-n \
    --bins-per-decade 200 --trials-per-cell 1000000 --kmax-gev 2e6 \
    --manifest /path/to/new-validated-manifest.json
```

The manifest is created exclusively after validation and never overwrites an
existing one. Passing this check establishes format/completeness, not physical
or statistical accuracy. Interrupted generation has no automatic resume.


## PSB background-grid convergence controls

The scenario runner accepts `--background-step-decades` (default 0.01),
which sets the requested logarithmic photon-energy grid spacing in
ContinuousBackgroundIntegral. The constructor adjusts that spacing to cover
the background bounds. The output records the requested spacing separately
from `--accuracy`, the quadrature relative tolerance. These controls change
numerical resolution, not the frozen-background physical assumptions.

For the extended static CMB+Inoue12 PSB checks, use Fe56, L=1 Mpc,
E/A=1e18 or 3.1622776601683794e18 eV, 20 batches of 10000 histories,
accuracy=1e-4, and spacing=0.01. Rates are additionally probed with spacing
0.005/0.0025 and tolerance 1e-5. The tighter tolerance succeeds for Fe56 and
A30, but still raises GSL_EROUND for A45 at the higher energy on macOS
and Linux amd64; it is not validated for the entire daughter chain.

The all-daughter integration regression uses the working tolerance 1e-4,
compares K/S rates and spacing 0.01/0.005 to within 0.05%, and additionally
checks tolerance 1e-5 for the specific Fe56/A30 regression points. Read
psb_validation_results.md for statistical limitations, reference-grid
controls, and raw-run locations. A zero observed rare yield and zero empirical
SE do not establish relative precision.
