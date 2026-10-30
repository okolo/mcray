# Testing mcray

The application and library require C++17. C++ unit and integration tests
use GoogleTest; CTest is the common test runner for local and Docker builds.
GoogleTest 1.18.0 is fetched only when BUILD_TESTING is enabled (the default).
Its release archive and SHA-256 are pinned in CMakeLists.txt.
CMake 3.16 or newer is required for test configuration. `bin/build.sh` uses
`--fresh` when available and clears the CMake cache on older versions,
including Ubuntu 22.04's CMake 3.22.

## Run tests

For fresh remote checkouts and archived results on an external Linux amd64
server, see [the server runner](server/README.md). Never launch calculations
expected to exceed 30 minutes on the workstation.

From the repository root:

```bash
bash bin/build.sh Release
ctest --test-dir build/Release --output-on-failure
ctest --test-dir build/Release -L unit --output-on-failure
bash bin/test-docker.sh
bash bin/test-docker.sh -L unit
```

The short CTest suite includes nuclear,
EM-secondary, configuration, background/convolution and no-interaction E2E
checks. Use `ctest --test-dir build/Release -N` for the exact inventory of
the selected revision. Unit tests exercise the production library:

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

## Grid magnetic fields

The `grid_magnetic_field_tests` target covers validated Cartesian grids,
trilinear interpolation, boundary extensions, component ownership, concurrent
reads and existing deflector clones. It uses only temporary synthetic input
and runs through CTest with prefix `unit.magnetic.` and label `unit`.
See [the public API contract](../docs/usage/grid_magnetic_fields.md).

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

## Installer orchestration checks

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests/installer -v
```

These short checks exercise the real shell installer with isolated command
and generator substitutes: no apt changes, network downloads, C++ builds or
physical calculations run. They check preserved/corrupt packages, pinned
source selection, interruption boundaries, locking, destination races,
paths with spaces and repeat installation. Temporary fixtures stay under
build/installer-tests in this checkout. They do not establish a native Ubuntu
installation or physical validation of prepared tables.

## PSB propagation against TransportCR

See [the PSB validation results](../docs/validation/psb_validation_results.md)
for the tested scope, comparison method and remaining limits.
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

See [the GZK method and pilot results](../docs/validation/psb_validation_results.md#approved-gzk-spectral-method-and-first-pilot).
Seven secondary species in four
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
~/propagation/sources/sophia. See ../docs/validation/sophia_table_generation.md for source
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
../docs/validation/psb_validation_results.md for statistical limitations, reference-grid
controls, and raw-run locations. A zero observed rare yield and zero empirical
SE do not establish relative precision.


## CRbeam nuclear source input

CRbeam accepts `--particle nucleus --A 56 --Z 26`. Both nuclear numbers are
unsigned integers with defaults 0. `--A 56` alone, or `--A 56 --Z 0`, chooses
Z=26 through the PSB-chain representative table. Source redshift still uses
lowercase `--z`. Negative integers, fractional/scientific notation, trailing
characters, missing values, and unsigned overflow are rejected.

The first-stage source is limited to Fe56. A nucleus requires A>=2; explicit
Z must not exceed A, and automatic selection requires a PSB-table entry.
The resolved A,Z and Z-selection method are printed. For other particles,
A,Z values have no physical effect, although their unsigned syntax is checked.
Existing numeric `--primary` remains supported. Named `--particle` accepts
existing Particle::Name values such as proton, neutron, photon, electron and
positron. Supplying both selectors is an error.

The earlier CLI-only stage was verified with 76 CTest entries on macOS and
Docker Linux amd64. The unified migration now removes its nuclear run guard
and connects separate PSB and TALYS processes; see ../docs/usage/talys_usage.md.

## Unified nuclei and TALYS checks

Short GoogleTest cases include identity/mass copying and mutation errors,
PSB/TALYS mass-model separation, isotope output and trajectory columns,
synthetic data/network validation, exclusive-channel selection, conservation,
clone behavior and failed-cascade isolation. CRbeam integration cases exercise
actual Off and PSB nuclear runs, along with CLI validation.
External real TALYS data are not bundled into the short suite or Docker image.

Prepare and verify the pinned data as described in ../docs/usage/talys_usage.md. From bin,
run an isolated frozen-background scan and independent Planck convolution:

```sh
../build/Release/tests/nuclear_scenario --process talys \
  --talys-data ../../validation/talys-fe56-34b10e2 --scan \
  --output ../../validation/talys-cmb-scan.tsv
python3 ../tests/talys/check_cmb_rates.py \
  --data ../../validation/talys-fe56-34b10e2 \
  --scan ../../validation/talys-cmb-scan.tsv \
  --output ../../validation/talys-cmb-independent.tsv
```

The independent script uses analytic CMB density and composite Simpson
quadrature in log(K), with its own log/zero interpolation. It reports
numerical differences and reference convergence without declaring a physical
TALYS acceptance threshold. mcray's existing unit constants are held fixed.

An isolated E2E control and GoogleTest verification, from bin:

```sh
../build/Release/tests/nuclear_scenario --process talys \
  --talys-data ../../validation/talys-fe56-34b10e2 --distance-mpc 1 \
  --energy-per-nucleon 3.9810717055349857e18 --primaries 1000 --batches 1 \
  --output ../../validation/talys-cmb-e2e.tsv
../build/Release/tests/nuclear_comparison --process talys \
  --mc ../../validation/talys-cmb-e2e.tsv --primaries 1000 --batches 1 \
  --energy-per-nucleon 3.9810717055349857e18
```

Use --background cmb-inoue for the agreed frozen Inoue12 background. The
TALYS E2E checks assert per-history A,Z and energy conservation, complete
history counts, and Fe56 survival against the static attenuation law using
the existing four-standard-error check. They do not compare with TransportCR
(which implements PSB) or establish rare-product statistical precision.

Real-data controls and migration comparison artifacts are outside Git in
../validation/unified-nucleus-20261008. See ../docs/usage/talys_usage.md for results,
physical limitations, data checksum verification and future comparisons.

Final migration verification on 2026-10-08: all 88 short CTest entries pass
on macOS arm64 (28.87 s) and Docker Linux amd64 emulation (94.90 s).
The real-data TALYS CMB and CMB+Inoue12 controls pass all three applicable
GoogleTest E2E checks on macOS; the CMB control and two-thread CRbeam TALYS
run also pass on Linux amd64. These results retain the physical-validation
limitations stated above and in ../docs/usage/talys_usage.md.

## Nuclear decay controls

The current short suite contains 96 CTest entries. NuclearDecay tests cover
lifetimes, branches, beta leptons, compound channels, emission bookkeeping,
preflight errors and exponential survival in the transport engine. See
../docs/usage/nuclear_decay.md for external real-data scenarios and data preparation.
nuclear_decay_comparison runs GoogleTest assertions on isolated scenario
outputs, including charge with leptons; exact beta event energy conservation
is not asserted under the approved CRPropa kinematics.


## Nuclear PPP

`NuclearPairs` unit cases exercise equal-gamma scaling using tabulated masses,
secondary energy/weight and ancestry, unchanged parent identity/energy,
thresholds, unsupported species and clone behavior. CLI tests validate the
Off/Losses/Full modes, tabulated-mass requirement and EM/PSB exclusions.
`nuclear_scenario --process nuclear-pair-cel|nuclear-pair` provides isolated
loss/secondary diagnostics without redshift losses or EM transport. See
[the public contract](../docs/usage/nuclear_pair_production.md) and
[regression results](../docs/validation/nuclear_transport_baseline.md).

For a one-step GZK spectator data inventory (not full reachability), run:

```sh
python3 tools/audit_nuclear_gzk_residuals.py --data /path/to/combined-bundle \
  --output /tmp/gzk-spectators.json
```

The output must not already exist. Exit 1 means coverage gaps were found;
it does not authorize stability assumptions, missing-table exceptions or
artificial dripping. The report hashes its four input tables.


## Nuclear GZK and policy comparison

`NuclearGzk` cases check modifiers at A=8/9, equal-gamma target rates,
free-nucleon regression, clones, pure-spectator breakup, both mass/decay
policies, missing data and unsupported modes. Two hundred sampled events
check baryon number, charge including leptons, energy (relative tolerance
1e-6), inherited weights and target frequencies (five binomial standard
errors). These are short implementation checks, not spectral acceptance.

`nuclear_scenario --process nuclear-gzk --gzk-residuals Strict|CRPropa`
provides isolated GZK plus nuclear decay, with redshift losses, EM transport
and free-neutron decay disabled. It accepts `--a`, `--charge` and the common
`--talys-data` bundle. Add `--preflight` to write the reachable identity count
and metadata without transport. Counts include pure-spectator markers and
free nucleons. Use identical primary histories/configurations for effect
studies. See [the public GZK contract](../docs/usage/nuclear_gzk.md).


## Nuclear photons

`NuclearPhotonsTest` covers native conditional convolution, Strict/CRPropa
saturation (including a zero denominator), Bernoulli line multiplicities,
rest-frame boosts, weights/ancestry, clone settings, missing coverage,
explicit exceptions, compound-decay line association and photon-free RNG
regression. CLI checks cover defaults, process requirements and creation with
`--noEMcascade` filtering. Statistical checks use the stated five-standard-error
bounds; they do not establish astrophysical spectral acceptance.

`nuclear_scenario` accepts `--photodisintegration-photons On|Off`,
`--decay-photons On|Off` and `--pd-photon-yields Strict|CRPropa`. Its photons
are retained without EM propagation. `nuclear_photon_probe DATA A Z E_PER_A_EV
Strict|CRPropa` reports native CMB conditional yields for one isotope, without
full-network transport or transition-frequency weighting. Both use the
[public photon data contract](../docs/usage/nuclear_photons.md). Failed
propagation and incomplete probe outputs must not be used as complete results.

## Integration performance comparison

`benchmark_nuclear_integration.py` compares two frozen executable directories
containing `CRbeam` and `nuclear_scenario`. Use the stable knot-split baseline
and optimized implementation built with the same Release toolchain. CRbeam
must support `--seed` in both variants. Copy binaries outside Git before
rebuilding; freeze the corresponding source patches and hashes separately.

```bash
python3 tests/benchmark_nuclear_integration.py \
  --before-dir /results/frozen-before --after-dir /results/frozen-after \
  --data /data/photon-enabled-nuclei --output /results/new-benchmark
```

Four cases cover Fe56 TALYS+decay with both photon sources on CMB and
CMB+Inoue12: complete nuclear histories at 1 Mpc (no EM interactions), and
production CRbeam at z=2.5e-5 with its EM cascade enabled. CRbeam uses a
1e11 eV tracking cutoff and thinning off. Nuclear PPP/GZK remain off.
Each variant receives a warmup, then three paired 64-primary runs with
identical pair seeds; the execution order alternates. OMP_NUM_THREADS=1,
data/binary hashes, commands, timings and species counts/weighted energies
are recorded. Times include setup, propagation and output. The script uses
a 240-second job timeout and 1800-second total local limit; a timeout retains
an incomplete result rather than silently reporting a speedup. These are
performance diagnostics, not a statistical spectrum-acceptance campaign.
Avoid concurrent builds or calculations when collecting timings.

## Independent deterministic PD reference

The [deterministic controls](deterministic/README.md) propagate PD-only Fe56
first moments with independent, production-exported and official coefficients,
including analytic secondary-photon rewards. Run their Python analytic/rule
checks with `python3 -m unittest discover -s tests/deterministic -v` in the
NumPy/SciPy environment. The validation-only `deterministic_pd_export` CMake
target exports production coefficients without sampling histories. These tools
also provide [combined PD+decay energy-state transport](deterministic/COUPLED.md)
with explicit grid and stiff-solver controls. Neither tool issues acceptance
verdicts.

## Python postprocessing checks

Run `PYTHONDONTWRITEBYTECODE=1 python3 tests/test_mc_rotate.py` from the
repository root with NumPy installed. These 18 deterministic checks compare
sphere-extrapolated delays with analytic and 100-digit Decimal references,
exercise tiny angles and output precision, compare vectorized geometry with
scalar matrices, verify block-size-independent output, and diagnose invalid
input with its actual file line number. Test
temporary files remain inside build/mc-rotate-tests. See
[the results and scope](../docs/validation/mc_rotate.md). They do not require
a C++ build or PyYAML and do not replace propagation validation.

## Worker lifetime regression

`unit.concurrent.PropagationLifetime.WorkerCleanupSerializesSharedReferenceRelease`
uses two and four OpenMP workers over repeated empty-stack batches. A cloned
interaction retains a shared intrusive token and observes overlapping clone
destructors. The test checks serial cleanup, restored reference counts after
each batch, and exactly one token destruction after the final owner releases
it. The empty stack isolates lifetime handling from simulation physics/RNG.

Run `ctest --test-dir build/Release -R unit.concurrent. --output-on-failure`.
The automatic [Linux build/test workflow](../docs/usage/ci.md) runs it as part
of the complete short suite on every PR.

The `unit.thinning.*` suite checks optional photon-threshold thinning,
including fixed Bernoulli moment controls, RNG bypass and fatal numerical
errors. Its contract and frozen short-check settings are documented in
[photon thinning](../docs/usage/photon_thinning.md). The configuration suite
checks CLI/YAML/snapshot selection and preservation of the energy policy.
