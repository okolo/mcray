# mcray project task: nuclei and TALYS

## Language policy

Everything committed to this repository must be in English. This includes
documentation, Markdown files, source-code comments, identifiers and
human-readable text introduced in code, scripts, configuration descriptions,
test messages, and commit messages. Do not add non-English explanatory text
to files intended for commit. Preserve scientific data, equations, file formats,
and externally defined identifiers when translating existing text.
Conversation with the user may remain in Russian.

## Commit policy

Use Conventional Commits for every commit:
`<type>[optional scope][!]: <description>`.
Write the subject, body, and footer text in English. Source-code comments
must also be in English, as required by the language policy above.
Use an appropriate type such as `feat`, `fix`, `docs`, `build`, `test`,
`refactor`, `perf`, `ci`, `chore`, or `revert`.
Keep the description concise and imperative, for example:
`fix(psb): validate the multinucleon branching index`.
Mark breaking changes with `!` before the colon or a
`BREAKING CHANGE:` footer explaining the change in English.
Keep commits focused and exclude unrelated changes and generated artifacts.

## Objective

Implement the propagation of atomic nuclei through intergalactic space in
mcray, with photodisintegration based on TALYS.
First reproduce the selected physical assumptions and verifiable results of
CRPropa; consider improvements after validating the baseline implementation.
TALYS provides photonuclear cross sections, not the complete transport model.
The full task includes isotopes, residual nuclei and secondary particles,
photon backgrounds, cosmology, decays, and associated energy losses.

User clarification dated 2026-10-01: retain the PSB approach for comparison
with TALYS. The existing implementation is insufficiently validated and may
contain errors; do not treat it as a numerical reference before validation.
Preserving the physical approach does not require preserving implementation bugs.
The user selected validated PSB with fixes and the original version preserved
in history. Do not plan a separate compatibility mode for the old behavior.
The original version is recorded at commit
`f62e40fcb1a62fd26377a2c416b6b9332499e8a9`; document changes in results
caused by fixes.
The agreed first stage covers PSB validation and fixes, and TALYS
photodisintegration. Add decays, nuclear pair production, and photopion
processes later; this intermediate result is not a complete nuclear transport model.
Primary nuclei in the first version are limited to Fe-56, with the required
daughter isotopes. Initially use raw CRPropa3-data cross sections and
convolution with mcray backgrounds.
Current agreed direction (D27, 2026-10-08): one Nucleus ParticleType with
explicit A,Z, superseding the old per-A enum compatibility requirement.
Free p,n retain their types. Set isotope identity through a method that keeps
A,Z,Type consistent. Select one exclusive TALYS channel per event. Use
tabulated CRPropa masses for TALYS nuclei and E_i=E_parent*A_i/A_parent.
PSB retains weighted branching, its chain indexed by A, and approximate masses;
fix implementation errors. PSB and TALYS are separate processes sharing
particle helpers and background integrals.
For TALYS, obtain the total rate from the full cross section and channel
probabilities P_c=R_c/sum(R_c) from convolved exclusive cross sections.
Terminate the parent particle; create the residual and emitted particles as
products with the parent's weight, without a probability multiplier for an
already selected channel.
Take light nuclei from PD_external in the selected CRPropa3-data version.
Treat missing required masses/cross sections and inconsistent channels as
calculation errors. Validate the reachable network before a run and diagnose
runtime violations; do not use approximate substitutions or silent omissions.
Initially compute rates through direct BackgroundIntegral convolution;
consider rate tabulation for acceleration after validation.
Interpolate log(sigma) against log(K) where both nodes are positive.
On intervals with a zero node, interpolate sigma against log(K), without
epsilon substitutions. Two zero nodes imply sigma=0.
Reuse the existing interpolation implementation: LinearFunc with LogScale on
both axes supports the positive intervals.
At the user's instruction, no new class was created: LinearFunc received an
optional aHandleLogZeros=false parameter. Setting it to true enables the
agreed zero-interval handling; Clone preserves the parameter. Legacy calls
retain their behavior. The global semantics of LogScale were not changed.
The unified identity/mass contract, pinned data and zero-padded grid policy
were approved on 2026-10-08; see D29-D31 and docs/usage/talys_usage.md.

Use the predecessor added by the user, `../transportcr`, for PSB validation.
It solves one-dimensional transport equations; baseline commit:
`00d5dd5c3b69f00fdc51330f6f81fa4b9af71791`.
Validation is described in https://github.com/okolo/mcray_dev/blob/main/docs/reference/transportcr_psb_validation.md.
TransportCR is useful for checking the port, but its PSB formulas are not
independent of mcray: check shared defects separately.
At the user's explicit instruction, the getRate range check was fixed in
both codes (`||` replaced with `&&`).
The user superseded the earlier TransportCR commit prohibition on 2026-10-08:
its existing fixes were committed by the user on master, and future authorized
work uses the mcray-validation branch of the transportcr workspace submodule.
Read its AGENTS.md before changes; do not modify physics without agreement.

## Working baseline and locations

- Development repository: the current `mcray` directory.
- Working branch: `talys`, created from `main` at commit
  `f62e40fcb1a62fd26377a2c416b6b9332499e8a9`.
- Read-only reference: sibling `../CRPropa3` at commit
  `73fa865b3304113da2059755a83f5939c1039d68`.
- `https://github.com/okolo/mcray_dev/blob/main/docs/reference/spec.md`: mcray internals and limitations of the original model.
- `https://github.com/okolo/mcray_dev/blob/main/docs/reference/crpropa_spec.md`: CRPropa specification, outside both Git repositories.
- `https://github.com/okolo/mcray_dev/blob/main/docs/design/mcray_crpropa_comparison.md`: comparison and planning questions.
- `https://github.com/okolo/mcray_dev/blob/main/docs/design/talys_implementation_plan.md`: interactive draft plan and decision log.
- `https://github.com/okolo/mcray_dev/blob/main/docs/reference/transportcr_psb_validation.md`: predecessor analysis and PSB checks.

## Workflow

1. Study both codebases and prepare specifications.
2. Compare architecture, physical models, and data.
3. Develop the implementation plan from the comparison.
4. Implement and validate baseline TALYS support.
5. Only then investigate improvements to the physical model and numerical methods.

Steps 1–2 are covered by documents dated 2026-10-01. An interactive plan
draft has been started; the PSB getRate range check has been fixed.
At the user's separate instruction, the LinearFunc zero-handling parameter
was implemented and checked with a standalone C++ test.
Before continuing, read the documents and decision log.
Develop the plan with the user and ask clarifying questions.
Do not make architecture decisions without approval: record options and
recommendations as proposals, and user answers as agreed decisions.
Do not interpret silence as agreement or begin dependent changes.

## Implementation rules

- Existing untracked files, IDE directories, builds, and calculation results
  belong to the user; do not delete them or include them in changes.
- Study CRPropa as a reference; implement changes in mcray.
- Tie conclusions to code and data versions; the external CRPropa3-data
  `master` branch is not a pinned table set.
- First define the contract for identity `(A,Z)`, masses, charges,
  secondaries, weights, and output. The current `A02`…`A56` enum does not
  distinguish isotopes.
- Audit enum arithmetic and `ParticleTypeEOF` arrays.
  Do not add nontrivial Particle fields while retaining memcpy/memset.
- Separate raw TALYS cross sections from rates already convolved with a
  particular background. Do not apply rates for one EBL to a different
  EBL or an arbitrary background.
- Explicitly distinguish mcray natural units from CRPropa SI units,
  rates per time, physical and comoving lengths, and dn/dE versus E dn/dE.
  Account for c and (1+z) factors exactly once.
- Define behavior for missing data, out-of-grid values, and light nuclei;
  a missing table does not imply physical stability.
- Keep the total rate consistent with exclusive-channel probabilities.
- Do not run TALYS for each transport event: select a reproducible data
  preparation workflow and immutable runtime tables separately.
- Account for Clone(), independent RNGs, and OpenMP when implementing processes.

## Baseline validation

First compare isolated rates on identical backgrounds at z=0, then
cosmological scaling, single channels and isotope chains, and finally
post-propagation spectra and composition.
Check A and Z in exclusive photonuclear events, energy balance under the
selected kinematics, channel normalization, absence of negative/NaN rates,
and correct weights. For beta decays, check charge together with the leptons,
not only among nuclei.
Perform the first comparison without thinning, using identical data,
cosmology, and step controls. Equal seeds in different codes do not imply
equal trajectories; compare distributions with statistical uncertainties.

macOS build: `cd bin && bash ./build.sh Release`.
Output: `build/Release/CRbeam` and utilities in `build/Release/mcray/`.
The old `--test nucl_dis` does not replace new TALYS checks.

### Build verification through build.sh

`bin/build.sh` adapts the TransportCR script: the CMake project is
`src/app/crbeam`, with output in `build/Release` or `build/Debug`.
Both CMakeLists.txt files restrict the `--as-needed` linker flag to Linux.
On macOS the script selects Homebrew Clang for C/C++, gfortran for Fortran,
C++17, and Homebrew dependency paths. PATH ordering prevents the Miniforge
assembler from shadowing system tools. Linux retains GCC selection.

The former crash in Utils::default_tables_path() occurred when the GCC
application used libstdc++ while Homebrew Boost.Filesystem used libc++.
Switching to Clang removed the crash; rebuilding Boost was unnecessary.
The earlier GCC build also failed with the macOS 26 SDK Mach headers;
the current Clang build does not require an older SDK.

Verified on macOS arm64: Homebrew Clang 22.1.5 and gfortran 15.2.0.
`bash bin/build.sh Release` from the repository root built CRbeam,
mcray_test, and z2t without extra environment variables.
`build/Release/CRbeam --help` returned 0; OpenMP reported 10 available threads.
otool showed that the application and Boost use libc++; neither libstdc++
nor libgomp is linked, and LLVM libomp is used.
`build/Release/mcray/z2t 0.1` completed successfully.
When launched from bin, tables are found through the `tables/` fallback;
the warning about missing tables beside the executable remains.
Logs: `/tmp/mcray-clang-build.log` and `/tmp/mcray-clang-help-bin.log`.
Full nuclear transport calculations were not run; Debug was not verified.
Primary-platform Linux verification is described below.

### Primary platform: Linux amd64

The user specified Linux amd64 as the primary target platform.
The Ubuntu 24.04 Docker environment is defined in
`docker/linux-amd64.Dockerfile`.
Run `bash bin/test-docker.sh` from the root; see docker/README.md.
The build and run explicitly request linux/amd64.
On macOS arm64 this uses emulation rather than a native amd64 CPU.
CTest discovers GoogleTest unit cases for interpolation and PSB branching,
and includes crbeam_help and z2t_smoke. They do not replace physical PSB/TALYS validation.
Run these checks on the primary platform for subsequent changes;
do not automatically include long historical calculations in the short suite.

Before the GoogleTest migration, verification dated 2026-10-01 showed that
the full Release build and all four CTest tests
passed in Docker Linux amd64 (uname: x86_64), emulated on macOS arm64.
Toolchain: GCC/G++/gfortran 13.3.0, Boost 1.83, GSL 2.7.1, Ubuntu 24.04.
All four tests also passed natively on macOS arm64 with Clang 22.1.5.
A native run on physical Linux amd64 hardware has not yet been performed.

## C++ and test framework policy

The user selected GoogleTest and migration to C++17 on 2026-10-01.
Require C++17 in the application and library builds. Use GoogleTest for
new C++ unit and integration tests, with CTest as the common runner.
The test-only dependency is GoogleTest 1.18.0, pinned by archive SHA-256
in tests/CMakeLists.txt; production builds may disable BUILD_TESTING.
See tests/README.md for local runs, labels, and offline dependency setup.
Existing command-line smoke checks run directly through CTest.
The approved nuclear validation plan uses GoogleTest for end-to-end output
assertions, with scripts preparing and running isolated reference cases.
Do not treat smoke checks as validation of nuclear propagation.
Use explicit numerical tolerances; compare Monte Carlo distributions with
statistical uncertainties rather than requiring identical trajectories.

C++17/GoogleTest verification dated 2026-10-01: full Release builds and all
26 CTest cases (24 unit, two smoke) passed on macOS arm64 and Docker
Linux amd64 under emulation. BUILD_TESTING=OFF configured successfully
without downloading GoogleTest. Physical propagation checks remain pending.

## Current validation task

User instruction dated 2026-10-01: plan unit tests and TransportCR E2E
comparisons for PSB, GZK, and pair production on nuclei. Use isolated
interactions, a point source at a nearby fixed distance (e.g. z=0.01),
redshift losses disabled, and the TransportCR short-distance regime with
XML input. Include mcray Monte Carlo uncertainty. See
https://github.com/okolo/mcray_dev/blob/main/docs/validation/nuclear_validation_plan.md for the proposal and implementation gaps.
This expands validation planning; it does not approve new nuclear GZK/pair
physical approximations. Existing mcray GZK supports nucleons and pair
production supports protons; nuclear implementation requires agreement.

The user accepted CMB-only monoenergetic Fe-56 as the first PSB comparison,
then CMB plus EBL, and light-nucleus PSB unit checks. Validate existing
nucleon GZK/proton pair code first; nuclear extensions require later approval.

The user approved https://github.com/okolo/mcray_dev/blob/main/docs/validation/nuclear_validation_plan.md and authorized execution.
Plan commit: 327bd46. Read docs/validation/psb_validation_results.md for current regression
findings, supported E2E scenarios and numerical/statistical limitations.
Do not mark PSB validated merely because short tests or one energy point pass.

The user selected Inoue12 Baseline for the CMB+EBL validation stage. Use the
identical EBL_Inoue/EBL_proper_baseline.dat table in both codes, record its
SHA-256, and freeze the background at z=0. The user supplied finer SOPHIA
reference tables at /Users/ok/propagation/tables/sophia2; copy selected files
into isolated comparison runs and record hashes without modifying the source
tables or committing in TransportCR.

For the initial CMB+EBL PSB comparison, use the short-path control L=1 Mpc
without thinning. The user capped the ensemble at 200,000 histories and
accepted marking the rare neutron yield inconclusive if its relative standard
error exceeds 1%; do not relax that criterion or claim full validation.
The user approved preserving different pair-product models: TransportCR's
E^(-7/4) approximation and mcray's differential QED kernel. Compare energy
balance and independently integrate the mcray kernel; do not force equal pair
counts or shapes by silently changing either physical model.


### Independent SOPHIA generation discovered on 2026-10-07

The predecessor generator is at /Users/ok/propagation/sources/sophia.
Use tests/transportcr/build_sophia_generator_docker.sh to copy and hash its
sources and build in Linux amd64 outside Git. Source provenance, pilot data,
and pending full-generation statistics are in docs/validation/sophia_table_generation.md.
Do not modify the original generator or replace its sampling with mcray.
Diagnostic 100-event tables are loader tests only, not acceptance references.
Separate table-generation Monte Carlo uncertainty from grid discretization
before claiming the fixed 0.2% reference criterion. Full table generation and
new acceptance ensembles remain subject to the recorded statistics decisions.


### PSB CMB+Inoue12 energy extension on 2026-10-07

Two additional Fe56 controls at E/A=1e18 and 10^18.5 eV, L=1 Mpc, have
complete 200000-history ensembles outside Git in
../validation/psb-inoue-extension-20261007. GetRateK now uses GAUSS61 at
unchanged tolerances after an A30 GSL_EROUND failure discarded histories.
All daughter rates are checked at working tolerance 1e-4; tolerance 1e-5
remains unsupported for A45 at the higher energy with spacing 0.01.
See docs/validation/psb_validation_results.md. Neutron relative SEs are 6.919% and 1.396%;
keep both outputs statistically inconclusive under the unchanged 1% target.
Do not increase the ensemble cap or reinterpret a precision-guard failure
as a physical disagreement without a separately agreed next step.


### Unified nuclear identity and CLI approval on 2026-10-08

The initial A-group Type proposal (D24) was superseded by D27 after the user
removed the legacy nuclear API compatibility requirement. Target one Nucleus
Type with explicit A,Z; migrate PSB and nuclear output/filtering together.
Free p,n retain Proton/Neutron and existing mcray masses. The user subsequently
approved D25/D26: store numeric nuclear mass with A,Z and preserve it on copy.

D28 approves unsigned CRbeam --A/--Z with defaults 0, A>=2 for nucleus,
Z=0 selecting the PSB-chain representative, and ignored physical A,Z for
non-nuclear primaries. First-stage primary scope remains Fe56. The CLI input
contract is implemented with strict unsigned parsing and named --particle
selection alongside legacy --primary. The unified migration now connects
separate PSB/TALYS processes through --photodisintegration TALYS/PSB/Off,
default Off. See docs/usage/talys_usage.md for data preparation and limitations.


### Approved TALYS data and finite coverage, 2026-10-08

Pin CRPropa3-data 34b10e27375c0e9894d7c984ac15fb96e12975ef. Use
PD_Talys1.8_Khan full totals and thin exclusive channels, PD_external for A<12,
and tabulated nuclear masses. Prepare the Fe56 network outside Git with
checksums using tools/prepare_talys_data.py. Extend logarithmic grids with
zero cross sections to 513 nodes as in CRPropa; sigma=0 outside that grid.
Keep the agreed LinearFunc interpolation and zero handling.

The user explicitly approved a photodisintegration-only terminal whitelist
for missing tables: (A,Z)=(5,2),(5,3),(6,2),(8,4),(9,5),(12,5),(16,9),(20,11),
(39,20),(42,21),(46,23),(47,24),(50,25). Warn and record it in metadata.
These states are not asserted to be physically stable. All other missing
required data remain errors. Nuclear decays and gamma emission are disabled
in this stage. Existing free-neutron decay and other CRbeam processes remain.

A separate future task, https://github.com/okolo/mcray_dev/blob/main/docs/tasks/T1/README.md, must assess the
unknown high-energy tail and the zero-padding transition. Do not change
production extrapolation or add high-energy mechanisms without approval.

Unified Nucleus/TALYS implementation verification: all 88 short CTest entries
pass on macOS arm64 and Docker Linux amd64 emulation. Real-data isolated
TALYS and CRbeam controls are documented in docs/usage/talys_usage.md. They do not
establish full CRPropa agreement, rare-output precision or validity of the
omitted high-energy tail. TransportCR was not modified by this migration.

### Nuclear decay stage approved on 2026-10-08

After the unified Nucleus/TALYS commit 448b48b, the user authorized nuclear
decay implementation. This supersedes the earlier first-stage exclusion
of decays, but does not approve nuclear GZK/pair production.
Use a separate NuclearDecay process and --nuclear-decay On|Off, default Off,
with tabulated TALYS masses; keep PSB unchanged and free-neutron decay separate.
The user approved CRPropa beta kinematics (same daughter gamma, no recoil),
energy sharing by A for nuclear emissions, compound channels, e+/e- and
neutrino products, with decay photons disabled. Exact event energy balance
is not guaranteed for beta decays under this approximation.
The user separately approved fixing the pinned generator's inclusive/exclusive
branch normalization and adding Be11 B-P as beta-minus plus proton emission.
Preserve the unmodified reference outside Git. Do not interpret missing decay
data as stability or enable artificial unknown-isotope dripping.
Read docs/usage/nuclear_decay.md for the data contract, preparation and validation limits.

Nuclear decay implementation verification: all 96 short CTest entries pass
on macOS arm64 and Docker Linux amd64 emulation. External Be8, B12 and
TALYS+decay Fe56 controls pass applicable GoogleTest assertions on both;
two-thread CRbeam smoke runs complete. Data and original/corrected generator
outputs remain outside Git. Full CRPropa spectral agreement and rare-branch
precision remain unestablished.

Nuclear decay and combined TALYS comparison planning is recorded in
https://github.com/okolo/mcray_dev/blob/main/docs/validation/nuclear_decay_comparison_plan.md. Proposed dual CRPropa data references,
statistical budgets and ensemble cap require the user's new approval;
previous GZK/PSB approvals do not automatically cover these comparisons.
Do not launch large acceptance ensembles before freezing the observable
families, reference data, numerical controls and validation seeds.

### Comparison execution authorized on 2026-10-08

The user authorized starting the comparison checks and required checks longer
than 30 minutes to be provided as a separate external Linux server script.
Use tests/crpropa/run_external_linux.sh for long ensembles; keep local runs
below 1800 seconds. The user also approved identical tabulated photon density
in comparison adapters, without changing production constants/backgrounds.
Measure native-background effects separately. Short decay pilots, deterministic
rates and TALYS numerical diagnostics are recorded in
docs/validation/nuclear_decay_comparison_results.md. Acceptance observables and numerical
reference treatment still need to be frozen before large TALYS ensembles.
The shared 51-node beta-CDF bias is documented; do not silently replace the
approved sampler. Production physics was unchanged during these diagnostics.

## Documentation organization

The user approved the hybrid scheme on 2026-10-08. This repository owns code,
public usage and implemented contracts, test adapters, and confirmed validation
results and limitations, indexed in [docs/README.md](docs/README.md).
Keep README.md and AGENTS.md at the repository root.

The private mcray_dev repository owns design, research, decision history,
validation plans and the future task index. Its remote is
https://github.com/okolo/mcray_dev. In the local adopted workspace, its root
is the parent of this repository and this checkout is its mcray submodule.
See ../AGENTS.md and ../docs/README.md when working in that private workspace.
Do not require private access for public builds, tests or usage. Each document
has one canonical home; link across repositories rather than maintain copies.
Future tasks have stable T1, T2, and subsequent IDs in the private index;
registration does not authorize implementation. T2 records the shared beta
sampling bias. Preserve this checkout's current branch and uncommitted work.
Do not publish a private gitlink until its mcray commit is available remotely.

Private design/research/planning documents belong only in mcray_dev. Before
publishing an unpublished development branch, verify that such files are absent
from every commit being pushed, not merely deleted from its final tree. Preserve
an ignored external backup before any authorized history rewrite. Never rewrite
published history or force-push without separate explicit authorization.
