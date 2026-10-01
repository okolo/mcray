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
Agreed: add explicit A,Z to Particle while retaining ParticleType compatibility;
select one exclusive TALYS channel per event.
Legacy particles remain specified through ParticleType; explicit isotopes are
set through a method that keeps A,Z,Type consistent. Use tabulated CRPropa
masses for TALYS and E_i=E_parent*A_i/A_parent.
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
Field details and consistency checks, out-of-grid behavior, and class names
have not yet been agreed.

Use the predecessor added by the user, `../transportcr`, for PSB validation.
It solves one-dimensional transport equations; baseline commit:
`00d5dd5c3b69f00fdc51330f6f81fa4b9af71791`.
Validation is described in transportcr_psb_validation.md.
TransportCR is useful for checking the port, but its PSB formulas are not
independent of mcray: check shared defects separately.
At the user's explicit instruction, the getRate range check was fixed in
both codes (`||` replaced with `&&`).
Do not commit in TransportCR; obtain separate agreement for other changes there.

## Working baseline and locations

- Development repository: the current `mcray` directory.
- Working branch: `talys`, created from `main` at commit
  `f62e40fcb1a62fd26377a2c416b6b9332499e8a9`.
- Read-only reference: sibling `../CRPropa3` at commit
  `73fa865b3304113da2059755a83f5939c1039d68`.
- `spec.md`: mcray internals and limitations of the original model.
- `../crpropa_spec.md`: CRPropa specification, outside both Git repositories.
- `mcray_crpropa_comparison.md`: comparison and planning questions.
- `talys_implementation_plan.md`: interactive draft plan and decision log.
- `transportcr_psb_validation.md`: predecessor analysis and PSB checks.

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
C++14, and Homebrew dependency paths. PATH ordering prevents the Miniforge
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
CTest includes table_function_log_zeros, psb_branching_bounds,
crbeam_help, and z2t_smoke. They do not replace physical PSB/TALYS validation.
Run these checks on the primary platform for subsequent changes;
do not automatically include long historical calculations in the short suite.

Verification dated 2026-10-01: the full Release build and all four CTest tests
passed in Docker Linux amd64 (uname: x86_64), emulated on macOS arm64.
Toolchain: GCC/G++/gfortran 13.3.0, Boost 1.83, GSL 2.7.1, Ubuntu 24.04.
All four tests also passed natively on macOS arm64 with Clang 22.1.5.
A native run on physical Linux amd64 hardware has not yet been performed.
