# Photodisintegration physical-validation diagnostics

Date: 2026-10-11. Source: mcray 18f230f985c941d2acbea452d98be76ec1264d81
plus validation-adapter changes recorded by hashes. Production physics is
unchanged. This is the first network numerical gate, not completed physical
acceptance or independent validation of nuclear cross-section measurements.

## Scope and references

The combined Fe56 TALYS+decay network contains 158 bound isotopes and free
p/n. Each bound isotope is checked at gamma=1e9, 3.9810717055349693e9 and
1e10 (474 isotope/gamma cells), including the 13 approved PD terminals.
The scan records full rates, exclusive sums and every exclusive channel.
A failure invalidates that cell; it is retained explicitly and is not counted
as a successful zero or dropped from the future acceptance family.

Use CRPropa3 73fa865b3304113da2059755a83f5939c1039d68 and CRPropa3-data
34b10e27375c0e9894d7c984ac15fb96e12975ef. The verified photon-enabled
Fe56 bundle retains the previously approved raw nuclear cross sections and
corrected decays. Photons, GZK and PPP are not sampled by these rate checks.
Identical frozen z=0 density tables are used on CMB and CMB+Inoue12.
The CMB table is the existing 4097-node export; the combined field is the
existing matched export, not a substitute native CRPropa background.

Two distinct numerical references are retained:

- The pinned CRPropa convolution on its original 513-node rest-energy grid,
  plus its serialized 201-node gamma-table rate interpolation. The latter
  does not emulate nearest-node runtime branching.
- Independently implemented analytic integration of the agreed piecewise
  photon-density interpolant, followed by Gauss-Legendre rest-energy
  integration split at every cross-section node. Orders 8, 16 and 32 check
  observed convergence. These integrate the approved cross-section interpolant;
  they do not assess its physical accuracy or unknown high-energy tail.

## Native network results

Each setting below scans the same 474 cells. The logarithmic background step
and integration tolerance are changed only in the diagnostic adapter.
Relative-error maxima apply to successful rates above 0.1/Mpc and are
measured against the independent order-32 convolution. This resolution
threshold is diagnostic, not an approved absolute acceptance budget.

| Background | log10 step | Relative tolerance | Failed cells | Max full-rate difference | Max exclusive-rate difference |
| --- | ---: | ---: | ---: | ---: | ---: |
| CMB | 0.01 | 1e-4 | 1 | 0.04081% | 0.12057% |
| CMB | 0.005 | 1e-4 | 2 | 0.03022% | 0.08041% |
| CMB | 0.005 | 1e-5 | 217 | 0.00937% | 0.00935% |
| CMB+Inoue12 | 0.01 | 1e-4 | 7 | 0.02255% | 0.10325% |
| CMB+Inoue12 | 0.005 | 1e-4 | 5 | 0.00598% | 0.16757% |
| CMB+Inoue12 | 0.005 | 1e-5 | 387 | 0.00566% | 0.00583% |

The tighter-tolerance rows have different successful subsets. Their smaller
reported errors cannot demonstrate network convergence: most cells fail.
At default controls, CMB fails for Li7 at gamma=1e10. CMB+Inoue12 fails
for Li7 (two cells), O16, Ne20 (two cells), V50 and V51. Failures report
GSL_EROUND: roundoff prevents achieving the requested tolerance. No failed
integral is accepted and no GSL status is suppressed. Halving spacing or
merely tightening tolerance does not establish robust numerical convergence.

Default CMB yields 473 successful full rates and 5306 exclusive rows;
245 full and 658 exclusive rates exceed the diagnostic resolution threshold.
The largest resolved full-rate residual is N13 at gamma=1e10; the largest
exclusive residual is Mg25 channel 000001 at gamma=1e10. Default CMB+Inoue12
yields 467 full and 5226 exclusive rows, of which 253 and 657 are resolved.

On the successful default CMB subset, the original pinned CRPropa quadrature
has maximum differences of 1.34346% (full, Fe56) and 1.52977% (exclusive,
O16 channel 000002). For CMB+Inoue12 the corresponding maxima are 1.12187%
and 1.23988%. These are direct-quadrature residuals, before additional
runtime gamma interpolation or branching effects. They reinforce the need
to retain both official and independently converged reference controls.

For successful resolved default rows, the maximum relative order-16 versus
order-32 change is below 3.1e-11 on both backgrounds. This is an observed
rest-quadrature refinement result, not a rigorous complete error bound.
The analytic photon integral passes 36 comparisons with adaptive SciPy
quadrature across positive and zero-ended synthetic density intervals;
a steep-tail test checks finite nonnegative monotonic integrals.

## Verification, provenance and next gate

The changed native probe builds successfully. All 21 selected TALYS, decay,
shared-background and no-interaction CTest entries pass (0.31 s). The two
Python diagnostic tests pass, including the 36 quadrature comparisons.
The existing full 125-test suite was not repeated for this adapter-only work.
No new Linux run or statistical transport acceptance is claimed here.

Commands, source/data hashes, six complete job logs, all successful rates,
failed-cell records, independent/reference comparisons and dependency versions
are retained outside Git in ../validation/pd-physical-validation-20261011.
The earlier B12-root exploratory scan is retained separately and is not the
Fe56 network result. Reproduction commands and failure exit code 3 are in
[the adapter instructions](../../tests/crpropa/README.md#network-rate-convergence-diagnostic).

The initial numerical gate remained open. The follow-up below addresses
these failures. The initial next step was to isolate the failing full/exclusive
integrands and resolve their integration stability without changing the
approved cross-section/background physics or ignoring error statuses. Then
repeat both backgrounds, extend gamma coverage and quantify absolute errors
for unresolved rates. Only after this gate should transport step convergence,
photon spectra and integrated acceptance ensembles proceed. Nuclear photon
yield saturation, missing-line physics and energy conservation remain their
separately documented limitations; these rate checks do not settle them.

## Integration fix and repeated network checks, 2026-10-11

The user authorized locating the failing channels, fixing numerical integration
without changing the physical model, and repeating both backgrounds.
Component-by-component reruns of the original unsplit quadrature identify:

| Parent | Exclusive code | Emitted fragments | Failed gamma on default shared backgrounds |
| --- | --- | --- | --- |
| Li7 | 100000 | n | CMB: 1e10; combined: 1e9, 1e10 |
| O16 | 000001 | He4 | Combined: 1e10 |
| Ne20 | 001000 | D | Combined: 1e9, 3.9810717055349693e9 |
| V50 | 110000 | n+p | Combined: 3.9810717055349693e9 |
| V51 | 400000 | 4n | Combined: 1e10 |

The full cross sections succeed in these default failing cells; their exclusive
channels cause the failures. The component lists at refined controls are also
retained. A component failure does not establish a divergent physical rate.

ContinuousBackgroundIntegral::GetRateK now splits tabulated cross-section
integrals at the cross-section interpolation knots and both kinematic mappings
of the photon-integral knots. The second mapping is included when beta<1.
Each piece uses the same 61-point adaptive GSL QAG rule and requested relative
tolerance; the absolute allowance, when supplied, is distributed by interval
width. Compensated summation combines the pieces. GSL errors remain fatal to
the rate calculation. Non-tabulated cross sections retain the previous path.
The table accessor returns physical knots, including on scaled table clones.

The cross sections, zero-ended interpolation, photon-integral interpolation,
background construction, masses, kinematics and process-selection model are
unchanged. Splitting exposes known derivative breaks to quadrature rather
than requiring one global adaptive integration to discover them. The evidence
supports this as the source of the observed quadrature instability; it is not
an assertion that floating-point roundoff has been eliminated everywhere.

All six native repeat jobs succeed: all 474 isotope/gamma cells, 5310 exclusive
rows and 474 exclusive sums per job, with zero failed cells. No missing-cell
selection is used in the new comparisons. For rates above 0.1/Mpc:

| Background | log10 step | Relative tolerance | Max full-rate difference | Max exclusive-rate difference |
| --- | ---: | ---: | ---: | ---: |
| CMB | 0.01 | 1e-4 | 0.03766% | 0.03800% |
| CMB | 0.005 | 1e-4 | 0.00942% | 0.00951% |
| CMB | 0.005 | 1e-5 | 0.00942% | 0.00951% |
| CMB+Inoue12 | 0.01 | 1e-4 | 0.02234% | 0.02372% |
| CMB+Inoue12 | 0.005 | 1e-4 | 0.00558% | 0.00593% |
| CMB+Inoue12 | 0.005 | 1e-5 | 0.00559% | 0.00593% |

Halving the background step reduces the resolved residuals approximately
fourfold on both fields. At fixed step 0.005, tightening tolerance from 1e-4
to 1e-5 changes any resolved full/exclusive/sum rate by at most 1.50e-8
relative on CMB and 5.90e-8 on CMB+Inoue12. These observations separate the
remaining background-interpolation error from the former adaptive-quadrature
failures. The production background grid and default tolerance were not changed.

Two regression cases were added: a narrow zero-padded triangular tabulated
channel with an analytic nonzero rate, and physical interpolation knots after
scaled-table cloning. The former guards against a global quadrature missing a
narrow supported channel and checks an intermediate redshift and a background
clone. All 127 CTest entries pass natively (29.36 s) and in Docker Linux amd64
emulation (98.47 s). Two native 20-history Fe56 controls, PD photons and
PD+decay photons, complete and each produces 19 photon rows. These are wiring
controls, not statistical spectral acceptance.

The numerical fix costs additional quadrature work. The six native network
repeat scans took approximately 27–53 seconds each, whereas the original
ordinary scans took a few seconds. Concurrent test/container work affects
these wall times; they are not a controlled performance benchmark. No tables
were replaced by an accelerated approximation to recover speed.

All localization, repeat and reference-comparison artifacts are outside Git
in ../validation/pd-integration-fix-20261011. The original failed outputs are
preserved as before-*; after-* contains complete repeat scans. Both use the
same fixed data and gamma values. The run manifest freezes the patch, added
adapters, inputs, commands, binaries, container identity and output hashes.

These checks resolve the observed failures in this sparse numerical matrix.
A full gamma scan, native-background effects, near-zero absolute budgets,
competing-process step convergence, photon-spectrum comparisons and the frozen
statistical acceptance campaign remain separate required validation work.

Both default full-network repeats also succeed in Docker Linux amd64 emulation:
474 isotope/gamma cells, 5310 exclusive rows and 474 sums on each background,
with no failed cells or zero/nonzero disagreements against native macOS.
Maximum cross-platform relative differences across all nonzero rows are
6.67e-16 (CMB) and 8.89e-16 (CMB+Inoue12). The two Linux jobs took 92.23 s
and 112.71 s. This reproduces the sparse matrix; it is not a native-amd64
hardware performance benchmark or a complete energy-domain acceptance claim.


## Integration optimization, 2026-10-11

The optimized implementation retains all cross-section and photon-integral
breakpoints, compensated summation, requested tolerances and fatal GSL error
handling. Each smooth piece now uses adaptive QAG15 instead of QAG61. It also
skips intervals certified exactly zero by the explicit zero-ended interpolant,
and avoids evaluating the neighboring redshift plane when its interpolation
weight is exactly zero. No positive contribution is discarded by a threshold;
non-tabulated cross sections retain the whole-range QAG61 path. The input
tables, interpolants, physical approximations and default grids are unchanged.

All six full-network rate scans above succeed again, with 6258 rows each and
no failed cells or zero/nonzero disagreements against stable split-QAG61.
The maximum relative rate difference is 9.12e-14 on CMB and 2.11e-15 on
CMB+Inoue12, across full rates, exclusive rates and sums, including weak rows.
Comparison to the hash-verified independent convolutions retains the residuals
and background-grid convergence reported above. This preserves the achieved
accuracy in the tested matrix; it does not establish complete energy-domain
acceptance.

The regression suite now has 129 passing entries: 29.15 seconds natively and
99.42 seconds in Docker Linux amd64 emulation. The analytic narrow-channel
case checks exact redshift nodes, an intermediate redshift and clones with a
redshift-dependent background; callback counts verify the skipped plane.
A separate conservative-zero test guards against discarding a narrow positive
peak. CLI tests cover the optional explicit unsigned CRbeam seed, including
zero; omitted seeds preserve automatic seeding. The effective seed is recorded
in production output parameters.

Both default full-network rate scans also succeed in the optimized Linux
container. All 6258 rows match their native keys and zero/nonzero status;
maximum cross-platform relative differences are 1.45e-15 (CMB) and 1.00e-15
(CMB+Inoue12). These checks use emulation, not native amd64 hardware.

### Controlled complete-propagation timings

The baseline is the stable split-QAG61 implementation described above, with
identical optional CRbeam seed support in both binaries. It is not the original
unsplit implementation that failed some rate cells. Both variants were frozen
using the same Release configuration and toolchain before timing, on macOS arm64 with
OMP_NUM_THREADS=1. No builds or other calculation jobs ran concurrently.

Each case uses 64 monoenergetic Fe56 primaries at E/A=1e19 eV, TALYS and
nuclear decays, and both nuclear photon sources. Each executable receives an
eight-primary warmup, followed by three paired runs with seeds 2026101100,
2026101101 and 2026101102. Pair order alternates. Times include initialization,
complete propagation of tracked secondaries and output. Nuclear-only scenarios
use L=1 Mpc, retaining photons without EM interactions. Production CRbeam
uses z=2.5e-5 (approximately 0.110535 Mpc), EM cascades on, a 1e11 eV tracking
cutoff and thinning off. Nuclear PPP/GZK are off in both variants; production
free-nucleon and EM processes remain enabled. These background cases use native
mcray backgrounds, separately from the shared-density rate checks.

| Complete calculation | Background | Before median, s | After median, s | Speedup |
| --- | --- | ---: | ---: | ---: |
| Nuclear histories, no EM interactions | CMB | 84.40 | 10.00 | 8.44x |
| Nuclear histories, no EM interactions | CMB+Inoue12 | 96.18 | 10.24 | 9.39x |
| CRbeam with EM cascade | CMB | 16.40 | 10.73 | 1.53x |
| CRbeam with EM cascade | CMB+Inoue12 | 26.91 | 18.87 | 1.43x |

The ratios compare medians; different seeds produce different cascade workloads.
All 24 measured jobs and eight warmups complete successfully within the local
wall limit. All 12 measured pairs have identical species, counts, weights and
per-species summed output energies. All printed energy/weight rows match
pairwise. Nuclear-only non-comment particle rows are identical; CRbeam retains identical interaction counts, with
small differences in printed trajectory coordinates from floating-point
propagation (maximum absolute difference 1.01e-13 in dimensionless columns;
printed production redshifts match). Headers can differ in the last digits
of their rate diagnostic.
This is a paired regression and performance result, not a statistical
physical-acceptance claim. The larger nuclear-only speedup does not transfer directly to cascades, whose EM work is
largely unaffected. Other energies, path lengths, cutoffs and process mixtures
may have different speedups.

The reproducible driver is tests/benchmark_nuclear_integration.py, documented
in tests/README.md. Frozen before/after executable and source snapshots,
patches, rate/reference comparisons, commands, seeds, output checks and hashes
are outside Git in ../validation/pd-integration-optimization-20261011.
The independent reference files were checked against the previous frozen
manifest before reuse; residuals were recomputed from the optimized rates.


## Full common gamma-domain grid, 2026-10-11

Source: mcray e43a8e3b623fc2e16a865a0e82013d5f4341617b; the native probe
matches the previously frozen optimized binary by SHA-256. Only validation
adapters/documentation changed in this stage. Nuclear data and both shared
z=0 photon-density tables retain the pinned versions and checksums above.

The new grid spans gamma=1e6 through 1e14 in 0.02 log10 increments: 201
original CRPropa gamma nodes plus 200 logarithmic midpoints. This covers the
entire common tabulated gamma domain on a finite grid; it is not a bound
between every pair of sampled points. All 158 bound isotopes in the combined
Fe56 network are included, with the 13 approved PD terminals remaining zero
under their existing data contract. Decays define network closure; decay,
photon, PPP and GZK events are not sampled by these rate-only calculations.

All six native settings complete without GSL errors, negative/nonfinite rates,
failed cells, duplicate keys or missing/unexpected rows. Each setting contains
63358 isotope/gamma cells, 63358 full rates, 63358 exclusive sums and 709770
exclusive rates: 836486 rows, or 5018916 across the six settings. Coverage is
checked against the prepared network and channel identities, not only row
counts. The independent/reference comparisons also retain every requested row.

Both default matrices are compared to unchanged pinned official convolution,
its serialized 201-node gamma interpolation, and independent interval orders
8/16/32. Cached photon samples only remove repeated identical evaluations;
the nine-gamma pilot has exactly identical 17352 comparison rows with caching
on and off. Maximum order-16 versus order-32 changes for resolved rates are
below 6.9e-11 relative across both backgrounds, including exclusive sums.
This remains observed refinement evidence, not a rigorous error bound.

The resolution threshold remains the diagnostic independent rate >0.1/Mpc.
No new absolute acceptance budget is inferred from it. Maxima relative to the
independent order-32 reference are:

| Background | log10 step | Relative tolerance | Max full error | Max exclusive error | Max sum error | Max weak absolute error, 1/Mpc |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| CMB | .01 | 1e-4 | 0.04211% | 0.04097% | 0.04206% | 4.024582e-05 |
| CMB | .005 | 1e-4 | 0.01054% | 0.01025% | 0.01054% | 1.006312e-05 |
| CMB | .005 | 1e-5 | 0.01054% | 0.01025% | 0.01054% | 1.006136e-05 |
| CMB+Inoue12 | .01 | 1e-4 | 0.03445% | 0.03379% | 0.03514% | 3.507335e-05 |
| CMB+Inoue12 | .005 | 1e-4 | 0.00861% | 0.00844% | 0.00878% | 8.772369e-06 |
| CMB+Inoue12 | .005 | 1e-5 | 0.00861% | 0.00844% | 0.00878% | 8.769725e-06 |

Halving background spacing again reduces both resolved and weak absolute
residuals approximately fourfold. At spacing 0.005, tightening tolerance
from 1e-4 to 1e-5 changes any resolved full/exclusive/sum rate by at most
1.951e-08 (CMB) relative; 6.572e-08 (CMB+Inoue12) relative.

The largest default full-rate residuals occur for Cr52 at gamma=3.4673685e9
(CMB) and gamma=3.0199517e9 (combined field). This is consistent with the
remaining background-integral interpolation error rather than unresolved
adaptive-quadrature failures. Production defaults have not been changed.

Weak rates are retained in the full outputs and absolute-error summaries.
In the extreme CMB tail, mcray can return zero while the independent reference
is still positive; the largest such exclusive rate is 3.735e-289/Mpc
(full/sum: 2.879e-289/Mpc). Relative discrepancies can therefore reach 100%
without controlling the useful absolute scale. The combined field has no
such zero/nonzero mismatch in this matrix. Neither field has a positive
mcray rate with an exactly zero independent reference. These statements do
not substitute for an approved near-zero acceptance budget.

### Official-control differences at high gamma

The wider scan exposes a finite-upper-limit term in the pinned official
`calc_rate_eps` convolution, in addition to its rest-grid quadrature error.
For the agreed cross section, zero outside its grid, define

\[
F(K)=\int_{K_0}^{K}t\sigma(t)\,dt,\qquad
I(x)=\int_x^{\epsilon_{\max}}\frac{n(\epsilon)}{\epsilon^2}\,d\epsilon.
\]

The ultrarelativistic integration-by-parts identity is

\[
R=\frac{F(K_{\max})I(K_{\max}/2\gamma)}{2\gamma^2}
 +\frac{1}{\gamma}\int_{K_0}^{K_{\max}}
 \frac{F(K)n(K/2\gamma)}{K^2}\,dK.
\]

The pinned official formula computes the second term on its finite grid.
The first term need not vanish: sigma(Kmax)=0 does not imply F(Kmax)=0.
At large gamma, photons above Kmax/(2 gamma) can make this contribution
important even with the same zero-truncated cross section. This effect is
separate from uncertainty in the physical cross section beyond its grid.

For Be11 on CMB at gamma=1.2589254117941663e12, the independent rate is
0.1039611838865/Mpc and the missing endpoint term is 0.0082774704731/Mpc
(7.9621%). The unchanged official result is 0.0956954594435/Mpc, a -7.9508%
residual. On all resolved rows, maximum direct official discrepancies are
7.9508% on CMB and 7.9527% on the combined field. The serialized gamma-table
control has discrepancies up to approximately 7.93%; interpolation is a
separate effect, with an additional increment of up to approximately 5.30%
of the independent rate on resolved CMB rows near a steep onset.

The separate `diagnose_rate_boundary.py` control uses Be11, Be9 and Fe56 on
both backgrounds. It refines the agreed rest-grid interpolant by factors
8/16/32/64 and adds the analytic endpoint contribution only diagnostically.
For resolved control rows, corrected factor-64 rates agree with independent
interval quadrature within 1.06e-6 relative, with observed refinement
convergence. Official source/tables and production rates remain unchanged.
No nearest-node runtime channel sampling is validated by the gamma-table
rate interpolation control.

### Reproduction and remaining scope

`tests/crpropa/run_network_rate_scan.py` reproduces the full matrix, validates
coverage, runs both reference controls and produces separate weak/resolved
summaries. Seven adapter tests and five server-runner tests pass on macOS.
A complete nine-gamma driver pilot, including tests, checksum verification,
both backgrounds, reference comparisons, summaries and endpoint controls,
completes successfully in 19.77 seconds. The unchanged C++ regression suite
was not rerun for this adapter-only stage; its prior native/Linux result is
129 passing tests. No new full-domain Linux result is claimed here.

All native matrix jobs were bounded; none was expected to exceed the local
30-minute limit. A separate native Linux amd64 extended workflow is prepared:
801 gamma points, additional spacing 0.0025/tolerance 1e-6, and independent
orders through 128. It fetches the requested Git revision at startup, records
the exact commit, builds/tests, verifies transferred inputs/tools, executes
adapter tests and archives successes or failures into the configured test-data
directory. The extended calculation has not been executed. Commands and dependencies are documented
in [the adapter instructions](../../tests/crpropa/README.md#full-gamma-domain-network-scan).

The first external launch on 2026-10-12, run ID
`20261012T044047Z-8b742449`, verified the transferred payload and fetched the
requested source commit, but failed during build because the server CMake
did not support `--fresh`. CTest and the rate workload did not start. The
runner reported a failed archive in the configured test-data directory.
The adapter instructions now explain installing CMake >=3.24 in the server
virtual environment. No successful retry or Linux numerical result is yet
available.

Artifacts, commands, six coverage checks, complete rates/reference rows,
absolute-error summaries, endpoint diagnostics and source/input hashes are
outside Git in ../validation/pd-full-gamma-20261011. The final transferable
bundle is external-inputs-v2.tar.gz with a separate SHA-256 file; the earlier
bundle and pilot source snapshots are preserved. There are no random seeds
in deterministic rate scans.

This completes the finite-grid native numerical rate check across the common
gamma domain. Extended Linux and gamma/grid refinement, native-background and
redshift checks, near-zero acceptance budgets, transport-step convergence,
photon spectra and frozen statistical acceptance remain open. The unknown
physical high-energy tail and nuclear-photon energy/coverage limitations are
unchanged.
