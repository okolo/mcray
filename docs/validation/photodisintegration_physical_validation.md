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
