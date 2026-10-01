# PSB validation progress

Date: 2026-10-01. Plan: nuclear_validation_plan.md, approved and committed
as 327bd46. Scientific PSB assumptions remain weighted branching, an A-indexed
chain, and approximate masses. Implementation changes below are not yet committed.

## Reproduced defects and repairs

1. **Uninitialized array index offset.** AutoDeletePtrArray::iMin had no
   initializer. PSB secondary-map construction indexed arrays before any
   setIndexShift call, throwing std::out_of_range in new production-map tests.
   Default offset is now zero; explicit shifted indexing is preserved.
2. **Unsupported channels in event sampling.** Above K=30 MeV, the total
   convolution excludes the Gaussian groups, but SampleSecondariesWithK
   checked only their lower thresholds. Gaussian::f itself is not truncated;
   consequently event selection included out-of-support Gaussian tails.
   A fixed-K Fe-56 test failed at event 8 with neutron yield 0 instead of
   the expected high-energy yield 1.669. Sampling now checks both bounds,
   excludes zero cumulative-probability entries, and rejects empty support.
3. **Clone changes the physical sampling mode.** A PSB object configured for
   S sampling cloned with the default K mode. A routing test observed rate
   2 in the original and 1 in its clone. Clone now retains fSampleK and uses
   an independent BackgroundIntegral clone, as required by its thread-safety
   contract. Non-nuclear GetSecondaries requests now return false.
4. **Quadrature failure loses MC histories.** With CMB, Fe-56 E/A=
   2.238721138568347e18 eV and integration tolerance 1e-4, sampling a daughter
   A55 reproducibly raised GSL_EROUND. PropagationEngine catches such failures
   and discards the unfinished primary cascade. The first 20,000-primary pilot
   retained only roughly 60% of histories. Its output is invalid and excluded
   from physical comparisons. A specific daughter regression reproduces the
   failure independently of the engine. The logarithmic sampler now uses the
   61-point rather than 15-point Gauss-Kronrod rule at the same tolerance.
   The regression and subsequent full-history checks pass; no tolerance was
   relaxed. The validation runner fails if any injected history is missing.

The propagation engine's general exception/retry API has not been redesigned.
Shared TransportCR getRate bounds were repaired earlier. These new repairs
are in mcray; the TransportCR checkout remains unchanged.

## Checks implemented

- Gaussian/plateau integrals through independent Simpson quadrature.
- Identity/mass and weighted baryon bookkeeping across A=2..56.
- High-energy event/channel support and special light-nucleus multiplicities.
- K/S product equivalence, Clone mode, unsupported inputs and zero samples.
- Real CMB convolution compared to the analytical Planck integral with
  independent K quadrature; static-background and K/S rate checks.
- No-interaction propagation through the actual engine and output comparator.
- Real PSB propagation grouped by independent primary histories, plus a
  Gaussian-batch uncertainty estimate for fixed composition groups.
- Survival checked against exp(-R56*L), and weighted baryon number checked
  per primary history. No thinning or unrelated interactions are enabled.

Current short suite: 44 CTest entries (37 unit, three integration, two smoke,
and a two-step all-off E2E fixture). Full nuclear CLI support is not claimed.
The comparison uses existing library interfaces and GoogleTest assertions.
All 44 short tests pass natively on macOS arm64 and in Docker Linux amd64
under emulation after the reported repairs. The final container snapshot was
built through bin/test-docker.sh; no physical Linux amd64 host was available.

## Initial experiment

CMB temperature 2.7255 K, fixed at z=0. No redshift losses, background evolution,
neutron decay, magnetic fields, cascades or other interactions. Source is
monoenergetic Fe-56. Time/path convention: H0=67.8 km/s/Mpc, OmegaLambda=0.692,
z_source=0.01 gives L=c*delta_t=43.895883086391059 Mpc; both programs receive
the same elapsed interval/path. TransportCR CalcMode=1 and MomentSource=1.
Source DeltaSourceEnergy is total energy in MeV, not E/A.

Moderate attenuation case: E/A=2.238721138568347e18 eV. MC event rate is
0.027209692731076717 /Mpc; analytical Fe survival is 0.3028876025594345.
Independent validation seeds start at 20261011; 20 batches of 10,000 primaries
were chosen after a separate 20,000-primary pilot lacked sufficient precision
for neutron yield. All 200,000 histories complete after the numerical repair.

| Observable per injected Fe | mcray mean | MC standard error | TransportCR, resolution 16 |
| --- | ---: | ---: | ---: |
| Surviving Fe-56 | 0.302565 | 0.00126893 | 0.302622 |
| A=41..55 residuals | 0.697435 | 0.00126893 | 0.697381 |
| Neutrons | 0.176320 | 0.000930863 | 0.176175 |
| Protons | 1.103520 | 0.00218787 | 1.10401 |

The six predefined groups also include A=2..20 and A=21..40, whose expected
population is negligible in this case. The comparator requires useful MC
precision for populated groups and tests discrepancies against a 1% numerical
allowance (absolute floor 0.001 per primary) plus four measured standard errors.
The four-error cutoff conservatively approximates a Bonferroni family-wise
alpha=0.01 for six approximately normal batch means. This is a compatibility
check, not proof that the true bias is below 1%; confidence intervals and
reference convergence must accompany the result.

## Reference convergence and limitations

TransportCR rates at this E/A:

| energyResolutionFactor | Fe rate /Mpc |
| --- | ---: |
| 4 | 0.0274296 |
| 8 | 0.0272780 |
| 16 | 0.0272297 |

The resolution-4 discrepancy is about 0.8%; resolution 16 approaches mcray
within about 0.074%. The 8-to-16 refinement changes the rate by about 0.177%.
Full source/final spectrum comparisons at resolutions 4, 8 and 16 passed
the provisional statistical criteria. Resolution 16 used microStep=0.5 Mpc;
halving it to 0.25 Mpc preserves Fe/p/n yields at printed precision and also
passes the comparison. The very small A21..40
yield is about 1.5e-10 and is not statistically resolved by the MC ensemble; do not declare all PSB validation complete.
TransportCR uses mean nuclear mass 939.0 MeV, whereas mcray uses the average
of its proton/neutron masses. Account for this small difference when comparing
at fixed gamma versus fixed E/A; do not silently change approved approximate
masses to force agreement.

Raising the grid minimum in a PSB-only run reproduces the broad-grid
resolution-4 Fe/p/n yields exactly at printed precision. The maximum remains
fixed by CFixedEnergyInjection::BaseEmax (total source energy), regardless of
the XML Emax value. The configuration helper documents this limitation.
The all-off reference spectrum is unchanged after the path; output normalization
uses its initial injected Fe population, not a raw flux comparison.

Additional weak/strong CMB comparisons use E/A=1.7782794100389228e18 and
2.511886431509572e18 eV. The weak case passes with 200,000 histories:
Fe survival 0.883860 +/- 0.000593 versus TransportCR 0.883888, proton yield
0.128355 +/- 0.000620 versus 0.128352. The tiny neutron yield is below the
1% relative-precision requirement's applicability threshold (yield 0.01).
The 200,000-history strong case also passes, but a more conservative final
500,000-history ensemble was selected from pilot variance for the small
surviving Fe population. Its batches have seeds 20261121..20261140; independent
batches ran in separate processes. All histories complete and all four
applicable checks pass: Fe survival 0.044692 +/- 0.000292 versus reference
0.044695; p yield 2.077904 +/- 0.001104 versus 2.07696; n yield 1.040750 +/-
0.000795 versus 1.04095. The Fe standard error is below the agreed 1% relative
precision requirement. These additional cases ran natively on macOS; the
moderate case was also validated on Linux amd64.

Native reference runs used the previously built TransportCR executable with
banner 00d5dd5 + local changes. The separate Docker validation built the
pinned committed reference 94eca0c886df80b1d0b1dd09dae95b0483c2273c using
its bin/build.sh. Linux amd64 under emulation passed all four applicable PSB
E2E checks with 200,000 complete histories. Composition agrees with the table
above at printed precision. Outputs are in build/transportcr-validation.NC1e0W;
the updated comparator also passed against these outputs. Large raw histories,
spectra and logs remain outside Git. This validates one moderate-attenuation
scenario, not the entire PSB model or a run on physical amd64 hardware.

## Existing proton pair and nucleon GZK checks

A first proton continuous-loss comparison uses the same CMB and path with
source energy 1e19 eV. mcray retains a proton energy fraction 0.968294;
TransportCR gives 0.968953 at resolution 8 and 0.968660 at resolution 16.
Both satisfy the provisional 1% moment tolerance. The refinement changes the
energy fraction by 0.000293. Proton number is conserved; these checks do not
yet validate the emitted electron/positron spectra or the nuclear model.
Process-specific unit tests cover SOPHIA event nucleon number, energy and
weights, GZK Clone and supported species, pair rate/weight scaling and Clone,
and the continuous-loss proton-only contract.

The initial GZK run at 1e20 eV failed during rate integration with GSL_EROUND
at tolerance 1e-4. A diagnostic run at 1e-3 completes, but is not an acceptance
run and does not justify relaxing the intended tolerance. A regression test
reproduced the failure. GetRateS now uses the 61-point quadrature rule at the
unchanged tolerance; the regression and PSB convolution tests pass. A fresh
20,000-history GZK pilot completes at 1e-4 and satisfies preliminary p/n count
and energy-moment criteria against the existing grid. TransportCR contains
secondary SOPHIA tables only at 20 bins per decade (resolution parameter 2).
Its executable reports a missing table directory for finer grids, sometimes
with exit code zero; check actual final spectra, not just process status.
The user supplied existing finer tables in /Users/ok/propagation/tables/sophia2.
Total p/n cross sections and the 20-bin secondary tables are byte-identical
to the repository data. The helper copies only the selected table density
into an isolated run and writes a SHA-256 manifest; source tables remain intact.
For proton injection, retained reference p+n improves from 0.992648 at 20 bins
to 0.999279 at 50 and 0.999952 at 100. At 100 bins, the proton energy fraction
is 0.486548 versus MC 0.489785 +/- 0.002973; the neutron energy fraction is
0.336110 versus MC 0.332689 +/- 0.002568. Both satisfy the provisional moment
criteria. A neutron-injection pilot also passes at 100 bins: proton/neutron
yields 0.426850/0.573150 versus reference 0.429784/0.570295. Further step
checks halve the TransportCR step to 0.125 Mpc in Docker Linux amd64 and
preserve the proton-injection moments at printed precision. A fresh Linux MC
ensemble of 20,000 histories also passes the count/energy-moment comparison;
its p/n yields are 0.574950/0.425050. Inputs, table hashes, and outputs are
in build/gzk-validation.pbdu10. Secondary spectra and additional energies remain.
No complete GZK physical agreement is claimed.

## Remaining work

- Add CMB+Inoue12 Baseline checks with the agreed frozen z=0 spectrum.
- Extend spectral comparisons beyond integrated composition, with covariance.
- Complete full proton pair/GZK spectral comparisons and additional energies.
- Agree nuclear GZK/pair approximations before implementing nuclear support.
