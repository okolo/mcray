# PSB validation progress

Date: 2026-10-01. Plan: nuclear_validation_plan.md, approved and committed
as 327bd46. Scientific PSB assumptions remain weighted branching, an A-indexed
chain, and approximate masses. The initial repairs and validation infrastructure
were committed as 7f4b45f. Subsequent CMB+EBL and spectral extensions are not
yet committed.

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

- Extend CMB+EBL controls to further energies/paths after accounting for the
  cost of weighted branching without thinning.
- Extend spectral comparisons beyond integrated composition, with covariance.
- Complete full proton pair/GZK spectral comparisons and additional energies.
- Agree nuclear GZK/pair approximations before implementing nuclear support.

## Subsequent CMB+EBL and spectral checks

The user selected Inoue12 Baseline. Both repositories have identical
EBL_Inoue/EBL_proper_baseline.dat with SHA-256
197a50e810b23538d506763be2664cfa9ffead2d8f0747df90665c466e2ea993.
The table starts at z=0.01. Production Inoue12 already enables
MatrixBackground's extend-to-zero behavior: both its evaluator and MinZ()
include z=0. CompoundBackground therefore retains it below the first table
plane. The earlier claim that this component was excluded was incorrect.
The test-only FrozenInoueBaseline fixture freezes evaluation at z=0; it does
not repair a production lower-bound defect. Direct production and frozen
integration regressions check the compound background at the start, midpoint
and end. Production EBL evolution was not changed.

Initial incomplete long-path pilots are excluded; the redundant MinZ override
did not change the production background's already-correct lower bound.
The 43.8959-Mpc static EBL pilot was stopped because weighted PSB branching is
too expensive without thinning; no acceptance result is inferred from it.
The planned thin-path control instead uses L=1 Mpc, the same E/A as the
moderate CMB case, and no thinning. The user capped it at 200,000 histories
(20 batches of 10,000; seeds 20261501..20261520), accepting an inconclusive
rare neutron yield rather than relaxing precision requirements.

| Observable | MC mean | Batch SE | TransportCR resolution 16 |
| --- | ---: | ---: | ---: |
| Fe-56 survival | 0.884800 | 0.000654 | 0.885330 |
| A41..55 yield | 0.115192 | 0.000654 | 0.114663 |
| Proton yield | 0.140148 | 0.000973 | 0.139327 |
| Neutron yield | 0.012049 | 0.000213 | 0.011559 |

All histories and weighted baryon checks pass; mean discrepancies meet the
provisional agreement criteria. The composition test still returns failure
for insufficient neutron precision, as intended. Its relative standard error
is about 1.84% of the reference yield, exceeding 1%. This scenario is not a
complete pass. The initial Fe rate is 0.1218120771/Mpc in mcray, versus
0.1217942844/Mpc inferred from TransportCR survival. At resolution 8 the
reference rate is 0.1218573452/Mpc (8-to-16 change about 0.052%).
Halving microStep from 0.01 to 0.005 Mpc preserves Fe, p and n yields at
printed precision. These reference refinements pass numerical checks despite
the separately inconclusive neutron MC precision.

The runner now optionally retains individual terminal particles, including
their primary-history membership. A fresh GZK proton ensemble of 150,000
histories (20 batches of 7,500; seeds 20261401..20261420) passes nucleon and
all seven EM/neutrino secondary count/energy-moment checks. The reference uses
100-bin SOPHIA tables, microStep=0.125 Mpc, and a lower energy bound of 1e14 eV
to retain secondary products. The electron yield is 0.131927 +/- 0.000949
versus 0.131596; photon energy fraction is 0.074486 +/- 0.000298 versus
0.074350. Diagnostic spectral bins are written to
/tmp/mcray-gzk-acceptance-spectra.tsv. A global covariance-aware spectral
acceptance decision is still pending; passing moments is not full shape validation.

Proton pair products differ physically between the codes. TransportCR
Ppp.cpp::P_PPP uses an E^(-7/4) approximation normalized by continuous losses;
mcray uses PPP's differential QED kernel. The 200-history pilot produces
electron/positron yields 90.59/92.19 versus 83.624/83.624 in TransportCR.
Do not repair this difference by silently replacing either model. The user
approved preserving both, comparing energy balance, and independently
checking the mcray sampler. A fixed rest-frame photon test at k=3 electron
masses integrates PPP::DiffSigma over electron energy and angle with fixed
Gauss-Legendre quadrature, then compares the first two lab-energy moments
with 2,000 generated samples. Quadrature refinement 64-to-128 nodes and
sampling checks pass. This checks the coordinate transformation and sampler
against the existing kernel; it does not independently derive that kernel or
validate every high-k/recoil approximation. A separate 2,000-history ensemble
(20 batches of 100; seeds 20261621..20261640) passes pair energy balance:
fraction lost by protons 0.0316968, generated e-/e+ energy fraction 0.0316167,
batch standard error of their difference 0.0001833. TransportCR gives losses
0.031340 and generated energy 0.031127 at resolution 16. The latter discrepancy
is below 1% of the loss and reflects grid/normalization approximations. Four
applicable pair E2E checks pass, including complete histories, primary moments
and the diagnostic spectral report. Equal e-/e+ counts or shapes are not
asserted because the user explicitly retained the different physical models.

Refining the pair reference to resolution 32 and microStep=0.125 Mpc also
passes the four applicable checks: proton loss 0.0314877 and emitted energy
0.0313733. The reference's emitted-versus-lost difference decreases from about
0.68% to 0.36% of the loss. The emitted-energy change between resolutions
16 and 32 is about 0.79%; further refinement remains necessary for the
proposed 0.2% numerical target on secondary energy. These initial passes do
not establish final high-precision pair-spectrum/loss convergence.

The expanded short suite has 47 CTest entries: 38 unit, five integration,
two smoke and a two-step all-off E2E fixture. Native macOS checks pass; the
independent QED test and the expanded suite also passed in Linux amd64 in
a container with the updated test sources, including the native Inoue12
extension-to-zero regression. Full GZK secondary and pair-energy
acceptance ensembles above were run natively; no Linux spectral-acceptance
ensemble is claimed.
