# PSB validation progress

Date: 2026-10-01. The user approved the validation scope described below.
Scientific PSB assumptions remain weighted branching, an A-indexed
chain, and approximate masses. The initial repairs and validation infrastructure
were committed as 6d8f1c4. Subsequent CMB+EBL and spectral extensions are not
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
mcray uses PPP's differential kernel. The 200-history pilot produces
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
independent PPP-kernel test and the expanded suite also passed in Linux amd64 in
a container with the updated test sources, including the native Inoue12
extension-to-zero regression. Full GZK secondary and pair-energy
acceptance ensembles above were run natively; no Linux spectral-acceptance
ensemble is claimed.


## Multi-energy pair-sampler validation

The independent rest-frame PPP-kernel quadrature check now covers k=2.2, 3, 5,
10, and 20 electron masses, including near-threshold photons. Each case uses
10,000 samples at primary proton energy 1e19 eV. Seeds are
20261601 + unsigned(100*k), independently fixed for each energy.
The first two lab-energy moments are compared with a two-dimensional
rest-frame energy/angle integration using 256 and 512 Gauss-Legendre nodes
per dimension. The largest relative quadrature change is 2.78e-8, below the
unchanged 1e-3 convergence requirement. Sampling comparisons retain the
existing 4*SE + 1% recoil-approximation allowance; no tolerance was widened.
All five energy cases pass natively and on Docker Linux amd64.

| k / electron mass | Integrated mean E_secondary/E_primary | Sample mean | Sample SE |
| --- | ---: | ---: | ---: |
| 2.2 | 0.000530991 | 0.000531597 | 0.000001189 |
| 3 | 0.000479482 | 0.000478933 | 0.000002464 |
| 5 | 0.000384079 | 0.000390015 | 0.000003523 |
| 10 | 0.000264988 | 0.000261088 | 0.000003779 |
| 20 | 0.000175160 | 0.000177049 | 0.000004121 |

The second moment at k=20 still has a relative sampling SE of about 9%;
passing this check is not a claim of percent-level tail accuracy. These
checks validate the sampler against the existing differential kernel, not
an independent derivation of the PPP differential cross section. They do not cover the separate asymptotic
sampling branch at k>=20000 or all recoil regimes.
The updated short suite has 51 CTest entries (42 unit, five integration,
two smoke, and a two-step all-off E2E fixture); all pass natively and in Docker
Linux amd64. Logs: /tmp/mcray-qed-multienergy-final.log,
/tmp/mcray-qed-native51.log, and /tmp/mcray-qed-linux51.log.
Production sources and TransportCR sources were not changed.


## Pair-reference grid and step convergence

Further isolated TransportCR runs use the same 1e19-eV proton source,
L=43.895883086391059 Mpc, static 2.7255-K CMB, pair production only, and
Emin=1e13 eV. Resolution denotes ten times as many energy bins per decade.
The reference checkout is 94eca0c886df80b1d0b1dd09dae95b0483c2273c;
the existing native executable SHA-256 is
a931520439aa1758f1fb2260a189e4bed099a421524b2b3f5d31c167c01a2a42.
Its compiled revision banner predates that checkout revision; the binary is
recorded explicitly rather than inferred from the current checkout label.
No TransportCR source or binary changes were made during these checks.

At resolution 32, reducing microStep from 0.125 to 0.0625 Mpc changes the
proton energy loss from 0.0314877288684 to 0.0314876573265, a relative change
of 0.0002273%. The emitted pair energy stays 0.0313733181222 at output
precision. The three applicable pair E2E assertions pass against the existing
2,000-history mcray ensemble; the optional diagnostic-report assertion was
not requested for this check.

To avoid computing empty bins far above injection at resolution 128, Emax
is reduced from 1e25 to 1e20 eV, retaining an extra decade above injection.
TransportCR anchors its grid to integer decades (Ranges.cpp::adjustRanges),
so the retained grid centers do not move. A separate resolution-32 control
with this upper limit produces a byte-identical final spectrum to the full
range. Input directories are /tmp/mcray-transport-pair-spectrum32-step,
/tmp/mcray-transport-pair-spectrum32-bounded,
/tmp/mcray-transport-pair-spectrum64, and
/tmp/mcray-transport-pair-spectrum128-bounded.

At resolution 64 and microStep=0.125 Mpc, the proton energy-loss fraction is
0.0315621582751 and the emitted pair-energy fraction is 0.0315063033646.
The 32-to-64 relative changes are 0.23582% for loss and 0.42209% for emitted
energy, using the finer result as denominator. The energy-balance residual
is 0.17697% of the proton loss. The three applicable pair E2E assertions
pass against the existing mcray ensemble; log:
/tmp/mcray-pair-energy-comparison64.log. This refinement alone does not meet
the 0.2% successive-grid target for secondary energy.

The reusable tests/transportcr/check_pair_convergence.py verifies compatible
physical XML settings and reports loss/emission changes as JSON. It ignores
only resolution, spatial step, upper grid bound, and the two derived source
normalizations written back by TransportCR. Native and Docker Linux amd64
checks correctly pass the step/range controls and reject 32-to-64 grid
convergence at threshold 0.002. Linux tool log:
/tmp/mcray-pair-convergence-tool-linux.log.

The pair E2E assertion now also directly compares the transferred energy
between codes, rather than relying only on the much larger retained proton
energy. It uses the same 1%-of-loss budget plus four batch SEs for both loss
and emitted energy, and checks each SE against 1% of the reference loss.
The resolution-64 comparison passes: loss SE is 4.45245e-9 and emission SE
is 0.000183261 (about 0.58% of the reference loss). The mean proton loss differs
by about 0.43% between codes. Counts and spectral shapes remain unconstrained
across the intentionally different pair-product models.
Native log: /tmp/mcray-pair-energy-comparison64-strict.log.
The strengthened check also passes in Docker Linux amd64 when analyzing
the same native physical ensembles/reference files; no Linux simulation
ensemble is claimed here. Log: /tmp/mcray-pair-flux-comparison-linux64.log.
It correctly rejects the coarser resolution-16 reference: the proton-loss
difference is 0.0003567621422, about 1.14% of the reference loss, exceeding
the 1% loss budget. Earlier resolution-16 passes describe the older retained-
energy and internal-balance checks, not this new direct-transfer assertion.

At resolution 128, microStep=0.125 Mpc, and Emax=1e20 eV, proton loss is
0.0315984378597 and emitted pair energy is 0.0315699614104. The 64-to-128
relative changes are 0.11481% and 0.201641%, respectively. The secondary
energy change narrowly exceeds 0.2%; it is not rounded into a pass.
The reference energy-balance residual decreases to 0.09012% of loss.
The strengthened pair E2E comparison with mcray passes all three applicable
cases; log: /tmp/mcray-pair-energy-comparison128-strict.log.
The native calculation took 5466 seconds, including 805 seconds of propagation.

An additional resolution-32 control with Emax=1.1e19 eV also reproduces the
full-range final spectrum byte for byte. A further resolution-160 refinement
therefore uses this tighter upper limit, retaining Emin=1e13 eV and
microStep=0.125 Mpc. The density increment is 25%, not another doubling;
report its actual successive changes explicitly and do not interpret them
as a rigorous absolute-error bound. Directory:
/tmp/mcray-transport-pair-spectrum160-tight.


A sixth fixed-photon PPP case at k=50 extends the sampler check to the more
energetic CMB photons relevant to a 1e19-eV proton. It retains 10,000 samples
and the same tolerances. The 256-to-512-node changes are 2.17181e-5 and
2.12044e-5 for the two moments, below 1e-3. Integrated mean energy ratio is
9.79214e-5; the sampled mean is 9.80702e-5 with SE=4.48238e-6. The second
moment is 1.77525e-7 by quadrature and 2.10515e-7 +/- 5.94789e-8 in sampling.
Both comparisons pass natively, but the latter has about 28% relative sample
SE; this does not establish precise tail accuracy. The high-k asymptotic
sampling branch remains outside these tests. Native log:
/tmp/mcray-qed-k50-native.log.

## Near-threshold PPP sampling defect and repair

A fixed-photon diagnostic at k=2.0005 electron masses exposes an existing
PPP defect on native macOS and Docker Linux amd64: the approximate total
cross section is positive (diagnostic rate 4.8667e-17 in fixture units), but
DifSigma_r::SetS contains the square-root factor k*(k*(m-1)-2*m)*(m-1).
For this k and m=Mp/Me, the factor is negative. Sampling throws GSL error 5
("could not integrate function"). Diagnostic logs:
/tmp/mcray-pair-threshold-native.log and /tmp/mcray-pair-threshold-linux.log.

PropagationEngine already permits GetSecondaries to return false near a
process threshold. PPP::SampleS now uses this existing null-interaction
contract when k*(m-1)-2*m is not positive, before evaluating the recoil
bounds. The total cross-section approximation and continuous-loss model
are unchanged. The original rate remains a proposal rate; proposals in
the rejected interval produce no secondary and leave the primary intact.
Consequently this repair does not fully reconcile the cross section, CEL,
and sampled secondary production near threshold. The boundary follows
the existing recoil formulas, not a newly derived exact physical threshold.
The contribution of the rejected interval has not been quantified for
arbitrary backgrounds and primary energies.
Previously exceptions triggered recovery attempts and could discard the
cascade after unsuccessful recovery.

The regression is now enabled and covers k=2, k=2.0005, points just below,
at, and above the sampling boundary 2*m/(m-1), and k=2.01. It verifies finite
rates, finite generated energies, empty output for rejected events, and
unchanged primary energy/weight. The regression passes natively after the
repair. This is a new production change in src/lib/PPP.cpp after commit
ce27bd5; earlier statements of unchanged production code refer to the
preceding validation stages.


## Completed PPP refinement and final validation

The original, unmodified TransportCR calculation at resolution 160 finished
in 10220 seconds, including 1342 seconds of propagation. It uses the previously
verified upper bound Emax=1.1e19 eV, Emin=1e13 eV, and microStep=0.125 Mpc.
No TransportCR source or binary was modified for this refinement.

| Resolution factor | Proton energy-loss fraction | Emitted pair-energy fraction | Balance residual / proton loss |
| --- | ---: | ---: | ---: |
| 64 | 0.0315621582751 | 0.0315063033646 | 0.17697% |
| 128 | 0.0315984378597 | 0.0315699614104 | 0.09012% |
| 160 | 0.0316061396434 | 0.0315836396987 | 0.07119% |

From 128 to 160, the relative changes using the finer result as denominator
are 0.024368% for proton loss and 0.043308% for emitted pair energy. Both
pass the unchanged 0.2% successive-refinement criterion. This increment is
25% in energy-bin density; the result is not an absolute discretization-error
bound. Machine-readable check: /tmp/mcray-pair-grid160-convergence-check.json.

The repeated mcray ensemble has 20 independent batches of 100 proton
histories, injection energy 1e19 eV, static CMB, and distance
43.895883086391059 Mpc. It passes all three applicable E2E assertions against
the resolution-160 reference; six assertions for other processes or optional
spectral output are skipped. Log:
/tmp/mcray-pair-energy-comparison160-fixed.log.
Mean proton loss is about 0.0316968 and emitted pair energy about 0.0316167.
Their batch SEs are 4.45245e-9 and 0.000183261, respectively. The loss differs
from TransportCR by about 0.287%; the emitted-energy comparison retains its
1%-of-reference-loss plus four-SE budget. Pair multiplicities and spectral
shapes are not required to agree between the distinct secondary models.

After the threshold repair, all 53 CTest entries pass on native macOS and
Docker Linux amd64 (44 unit, five integration, two smoke, and two all-off E2E
fixture entries). Logs: /tmp/mcray-pair-threshold-final-native53.log and
/tmp/mcray-pair-threshold-final-linux53.log. These are short-suite runs;
physical acceptance ensembles were generated natively.

Repeating the 2000-history physical ensemble with identical seeds after the
repair produces byte-identical physical data rows. Their SHA-256 is
2ba5e41167cce9065b30a07d7e65bd62c71e9da3d702d36e827438356d60037e.
The whole TSV files differ because the updated harness adds primary-species
metadata. This controls the selected scenario only; it does not establish
that threshold rejection is negligible for all inputs.


Following the terminology clarification, test descriptions name PPP directly
and the parameterized sampler suite is ProtonPppSampling. The final sources
were rebuilt and all 53 entries passed again natively and in Docker Linux
amd64. Logs: /tmp/mcray-ppp-final-native53.log (27.89 seconds) and
/tmp/mcray-ppp-final-linux53.log (86.46 seconds, following a fresh image build).
The TransportCR executable used for the resolution-160 run has SHA-256
 a931520439aa1758f1fb2260a189e4bed099a421524b2b3f5d31c167c01a2a42.


## Approved GZK spectral method and first pilot

On 2026-10-02 the user approved the following GZK spectral method: seven
secondary species, four fixed intervals from 1e14 to 1e20 eV, simultaneous
99% Bonferroni Student intervals, a 5% spectral budget, and three outcomes.
The old 1% integrated moment criteria remain unchanged. No new physical
simulation or retrospective spectral acceptance is claimed here.

Test-only shared statistics and five unit cases are implemented. The
comparison runner writes a pilot report with full-history occupancy and
batch uncertainties. It does not assert global spectral acceptance.
For 20 batches and 28 observables the Student critical value is
4.334545147994069. Native short-suite validation passes all 58 entries
(49 unit, five integration, two smoke, and two all-off fixture entries).
Log: /tmp/mcray-gzk-statistics-native58.log.

Two original TransportCR controls were run with the same injection,
background, path, Emin=1e14 eV and Emax=1e26 eV:
/tmp/mcray-transport-gzk-spectrum50 (50 bins/decade, step 0.125 Mpc),
and /tmp/mcray-transport-gzk-spectrum100-step (100 bins/decade,
step 0.0625 Mpc). Their total calculation times are 10 and 79 seconds.
XML comparison against the original spectrum100 case finds only energy
resolution plus rewritten source normalizations for the grid control,
and only microStep for the step control. SOPHIA table hashes are saved
in each case's sophia-manifest.json. TransportCR source and binary are unchanged.

The 150000-history pilot gives 28 fixed-bin outputs. With extraction by
whole-cell center energy, the absolute grid/step-change sum exceeds 0.2%
of reference yield for 23 observables. All 28 have positive observations
in all 20 batches, but this does not establish Student-interval validity.
Only five interval outcomes are conditionally compatible with the approved
agreement rule; the remaining 23 are inconclusive, not disagreements.
Logs: /tmp/mcray-gzk-wide-bin-pilot.log and
/tmp/mcray-gzk-wide-bin-pilot.tsv. The four applicable GoogleTest cases pass;
six unrelated cases are skipped. Reporting success is not spectral acceptance.

A separate read-only extraction diagnostic considers log-uniform population
within each reference cell and apportions a crossing cell by log-energy
intersection. The 50-to-100 grid changes then range from about 0.00085% to
0.49753%, instead of changes reaching several percent with center assignment.
Halving the spatial step changes these reconstructed yields by at most
0.00000179% at printed precision. Nine of 28 grid changes still exceed 0.2%.
This reconstruction is not yet adopted: it introduces an explicit
within-cell assumption and has been submitted for user approval.
Diagnostic script: /tmp/mcray-gzk-cell-diagnostic.py.

Further history-count estimates from the initial extraction are provisional:
several intervals have no available budget margin after residual and numerical
changes, so increasing MC statistics alone cannot resolve the comparison.
First settle extraction and reference convergence, then plan independent
acceptance statistics without searching for passing seeds.


Final short-suite runs pass all 58 entries on native macOS and Docker Linux
amd64. Logs: /tmp/mcray-gzk-final-native58.log (29.31 seconds) and
/tmp/mcray-gzk-statistics-linux58.log (88.91 seconds). The latest comparator,
including the guard rejecting incomplete convergence arguments, was separately
rebuilt in a disposable Linux container. Its all-off fixture and GZK pilot
report pass, and incomplete controls correctly return failure. Log:
/tmp/mcray-gzk-pilot-linux.log. The Linux pilot analyzes the native ensemble
and references; no new Linux physical ensemble was generated.
All changes in this stage are test tools and documentation. Production
mcray and TransportCR source files remain unchanged after 94322d7.


## Approved GZK cell reconstruction and repeated pilot

The user approved constant dN/dlog(E) within each reference cell on
2026-10-02. The comparison runner now shares a crossing cell by its logarithmic
intersection with the four unchanged intervals and preserves both tails.
Only the seven secondary species are reconstructed. The printed reference
grid is validated and its bins-per-decade density is inferred; no TransportCR
source, tables, or physical model were changed.

Three additional unit checks cover half/quarter crossing-cell fractions,
count conservation including tails, and grid rounding versus missing cells.
All eight statistical/reconstruction unit checks pass natively. The pilot
was repeated on the same 150000 histories and the same original 50/100-bin
and spatial-step controls. Logs: /tmp/mcray-gzk-cell-pilot.log and
/tmp/mcray-gzk-cell-pilot.tsv. Reporting passes all four applicable checks;
this is not a new acceptance ensemble or a global spectral agreement result.

Thirteen of 28 conditional interval outcomes are agreements; fifteen are
inconclusive and none are disagreements. Nine observed grid/step-change sums
still exceed the unchanged 0.2% numerical threshold:

| Species | Energy interval, eV | Refinement change / reference |
| --- | --- | ---: |
| Photon | [1e14,1e17) | 0.44520% |
| Electron | [1e14,1e17) | 0.32716% |
| Electron | [1e17,1e18) | 0.20399% |
| Electron | [1e19,1e20] | 0.40546% |
| Positron | [1e14,1e17) | 0.47960% |
| Electron antineutrino | [1e14,1e17) | 0.25988% |
| Electron antineutrino | [1e19,1e20] | 0.49752% |
| Muon neutrino | [1e19,1e20] | 0.45741% |
| Muon antineutrino | [1e19,1e20] | 0.20101% |

Values just above 0.2% are not rounded into passes. The available independent
SOPHIA secondary tables have 20, 50, and 100 bins per decade; denser tables
or the original generator have been requested before further refinement.

For the upper electron interval the pilot's simultaneous half-width is
about 14.71% of reference; its observed residual is about -2.60%. A heuristic
SE proportional to N^(-1/2), keeping 20 batches and the observed residual/
refinement envelope fixed, gives a total-history factor about 54.43, or
8.17 million histories. This is an extrapolation, not an approved run size,
a confidence guarantee, or a promise that acceptance will pass. Further
reference refinement may change the margin. No new acceptance ensemble
has been launched while reference convergence remains unresolved.


The final reconstruction sources pass all 61 CTest entries (52 unit, five
integration, two smoke, and two all-off fixture entries) on native macOS
and Docker Linux amd64. Logs: /tmp/mcray-gzk-cell-native61.log (27.52 seconds)
and /tmp/mcray-gzk-cell-linux61.log (91.84 seconds). A separate Linux
comparison-runner check also passes the reconstructed pilot report on the
same native input data: /tmp/mcray-gzk-cell-pilot-linux.log. This verifies
analysis portability, not a new Linux physical simulation ensemble.


## Independent SOPHIA generator recovered on 2026-10-07

The source was found at ~/propagation/sources/sophia and frozen outside Git
in ../validation/sophia-generator-20261007. Source-file hashes, image identity,
binary checksum, compiler versions, and pilot logs are retained there.
The new helper tests/transportcr/build_sophia_generator_docker.sh reproduces
the isolated Linux amd64 build without modifying the original checkout.

P/n diagnostic tables use 200 bins per decade, kmax=1e7 GeV, 100 trials per
cell, and seeds 2026100701/2026100702. They completed in 19.29/19.43 seconds
while run concurrently. All 18 secondary files passed format checks;
photon-energy row counts are 1564/1565. Repeating the proton run with the same
seed produced byte-identical nine secondary tables and the cross-section file.
These pilots do not validate physical spectra or numerical convergence.

Existing table statistics are 1.38 million, 2 million, and 1 million trials
per cell at 20, 50, and 100 bins per decade, respectively. Thus grid differences
also contain table sampling uncertainty. A full 200-bin p/n generation with
one million trials per cell and kmax=2e6 GeV has been proposed for approval;
rough runtime is 45--55 hours in the current Docker configuration, with no
completion guarantee. Full generation and acceptance ensembles are not running.
See sophia_table_generation.md for the statistical and provenance limitations.


The user approved the proposed full-generation parameters on 2026-10-07.
Both p/n Docker Linux amd64 jobs started at 20:33:59 UTC, with seeds
2026100711/2026100712. They are running as detached containers named
mcray-sophia200-p-20261007 and mcray-sophia200-n-20261007. Commands,
source/image/binary hashes, and completion markers are retained under
../validation/sophia-generator-20261007/full-generation.json and
full200-p/full200-n. This records job launch, not completed table generation.

The diagnostic TransportCR run with the 100-trial, 200-bin pilot tables
completed successfully in 149 seconds (130 seconds propagation), confirming
loader compatibility only. Its XML, table manifest, logs, and outputs are
../validation/gzk-smoke200-20261007. The original TransportCR remains unchanged.

The new validate_sophia_tables.py validates all 18 secondary tables and two
cross-section files and refuses non-successful/incomplete runs. The complete
pilots passed; the running full-generation directories were correctly rejected
without publishing a validation manifest. Table shape/range validation and
same-seed reproducibility do not establish the 0.2% reference accuracy.
Replica independence requires further review of the generator's warm-up-only
seed implementation; that implementation was preserved unchanged.


## Extended PSB CMB+Inoue12 energy checks, 2026-10-07

Two additional Fe56 source energies were checked at fixed L=1 Mpc:
E/A=1e18 and 10^18.5 eV. All other interactions, decay, redshift losses,
and thinning are disabled. CMB temperature is 2.7255 K and the approved
Inoue12 Baseline is frozen at z=0. Its SHA-256 remains
197a50e810b23538d506763be2664cfa9ffead2d8f0747df90665c466e2ea993.
No TransportCR source or original background tables were modified.

Raw XML, spectra, per-history outputs, failed diagnostics, comparison logs,
and JSON summaries are outside Git in
../validation/psb-inoue-extension-20261007. The completed-scenario manifest
records binary hashes, controls, seeds, and raw-file hashes. Each accepted
ensemble contains 20 independent batches of 10000 complete histories;
first seeds are 2026100781 and 2026100811. These were fixed before the final
ensembles. Diagnostic 2000-history pilots are excluded from final estimates.

### Rate-integration defect and bounded fix

With the previous 15-node quadrature in ContinuousBackgroundIntegral::GetRateK,
the high-energy ensemble encountered GSL_EROUND for the A30 daughter and
lost primary histories. Its partial output is retained and excluded. Direct
Fe56 probes at tolerance 1e-5 also failed at both energies.

GetRateK now explicitly uses GSL_INTEG_GAUSS61 on both redshift planes,
matching the existing GetRateS choice. Relative tolerances, error handling,
PSB cross sections, weighted branching, A-chain, and masses are unchanged.
This resolves the observed A30 failure without accepting failed integrals
or relaxing tolerances. Both final ensembles were rerun using the original
seeds with the corrected quadrature. All 200000 histories per scenario are
present; per-history weighted baryon conservation and static attenuation
checks pass. Low-energy grouped results are unchanged at the displayed
precision compared with its complete earlier ensemble.

Fe56 rates at the working tolerance 1e-4 and background spacing 0.01 are
0.05608426304/Mpc and 0.4483519928/Mpc. Switching quadrature changes the
previous rates by about 0.00040% and 0.00045%, respectively. Refining spacing
to 0.0025 with tolerance 1e-5 gives 0.05608752066/Mpc and
0.4483531450/Mpc. A30 also supports those controls; its high-energy rate
changes by about 0.0181% between the coarsest and finest probe.

The fix does not establish arbitrary-tolerance robustness: A45 at the higher
energy still raises GSL_EROUND with spacing 0.01 and tolerance 1e-5 on
both macOS and Linux amd64. Its
spacing-0.005 probe succeeds. The integration regression therefore covers
all A=2..56 at the actual working tolerance 1e-4, checks K/S equivalence and
spacing refinement within 0.05%, and checks the stricter 1e-5 tolerance for
Fe56/A30 specifically. The failed stricter all-chain diagnostic remains
recorded; no full-chain claim is made at tolerance 1e-5.

### TransportCR refinement

Reference runs use 40, 80, and 160 bins per decade, microStep=0.01 Mpc,
and a separate 160-bin run with microStep=0.005 Mpc. All four runs per
energy completed successfully. Delta-source deposition into adjacent
TransportCR total-energy cells is retained and included in grid refinement;
it is not treated as an exact Monte Carlo delta at a bin center.
Composition is normalized by initial Fe56 count from each run.

The largest relative 80-to-160 change across the six prescribed groups is
0.1095% at the lower energy and 0.1139% at the higher energy. Spatial-step
changes are below 0.000002% at printed precision. The 40-to-80 changes
reach 0.2674% and 0.3213%, showing why the coarser reference is inadequate
for all groups under the provisional 0.2% refinement target. These are
observed refinement changes, not rigorous upper bounds on reference error.

### Composition and precision

Values below are weighted yields per primary. SE is computed from the
20 complete-history batch means, preserving correlated daughter weights.
The reference uses 160 bins per decade and microStep=0.01 Mpc.

| E/A, eV | Group | mcray mean | Batch SE | TransportCR |
| --- | --- | ---: | ---: | ---: |
| 1e18 | Fe56 | 0.945080 | 0.000516089 | 0.945453434 |
| 1e18 | A41..55 | 0.0549191913 | 0.000516068 | 0.0545456184 |
| 1e18 | A21..40 | 8.08750e-7 | 1.93490e-7 | 7.10569e-7 |
| 1e18 | A2..20 | 0 | 0 | 1.74491e-13 |
| 1e18 | p | 0.0627884663 | 0.000597946 | 0.0622178312 |
| 1e18 | n | 0.0026777145 | 0.000181606 | 0.00262464915 |
| 10^18.5 | Fe56 | 0.638665 | 0.001027959 | 0.638644969 |
| 10^18.5 | A41..55 | 0.361302579 | 0.001027875 | 0.361324466 |
| 10^18.5 | A21..40 | 3.24204158e-5 | 1.64919730e-6 | 3.06189696e-5 |
| 10^18.5 | A2..20 | 2.09250e-10 | 2.09250e-10 | 5.61354e-11 |
| 10^18.5 | p | 0.459151275 | 0.001466146 | 0.458558261 |
| 10^18.5 | n | 0.0326608555 | 0.000446027 | 0.0319463301 |

All six discrepancies at each energy satisfy the existing exploratory rule
abs(MC-reference) <= max(0.001, 0.01*reference) + 4*SE. This is a discrepancy
check, not proof of a bias below 1%. Fe56, A41..55, and proton SEs satisfy
the 1% relative-precision requirement at both energies. Neutron relative
SEs are 6.919% and 1.396%; these yields remain statistically inconclusive.
The rare A21..40 and A2..20 groups also lack 1% precision. No sampled output
in the low-energy A2..20 group implies zero empirical SE, not a precision
or agreement guarantee for a nonzero reference yield.

The lower-energy GoogleTest comparison passes its four applicable checks.
The higher-energy comparison passes completeness, baryon conservation,
and static attenuation, but correctly fails the composition precision guard
for neutrons (reference yield exceeds 0.01). That failure is retained and
classified as insufficient statistics, not hidden or rounded into success.
Neither scenario constitutes full PSB validation at 1% precision for all
outputs. The agreed 200000-history cap was preserved; no acceptance-driven
seed replacement or extra histories were added.

The final native short suite passes all 62 CTest entries (52 unit, six
integration, two smoke, and two all-off fixture checks). Log:
/tmp/mcray-psb-extension-native62-final.log. Linux amd64 verification is
complete: all 62 entries also pass under Linux amd64 emulation (88.65
seconds). Log: /tmp/mcray-psb-extension-linux62-final.log. The Docker
image contains the corrected production integrator; its final integration
test source was mounted, copied into the isolated container, and rebuilt
before running the full suite. The initial stricter all-chain diagnostic
failed on both platforms and is preserved separately.

The Linux comparison runner also reproduces the native outcomes on the
same saved physical inputs: low-energy comparison passes, high-energy
comparison fails only the neutron precision guard. This verifies analysis
portability, not a separate Linux physical ensemble. Logs are
low/comparison-linux.log and high/comparison-linux.log in the raw-run
directory.
