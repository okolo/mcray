# Nuclear transport and secondary-photon pilot results

Analyzed 2026-10-17. These are audited pre-acceptance pilots, not final
statistical acceptance. Production source and physical models are unchanged.
The workflow is documented in [the server instructions](../../tests/crpropa/TRANSPORT_PHOTONS.md).

## Run and integrity

The native Linux run `20261016T102021Z-03eb6c57` requested and checked out
mcray 82df04adc32420fd68c12656e23873d7977d17fa, with CRPropa3 pinned at
73fa865b3304113da2059755a83f5939c1039d68 and CRPropa3-data pinned at
34b10e27375c0e9894d7c984ac15fb96e12975ef. Its archive SHA-256 is
`dc0318f5f21cc4cec1b64bd62a747768e3b7ef9123a6afebd81b5f4909dca07d`.
The archive, sidecar/embedded manifests and all 5055 payload hashes agree.
Frozen adapter sources, the original 86-file input/tool bundle and both
prepared photon-reference payloads also pass their checksum checks.

Recorded platform: native Linux x86_64, kernel 5.15.0-191, glibc 2.35,
Python 3.10.12, NumPy 2.2.6, SciPy 1.15.3 and GitPython 3.2.0. The run
reports 256 host CPUs; each workload process is configured with one OpenMP
thread. Release build, all 129 CTest cases and seven workload adapter tests
pass. The recorded native probe SHA-256 is
`ec09daf429b0ca8eb57bc72c1d63274e9fd75b845c592a1c39d1e03aeedafc43`;
the reference adapter SHA-256 is
`2aeea9808a3e1c18a6e77df130509d2d6dd8d45df74bfb5f20a516998c148479`.

The completed initial matrix has 12 cases, two programs, three numerical
levels and 20 batches of 1000 primaries: 1440 jobs and 1440000 history
executions. Each program/case/level has 20000 primaries. Refinement levels
reuse seeds inside a program/case; the 1.44 million executions must not be
treated as independent observations when comparing numerical levels.

All raw TSVs were independently checked: complete history IDs and completion
footers, unit weights, finite nonnegative energies, valid species, baryon
number 56 and total charge 26 including charged leptons, and zero photons
when both sources are Off. All 2709995 particle rows, including 330314 photon
rows across all settings, pass. Counts reconstruct exactly; independently
summed batch energy means agree within 1.28e-14 relative rounding. Recomputing
statistics from the verified archived batch means reproduces all 4709
observable records, 2509 program comparisons and 3225 paired-refinement
records. No discarded histories or incomplete jobs were found.

## Physical scope

Fe56 is injected at E/A=10^18.6 eV and propagated for 1 physical Mpc on each
shared static z=0 background: CMB and CMB+Inoue12. Cases are TALYS with PD
photons Off/On and TALYS+decays with both photon sources Off, PD only,
decay only and both On. There are no EM interactions, cosmological losses,
thinning, free-neutron decay, GZK or PPP in this fixture.

The approved CRPropa photon saturation policy, eight explicit missing-photon
transitions, corrected exclusive decay data and non-debited photon energy
approximation remain unchanged. Total energy is retained diagnostically;
exact event energy conservation with beta decays/photons is not asserted.
This does not validate propagated EM cascades or T5/T6 alternatives.

The finest mcray setting is background log10 spacing .0025 and tolerance
1e-6; finest CRPropa maximum propagation step is L*.003333333333333333.
Coarser settings follow the frozen workflow; no setting was selected by
repeating samples until agreement.

## Survival: separate reference rates from transport sampling

The source Lorentz factor is 4.2799186583774633e9. A separate bounded
calculation evaluates the unchanged official convolution and independent
interval Gauss-Legendre quadrature at this exact gamma, retaining the shared
background and agreed cross-section interpolation. Parent-rate metadata
from the finest native jobs supplies the mcray comparison; no new propagation
ensemble or production binary was run locally.

| Background | mcray finest, 1/Mpc | Independent order 64, 1/Mpc | Official direct convolution, 1/Mpc | Official runtime gamma-table interpolation, 1/Mpc |
| --- | ---: | ---: | ---: | ---: |
| CMB | 0.5661016388 | 0.5661137858 | 0.5726997403 | 0.5831613545 |
| CMB+Inoue12 | 0.6460028158 | 0.6460078238 | 0.6526894470 | 0.6631022901 |

Against the independent rate, official direct quadrature is higher by
1.16336% on CMB and 1.03429% on the combined field. Gamma-table interpolation
adds 1.84797% and 1.61188% of the independent rate, respectively. The total
runtime reference differences are therefore 3.01133% and 2.64617%.
The finest native rates differ by -0.002146% and -0.000775%.

This is **not** the approximately 8% high-gamma missing-boundary effect.
At this pilot gamma the endpoint contribution is zero on CMB and only
2.30e-10 of the independent rate on the combined field. Refining the official
rest grid by factor 64, including the endpoint diagnostically, leaves only
4.12e-6 and 3.65e-6 relative residuals. The official code/data remain intact.
See [the earlier rate report](photodisintegration_physical_validation.md#official-control-differences-at-high-gamma)
for the different high-gamma regime.

Fe56 is stable in the selected decay model, so its survival is exp(-R*L)
regardless of the photon toggles. The two runtime rate representations
predict different survival before any MC sampling:

| Background | Independent-rate survival | Official-table survival | mcray pilot, photons/decays Off | CRPropa pilot, photons/decays Off |
| --- | ---: | ---: | ---: | ---: |
| CMB | 0.5677275 | 0.5581311 | 0.56805 +/- 0.00424 | 0.55645 +/- 0.00457 |
| CMB+Inoue12 | 0.5241340 | 0.5152504 | 0.52550 +/- 0.00348 | 0.51450 +/- 0.00306 |

Errors in this table are batch standard errors, not simultaneous acceptance
intervals. The deterministic relative survival difference is about 1.72% on
both fields, already larger than the agreed 1% integrated budget. Across all
72 case/program/level groups, the largest standardized survival residual
against each program's own recorded rate is 2.66 batch standard errors.
This is descriptive consistency with the lifetime law, not a frozen-family
acceptance verdict. Interprogram differences must not be attributed to a
transport defect without accounting for these reference representations.

## Photon yields and energies

The finest-setting photon number per primary is:

| Background | Nuclear processes / photon sources | mcray | CRPropa | Difference relative to CRPropa | Combined batch SE, photons/primary |
| --- | --- | ---: | ---: | ---: | ---: |
| CMB | TALYS / PD | 0.41070 | 0.42545 | -3.47% | 0.00960 |
| CMB | TALYS+decay / PD only | 0.41780 | 0.44285 | -5.66% | 0.00918 |
| CMB | TALYS+decay / both | 0.41555 | 0.43475 | -4.42% | 0.00893 |
| CMB+Inoue12 | TALYS / PD | 0.48175 | 0.49995 | -3.64% | 0.00950 |
| CMB+Inoue12 | TALYS+decay / PD only | 0.48425 | 0.50505 | -4.12% | 0.00904 |
| CMB+Inoue12 | TALYS+decay / both | 0.48320 | 0.49910 | -3.19% | 0.00881 |

With both sources enabled, mean photon energy summed per primary is
1.70106e15 versus 1.82055e15 eV on CMB (-6.56%, combined SE 4.17e13 eV),
and 2.00111e15 versus 2.02558e15 eV on the combined field (-1.21%, combined
SE 4.55e13 eV). Complete primary histories, rather than individual photons,
are the independent statistical units. These offsets mix rate/branch/photon
probability interpolation and sampling uncertainty. They do not by themselves
identify a photon-emission bug, establish a final disagreement verdict or
justify attributing the photon offsets to the 8% endpoint effect. The raw
exploratory spectra remain available for fixed-bin design.

Decay-only photon cases are particularly sparse. At the finest level:

| Background | mcray photons / photon-producing histories | CRPropa photons / photon-producing histories |
| --- | ---: | ---: |
| CMB | 33 / 33 | 32 / 32 |
| CMB+Inoue12 | 34 / 31 | 34 / 32 |

Each entry represents 20000 primaries. Between one and five of the 20 batches
contain no photons. Count relative standard errors range from 11.9% to
18.7%; energy errors are also large. Similar counts do not establish
percent-level yield or spectral agreement. Zero or unobserved pilot bins
cannot be dropped from a future family or assigned exact-zero uncertainty.

## Numerical refinement and precision limits

The native refinements change few sampled histories; for both-photon cases,
level 1 to level 2 photon counts/energies are identical in these samples.
This is finite-sample evidence, not proof of a zero numerical envelope.
CRPropa step refinement remains noisy: for both sources, the last count
change is +0.00200 +/- 0.00553 photons/primary on CMB and
+0.00825 +/- 0.00617 on the combined field. Here errors are paired batch
standard errors. Equal seeds do not preserve identical event draws after
step/module changes. The pilot does not establish a sub-1% reference
transport envelope; apparent refinement differences cannot be treated as
pure numerical errors without their uncertainty.

A planning-only projection holds the observed variance fixed and increases
20000 to 200000 histories per program. Even for one observable, zero true
program difference, zero numerical envelope and a 99% Student interval with
19 degrees of freedom, predicted photon-count half-widths for both sources
are 1.86% (CMB) and 1.60% (combined); photon-energy half-widths are 2.07% and
2.03%. These already exceed the 1% integrated budget before a family
correction. They are rough pilot projections, not guaranteed precision or
approved final sizes. Rare decay-only outputs have a more severe deficit.
Do not increase the approved cap or weaken budgets; insufficient precision
must remain inconclusive. Broad spectral intervals need separate pilot-based
design under the unchanged 5% budget.

## Cost, retained evidence and remaining gate

The workload lasts 5 h 33 min 20 s. Summed native-job times are 19312.04 s
(5 h 21 min 52 s); reference jobs take 480.65 s (8 min 1 s). Setup/build and
preparation account for the remainder. These sequential timings describe
this campaign; direct runtime convolution and preconvolved tables perform
different work, so they are not a controlled code-speed benchmark.
Including the initial mcray build and CTest, the recorded run lasts about
5 h 37 min 21 s before archive compression.

The original archive is in ../mcray-test-data. Extracted results, complete
integrity/raw-history/statistics audits, source-rate decomposition, endpoint
controls and planning-only precision projections are retained outside Git in
../validation/transport-photon-analysis-20261017. Other session archives in
the transfer directory were not inspected or modified. Analysis documentation
was written on talys; no production changes, new ensembles, commits or pushes
are part of this analysis.

The independent pilot stage is complete for this initial matrix. Final
families, spectral bins, normalizations, fresh seeds, sizes and numerical
reference treatment are not frozen. The deterministic reference-rate gap,
reference step precision and sparse photon components must be resolved or
explicitly classified before any final acceptance claim. The 99%/1%/5%
budgets, 200000-history cap and physical-model limitations are unchanged.
