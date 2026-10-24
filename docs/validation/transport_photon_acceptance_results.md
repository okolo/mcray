# Frozen transport/photon acceptance results

Audited on 2026-10-21. Campaign `20261020T153324Z-d619dbcb` completed on
native Linux amd64 using 16 independent worker processes. **Execution passed;
physical acceptance is inconclusive under the unchanged frozen protocol.**
There are no MC comparisons classified as disagreement. An inconclusive
result does not establish the requested equivalence.

## Provenance and execution

- Public source: `172f25344fff97beed2e41941c32403c46adff53`, clean source tree.
- Protocol: [transport_photons_v1.json](../../tests/crpropa/protocols/transport_photons_v1.json),
  SHA-256 `e5f345a2773745790f7c59e6941b143cb70450e80a42af38a28aac11deef76a9`.
- Archive: `20261020T153324Z-d619dbcb.tar.gz`, SHA-256
  `75216d7cdaf110500d40c4ad2345410298e03f318e5ef4626c965b0492f0a948`.
- Pinned CRPropa source: `73fa865b3304113da2059755a83f5939c1039d68`;
  pinned data: `34b10e27375c0e9894d7c984ac15fb96e12975ef`.
- Retained CRPropa library SHA-256:
  `fae591c411199be359f5fd72cf4c337dbf0154c929460d0bfb338096a7dcfa42`.
- CTest: 144 tests, zero failures or skips. All 272 transport batches finished.
- New executions: 136000 per background/program, 544000 in total. Each group
  contains 96000 main histories and 20000 at each of two refinements.
- Wall time: about 36 min 7 s overall, including 54 s build, 4 min 10 s test
  stage and 30 min 59 s workload (2026-10-20 15:33:24–16:09:31 UTC).

The archive checksum, sidecar/manifest consistency and all 853 payload-file
checksums were verified. The raw histories of every batch were read again:
complete histories, weights, particle accounting, configuration and observable
means agree with the completed-batch records. Batch statistics, Student
intervals and decisions were recomputed from these means using independent
floating-point accumulation. All ten coefficient exports are byte-identical
to the frozen reference controls. No new transport histories were generated.

The statistical audit used NumPy 2.2.6, SciPy 1.15.3 and GitPython 3.2.0,
matching the server package versions. An initial audit with SciPy 1.18.1
found a roughly 1.6e-11 relative difference in the Student critical value
at 47 degrees of freedom; matching SciPy removes that audit discrepancy.
This is insignificant physically and is not a change to any acceptance
criterion. The archived report remains unchanged.

## Frozen decision family

The scope is Fe56, E/A=10^18.6 eV, 1 physical Mpc, frozen z=0 CMB and
CMB+Inoue12, combined PD, corrected decays and both photon sources. See the
[executable contract](../../tests/crpropa/ACCEPTANCE.md) for exclusions and
reference definitions. mcray MC is compared to deterministic M; official
CRPropa MC to deterministic C. Independent I supplies fixed normalizations.

| Numerical setting | Comparisons | Agreement | Inconclusive | Disagreement |
|---|---:|---:|---:|---:|
| Main | 128 | 24 | 104 | 0 |
| First refinement | 128 | 16 | 112 | 0 |
| Second refinement | 128 | 16 | 112 | 0 |
| Total | 384 | 56 | 328 | 0 |

For each background and program the main run accepts bound-nucleus count,
bound A sum, bound Z sum and bound nuclear energy, plus photon count and
energy in [10^15,10^16) eV. Both refinements accept the four bound-nucleus
moments; the smaller samples do not resolve the photon-bin comparisons.
The MC refinements therefore do not establish equivalence of every output
at the requested tolerance. Deterministic grid convergence and successful
execution cannot substitute for those missing statistical decisions.

## Precision of representative main outputs

Entries below are signed MC-minus-matching-reference differences, followed
by simultaneous confidence half-widths, both in percent of the fixed
independent-I mean. All these integral tolerances are 1%. The intervals
belong to the joint 99% family of 384 comparisons; they are not individual
99% intervals. The deterministic envelopes are negligible on this scale.

| Background | Program | Fe56 survival | Total photon count | Total photon energy |
|---|---|---:|---:|---:|
| CMB | mcray | +0.238 ± 1.187% | −0.909 ± 2.975% | −1.541 ± 3.325% |
| CMB | CRPropa | +0.175 ± 1.397% | −0.688 ± 3.157% | −0.190 ± 4.168% |
| CMB+Inoue12 | mcray | +0.332 ± 1.252% | +0.025 ± 3.321% | −0.977 ± 3.675% |
| CMB+Inoue12 | CRPropa | −0.290 ± 1.641% | +0.309 ± 3.136% | +0.460 ± 3.548% |

The lack of precision affects common outputs as well as rare channels.
Although many central estimates differ by less than 1%, acceptance requires
the entire confidence interval plus numerical envelope to fit inside the
budget. None of the rows above meets that requirement. A central difference
greater than 1%, such as CMB mcray photon energy, is also insufficient to
establish disagreement when its interval overlaps the allowed range.

Rare beta-minus outputs are much less resolved: the 96000-history main
samples contain zero electrons for both CMB programs and CMB+Inoue12 CRPropa,
and one electron for CMB+Inoue12 mcray. For CMB the deterministic expected
count is only about 0.21 electrons per main sample. Zero observations are
therefore not evidence of a missing process. Main positron counts are
308/355 on CMB and 473/509 on CMB+Inoue12 (mcray/CRPropa).
Histories, rather than individual correlated descendants, are the sampling
units; counting all emitted photons as independent would understate errors.

## Deterministic reference separation

All 64 M/I coefficient comparisons pass their frozen budgets. Of the 64 C/I
representation diagnostics, 24 agree and 40 disagree. These are differences
between the official CRPropa representation and the independent convolution,
not 40 failed MC implementation tests. The protocol deliberately preserves
those differences: CRPropa MC is tested against C, without correction, and
mcray MC against M. This run reproduces the frozen coefficient exports and
does not introduce a new integration-instability finding.

## Consequences and reproducibility

The campaign does not close its physical acceptance gate, establish the
full Release feature set, or validate other energies/distances. It also does
not demonstrate an implementation defect: every MC decision is either
agreement or inconclusive under the frozen rules.

The reservation accounts for 60000 previous pilot histories, 136000 new
histories and a conservative 4000 allowance for previous short controls per
background/program, filling the 200000 cap. This is reserved accounting,
not a claim that exactly 200000 histories were executed. Do not append seeds,
pool refinements into the main ensemble, relax tolerances, drop rare outputs
or repeat the campaign until it passes. A revised statistical method or
release decision requires separate agreement and explicit treatment of the
already observed samples; it is not authorized by this analysis.

Use [audit_transport_photon_acceptance.py](../../tests/crpropa/audit_transport_photon_acceptance.py)
to repeat the archive audit without simulations. Retain the archive, sidecar,
frozen reference inputs and the resulting audit JSON outside Git. The full
384-row statistics remain in the archived `acceptance-report.json` and the
recomputed audit JSON; this document summarizes their decisions.
