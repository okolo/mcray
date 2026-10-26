# Independent weighted-transport pilots: 2026-10-25 server results

Analysis date: 2026-10-26. Archive: `20261025T062534Z-0c98cfe7.tar.gz`.
Archive SHA-256:
`84b6c4d8485aeb9e2f7ad42cc028ae005cb4afe3a1cf789714aa0041afb9034e`.
Public source: `e2a325d21b10141ab0071eb5c51fe2d195098df7` on talys.
The [pilot contract](../../tests/crpropa/RARE_TRANSPORT_PILOTS.md) defines the
unchanged physical scope, proposals, observations and recovery rules.

## Outcome

The external campaign completed, and its data and summary calculations are
internally consistent. The current topology-guided full-Fe56 estimator is
not ready for a final rare-channel freeze. More observed secondary particles
did not provide enough independent weighted contributions. Some apparent
variance/cost gains are misleading because a small sample misses important
contributions. This is a failed efficiency/readiness gate, not a formal
rejection of the production physics or proof of an adapter bias.

No final physical acceptance is declared. The 99% simultaneous-confidence
requirement, 1% integral tolerance, 5% spectral tolerance, fixed independent
normalizations and approved history reservations are unchanged. These pilot
observations must not be pooled with future final samples.

## Integrity, inventory and controls

- Archive hash and all 3096 manifest-listed payload hashes match.
- Source, protocol and adapter-helper fingerprints match the published code;
  the source patch hash is that of an empty patch. CRPropa is pinned at
  `73fa865b3304113da2059755a83f5939c1039d68`.
- All 576 planned jobs have matching start and completion records, complete
  output footers and the planned settings. No unmarked, incomplete or
  recovered jobs occur.
- All 288000 primary executions were completed. There are 96000 primary
  executions at the main numerical level; the other 192000 are paired
  refinements with the same 192 batch seeds. Refinements are not additional
  independent replicas. The six profiles have independent seeds.
- The 48 seeded identity-control executions reproduce original product
  identities, weights, energies and order in all 12 paired controls.
- All 1024 finite-path Mn54 control executions pass their existing engineering
  check. Maximum absolute discrepancy in empirical SE units over eight
  lepton moments: mcray 0.753/1.195 and CRPropa 0.883/0.873 at eta 0.5/0.8.
  This check is not the final statistical method-check family.
- The complete reservation therefore used 289072 primary executions.
  Actual additional root-PD continuations: 329334, below the 720000 bound.
- All 155 Linux CTest entries pass. Their logged duration is 213.48 seconds.

The analysis re-read every raw pilot TSV, checked start/output hashes and
recomputed every stored batch summary and all 72 group summaries. The largest
relative floating-point discrepancy was 3.14e-16. A separate row aggregation
reconstructed per-primary contributions for tail and paired-refinement
diagnostics. Numerical analysis artifacts and raw data are retained outside
Git in the workspace's ignored validation area.

The logged run took 52 min 46 s before final archiving, including 48 min 29 s
in the workload phase, on 16 worker processes. Summed main-job elapsed time
was 42374.93 process-seconds, or 11.77 process-hours; this is not a CPU-time
measurement. Controls, reference compilation and other workload overhead
are not included in that sum.

## Electron and positron statistics

Each profile at each background/program/level contains 4000 primary histories
in eight batches. Define contribution concentration as
`N_contrib = (sum Y_h)^2 / sum(Y_h^2)`, where `Y_h` is the entire weighted
contribution from primary h. This is a concentration diagnostic, not a
replacement for primary count or a calibrated confidence interval.

For the beta-minus guide at eta=0.8, main-level electron counts are:

| Background | Program | Nonzero primary histories | N_contrib | Largest primary / total electron yield |
|---|---|---:|---:|---:|
| CMB | mcray | 5 | 3.23 | 49.6% |
| CMB | CRPropa | 4 | 4.00 | 25.0% |
| CMB+Inoue12 | mcray | 5 | 1.13 | 94.2% |
| CMB+Inoue12 | CRPropa | 8 | 1.32 | 86.7% |

At eta=0.5 the beta-minus guide has only 1--7 nonzero electron histories
per main-level cell, with N_contrib between 1 and 2. Across all main-level
profiles the largest electron N_contrib is 4. The finest-level beta-minus
controls remain similarly concentrated: N_contrib is at most 3.24. None of
the four ordinary main-level cells observes an electron. Zero observed
variance is unresolved, not evidence of an accurately vanishing yield.

For beta-plus, more observed histories do not imply more useful information.
At the main level on CMB, mcray ordinary transport has 17 nonzero positron
histories and N_contrib=17. The eta=0.5 plus guide has 296 nonzero histories,
but N_contrib=7.68. A small fraction of histories carries most of its total
weighted yield. Across the four main-level cells, the eta=0.5 plus guide
has N_contrib=7.68--19.08 despite 296--389 nonzero histories; eta=0.8 has
N_contrib=2.34--31.77 despite 400--571 nonzero histories.

The fixed mixture changes conditional channel probabilities after a physical
event occurs. It does not increase the chance that a slow beta decay happens
within the finite path, nor does a topology-only guide account for the full
probability and remaining-path value of producing a particular parent.
Those limitations remain possible causes of poor efficiency. The present
outputs do not decompose the loss of efficiency by individual nuclear channel.

## Cost and apparent gains

At the main level the native weighted profiles cost 2.19--2.91 times the
ordinary batch cost. Mean native ordinary batch times are 28.95 s on CMB
and 35.31 s on CMB+Inoue12; corresponding weighted times range from
67.11 to 98.35 s. Reference batches take approximately 0.89--2.07 s.
Short reference timings and eight batches do not support precise throughput
claims; distinct retained and new reference builds also remain explicit.

Root splitting alone is the most promising of the tested simple controls
for positrons, without establishing a final choice:

| Program/background | Observed variance/cost gain for positron count, split vs ordinary |
|---|---:|
| mcray / CMB | 1.37 |
| mcray / CMB+Inoue12 | 2.05 |
| CRPropa / CMB | 0.49 |
| CRPropa / CMB+Inoue12 | 5.80 |

These are unstable pilot estimates, not accepted efficiency factors. There
is no electron gain estimate against an ordinary sample with zero observed
variance. There is no uniformly successful topology-guide configuration.

In particular, the report's apparent positron-count gain of 25.65 for the
mcray CMB+Inoue12 minus-strong profile must not be used for a final budget.
Its sampled mean is 0.00025233, whereas its own deterministic coefficient
reference is 0.00517750, about 20.5 times larger. This is a diagnostic of
inadequate sampling and/or a need for further estimator audit, not a physical
acceptance comparison. A small empirical variance in an unrepresentative
sample is not a demonstrated variance reduction. Large deficits cannot be
silently dismissed as normal fluctuations or treated as proven estimator bias.

## Normalization and numerical refinement

Ordinary and split-only samples have exactly 56 weighted baryons and 26
units of charge per primary. For importance sampling, those sums need only
be conserved in expectation. The native weighted mean normalization residuals
are below 2 empirical SE throughout the matrix. A reference strong-guide
outlier requires explicit follow-up: on CMB+Inoue12 at level 2, minus-strong
has mean baryons 54.96366 +/- 0.22479 SE, a 4.61-SE deficit; charge has a
corresponding 4.61-SE deficit. This is a diagnostic alarm for weight-law and
tail checks, not a calibrated rejection after all comparisons.

Native ordinary and split-only per-primary observable contributions are
identical across all three background-integration settings on both backgrounds.
For the moderate guides, the largest paired mean change is below 6e-6 of
the fixed independent normalization. This supports numerical stability in
these pilots, but does not itself establish a general error bound.

One native CMB minus-strong history changes its continuation under refinement,
causing a 16.36% normalized proton-energy shift with 16.36% paired SE. Tiny
coefficient changes can change a discrete draw and its later RNG sequence;
pathwise agreement is not a requirement of convergent Monte Carlo. This
outlier does not establish a systematic integration error.

The reference step controls are not precise enough to certify the final
1% gate. For example, ordinary CMB photon counts change from 0.45875 to
0.44875 to 0.41200 per primary. The main-to-finest paired difference is
-0.04675 +/- 0.00846 empirical SE, or -11.17% +/- 2.02% of the fixed independent
normalization. The photon-energy shift is -9.46% +/- 2.96% in that normalization.
This signal must remain open: it cannot be described as established 1%
convergence. It also does not by itself prove a systematic step bias, given
eight batches, many comparisons and changed random trajectories. Other
profiles do not reproduce a common monotonic photon trend. Refinements must
be analyzed with their actual primary/batch pairing, not counted as independent
additional statistics.

## Remaining work before a final rare-method freeze

1. Audit the complete weighted genealogy and prefix accounting using small
   independently solvable PD-plus-competing-decay networks. Include moments,
   normalization and contribution tails; the current single-parent checks
   are necessary but insufficient.
2. Investigate the observed concentration and finite-path rarity. A candidate
   next method is conditional integration of rare decay contributions using
   original hazards, with explicit PD competition and daughter continuation.
   This is a proposal for review, not an approved implementation or a change
   to the physical model.
3. Resolve reference-step diagnostics and establish a valid statistical
   interval/tail-control method before freezing proposal parameters and final
   sample budgets. Do not extrapolate budgets from the misleading gains above.
4. Keep the approved ordinary long calculation separate. This campaign does
   not launch or replace it. Preserve all pilot observations and reservations;
   do not rerun the fixed campaign or pool these samples into final acceptance.
