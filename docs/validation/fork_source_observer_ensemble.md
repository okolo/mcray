# Frozen source/observer ensemble results

On 2026-10-31 the approved 48-million-history comparison completed locally
at clean revision `d3b7611aa3b69c46e72c41bc91aa0354625472d9`.
All 33 modern-minus-fork difference intervals contain zero. The implementations
are compatible for the frozen observable family at the reported resolution;
this is not proof of exact equality or complete physical equivalence.

The [approved protocol](../../tests/fork_compat/SOURCE_OBSERVER_ENSEMBLE.md)
fixes the physical model, 39-interval simultaneous family, counts, disjoint
seeds, bounds and failure rules before collection.
[Machine-readable evidence](fork_source_observer_ensemble.json) retains the
results, intervals, source/input/reference hashes and artifact hashes.

## Model and sampling

There are 8 million independent complete primary histories per implementation
and profile: +z pencil, axial 20-degree cone, and 20-degree cone rotated by
(theta,phi)=(0.4,0.7) radians. Source seeds 2004101 through 2004106 were used
once; pilot histories/seeds are excluded. Six isolated processes maintained
one uninterrupted source stream each, without reseeding at block boundaries.

Photons start at the origin at z=0.01, with energy 1e12 eV and weight 1.
The observer is a strict inside sphere of radius 0.01 Mpc centered 0.1 Mpc
along the declared source axis; the flight budget is 0.2 Mpc. H0=70 km/s/Mpc
and OmegaLambda=0.7 are explicit. There is no magnetic field or interaction.

The historical path executes hash-pinned original source/observer/free-motion
functions in a reduced interaction-free driver. It is not the complete old
application. Modern sampling uses SourceGeometry and propagation uses the
actual Result/PropagationEngine. Seventeen historical source hashes match
public snapshot `e1fc71f613efe7fecf3210e15339c34814245cd4`.

## Statistical results

| Source | Fork hits / 8000000 | Modern hits / 8000000 | Modern minus fork hit fraction |
| --- | ---: | ---: | ---: |
| +z pencil | 8000000 | 8000000 | 0 |
| Axial 20-degree cone | 665842 | 665387 | -0.000056875 |
| Rotated 20-degree cone | 665284 | 666163 | +0.000109875 |

The simultaneous family contains 33 difference intervals (hit indicator,
four unconditioned arrival indicators and six coordinate/direction means
for each profile), plus six individual hit-fraction intervals. Two-sided
Hoeffding bounds and the union bound give at least 95% simultaneous coverage
under the protocol's bounded-independent-history assumptions. Observables
from the same history need not be independent.

All 39 intervals were computed without clipping or adjusting the frozen
sample budget. The deterministic half-widths are:

- Hit/arrival-indicator differences: **0.000958674**.
- Individual hit fractions: **0.000677885**.
- Direction-component mean differences: 0.001917347.
- Coordinate mean differences: 0.000479337 Mpc, or 0.479337 kpc.

The accepted 0.001 hit-fraction precision target is met. The largest absolute
difference relative to its interval radius is the rotated-cone arrival
indicator at 0.092 Mpc: 0.000178875 / 0.000958674 = approximately 0.187.
Every difference interval contains zero. No bin, threshold, seed, count,
confidence allocation or failure rule was changed after observing final data.

## Validity and retained numerical differences

All 48 million histories completed with no reported failed histories,
replacement histories, nonfinite/bound violations or incomplete streams.
Initial/final weights remained 1. Modern final photon identity and energy
were checked. The runner verified six expected raw-record counts, six valid
summaries, unchanged adapter/production/reference/core hashes and each gzip
artifact's CRC and uncompressed SHA-256.

The historical minimum-step behavior remains visible:

| Source | Fork histories past end time | Largest excess, expressed as length | Modern histories past end time |
| --- | ---: | ---: | ---: |
| +z pencil | 0 | 0 | 0 |
| Axial cone | 1269 | 9.998733 pc | 0 |
| Rotated cone | 1397 | 9.993342 pc | 0 |

The fork's inside-sphere step floor is 0.001 R = 10 pc and can exceed the
remaining time. The modern cap was already corrected. These are separate
numerical diagnostics, not failed histories or a new acceptance tolerance.
A separate streaming audit of the two raw fork files confirmed that all
2666 overshooting histories were misses; none was counted as an arrival.
The statistical conclusion does not erase this known algorithmic difference.

## Execution, checks and artifacts

Native macOS arm64, Clang 22.1.5, six child processes on a host reporting ten
logical CPUs: **139.858 seconds** including compilation, collection,
compression and raw-archive read-back, before the final metadata bundle.
The explicitly authorized local four-hour limit was not approached.
The six compressed history files occupy **2.799 GiB**. Every raw history is
retained; the reducer uses bounded blocks and compensated sums rather than
loading millions of dictionaries into memory.

Prerequisites passed: 244/244 native C++ tests and 16 self-contained Python
tests, including child-process deadline cancellation. Native fixture records
matched all 6000 frozen pilot records, and independently recomputed sums
matched to the documented numerical roundoff tolerance. A separate 6000-history
fixture of the new C++ accumulator passed in the frozen offline Linux amd64
ASan/UBSan environment with leak detection and no sanitizer diagnostics.
That sanitizer check replays pilot seeds; it is not the 48-million-history
final ensemble. Optional missing-interaction-table startup warnings remain
in logs; this model loads no interaction tables.

Raw data, progress logs, executable, manifest, verified summaries, full report
and compressed metadata remain under ignored
`build/checks/source-observer-ensemble-final-native`. Compiler and fixture/
regression logs are separately retained and hashed in the evidence. Normal
CI runs the Python guards without spending these final reservations.

The final ensemble is native macOS evidence. The Linux sanitizer fixture ran
under amd64 emulation and supplies no native-Linux statistical/performance
claim. Native Linux amd64 reproduction, field coupling, other observers and
interacting cascades remain separate follow-up stages. This finite family
does not test every possible joint distribution or source/transport feature.
