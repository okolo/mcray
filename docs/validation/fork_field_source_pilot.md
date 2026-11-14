# 1 PeV electron/photon field/source cost pilot

The approved pilot completed on 2026-11-01: **12000 complete primary histories,
zero failed histories/processes**, from clean public revision
`f49bd273110da5165a4c5f21a8b39d837bf29e43`.
The [frozen protocol](../../tests/fork_compat/FIELD_SOURCE_PILOT.md) owns
species, energy, settings, counts, disjoint seeds and stopping/failure rules.
[Machine-readable evidence](fork_field_source_pilot.json) retains source/input/
reference/core hashes, raw/archive hashes, timing and descriptive counts.
This is cost/failure evidence, without a statistical-equivalence verdict.

## Run and provenance

Both implementations emit electrons and photons at 1 PeV, initial weight 1,
from the origin at z=0.01. The +z pencil and two 20-degree ordinary cones
use the previous axial/rotated geometry. They propagate through the existing
3x4x5 tabulated field with zero extension, accuracy 20, common cosmology and
strict inside sphere at distance 100 kpc/radius 10 kpc, up to 0.2 Mpc flight.
There are no random/CEL interactions, energy losses or cascades.

One immutable field per implementation is loaded before timing. The actual
old cone initializer/RNG, deflector, move and observer functions execute in
the reduced historical driver. The modern path uses SourceGeometry and the
actual constant-rate engine. The shared modern propagation/output helper was
split from the fixed-state main to avoid duplicating its implementation.
Only validation-tool files changed; production algorithms are unchanged.

The native macOS arm64 process took **3.412 seconds** including output/startup.
Sampling/propagation timers sum to **3.254 seconds**. The whole run including
compilation and reduction took **10.400 seconds**, below the 900-second limit.
These timings measure this workstation/driver configuration, not native Linux
performance or equivalent engine overhead. The historical reduced driver and
modern full engine have different setup/diagnostic costs.

## Descriptive results

Every row below has 1000 primaries per implementation. Times include source
sampling, per-primary engine/deflector setup, propagation and counters; exclude
grid loading, source-stream setup, JSON output, compilation and startup.

| Profile | Species | Fork hits | Modern hits | Fork compute, ms | Modern compute, ms |
| --- | --- | ---: | ---: | ---: | ---: |
| +z pencil | electron | 0 | 0 | 611.570 | 505.730 |
| +z pencil | photon | 1000 | 1000 | 1.161 | 6.102 |
| Axial 20-degree cone | electron | 0 | 0 | 605.548 | 501.544 |
| Axial 20-degree cone | photon | 83 | 74 | 2.837 | 9.330 |
| Rotated 20-degree cone | electron | 0 | 0 | 547.812 | 448.811 |
| Rotated 20-degree cone | photon | 70 | 92 | 3.251 | 10.521 |

Cone streams use independent source seeds and different RNGs. Their unequal
counts are descriptive pilot results, not a significance/acceptance test.
The pencil directions are deterministic; repeated pencil histories do not
measure source-distribution uncertainty. All emitted/hit weights retain their
ordinary-cone meaning, with no zero or negative emitted weights. Maximum
initial direction-norm residual is 2.2205e-16.

No electron reaches the observer in either implementation in this pilot.
Hit fractions alone would therefore be an uninformative test of the sampled
electron endpoints/directions at this precision. Zero observed hits does not
prove zero hit probability for the whole cone. A later comparison needs
pre-agreed endpoint/direction observables in addition to hit/time quantities;
it must not invent acceptance criteria from these pilot outcomes.

The fork has one end-time overrun in the axial-cone photon stream, by
4.8227468791e-6 Mpc of flight parameter, approximately **4.823 pc**. Modern has
none. This is the retained historical minimum-step behavior documented in the
[observer comparison](fork_observer_comparison.md), not a replaced/failed
history. Maximum engine moves are 4 for electrons and 2116 for photons; these
count outer moves, not internal magnetic integration stages. No historical
zero-scale calls occur. Energy/identity invariants and modern output weights
pass their per-history checks.

For the identically initialized +z pencil, all 1000 records within each
implementation repeat the same state. Old/new electron position differs by
8.3526e-16 Mpc, direction by 7.1254e-14 and flight time by zero on this platform.
Photon pencil states agree exactly. This extra fixed-input diagnostic is not
a pairing of the independent cone streams or a continuum-accuracy claim.

## Software checks and retained artifacts

- Twenty-one self-contained Python guards pass and use the existing Release
  CI entry point. They cover the revised model/reservations and rejection of
  partial/duplicate/nonfinite/invalid-weight/final-vector records.
- Separate 24-history software fixtures pass on native macOS Release,
  emulated Linux amd64 Release and Linux ASan+UBSan with leak detection.
  Their seeds 2006101–2006112 are disjoint from this pilot. No sanitizer
  diagnostics occur; these small fixtures are not statistical pilot evidence.
- Repeated pencil records in all three fixtures have identical per-history
  states/counters, checking sequential reuse of the immutable field.
- After splitting the common adapter, all 729 native fixed-state JSONL records
  remain byte-identical to the frozen f53b533 comparison. Its earlier numerical
  findings are unchanged; the frozen report has not been relabelled.

The pilot uses seeds **2005101–2005112**, now spent. Do not rerun or reuse them
in a final ensemble. Inputs, original extracted fragments, compiler command/
version/logs, manifest, executable, raw JSONL and compressed archive remain
under ignored build/checks/field-source-cost-pilot-native. The gzip archive
contains 7965566 raw bytes compressed to 933225 bytes. Full decompression
reproduces SHA-256
`0009799293a202e1748be6b259f64e8d224fa55c7649276b78cf8e80e2daa6be`.
Startup missing-interaction-table warnings are retained; this interaction-free
probe does not read those tables. The linked evidence records the compiler,
source/adapter/core hashes, execution prefixes and Linux image pin.

## Cost implications and remaining decisions

A linear compute-only extrapolation on this workstation gives roughly
7.5–10.2 minutes per million electron histories per stream, and 1.2–10.6 seconds
per million photons per stream. Repeating this entire twelve-stream matrix
at one million histories per stream would take about 54.2 minutes serially,
excluding output/startup and without a scaling guarantee. This is a cost
estimate, not a selected final sample size or server benchmark. Longer work
belongs on an external Linux amd64 server under the project rules.

The next statistical stage still needs agreed observable bounds/family,
confidence method and target precision, final counts, new independent seeds,
numerical refinement and stopping/failure rules. This pilot does not establish
physical equivalence, validate interacting cascades or approve a new ensemble.
The original photon final comparison and vMF point controls retain their
separate frozen scope.
