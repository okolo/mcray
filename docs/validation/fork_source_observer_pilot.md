# Historical source/observer cost pilot

The separately approved pilot completed on 2026-10-31 at clean source revision
`0f8d158`. This is cost and failure evidence, **not statistical acceptance**.
The [frozen protocol](../../tests/fork_compat/SOURCE_OBSERVER_PILOT.md) defines
the geometry, six independent source seeds, 1000 primaries per case and
implementation, retained raw records and failure policy.
[Machine-readable evidence](fork_source_observer_pilot.json) retains the
design, reference/source/tool hashes, descriptive totals and artifact hashes.

## Model and implementation

Photons start at the origin at z=0.01, with energy 1e12 eV and weight 1.
H0=70 km/s/Mpc and OmegaLambda=0.7 are common. There are no interactions or
magnetic fields. The strict inside observer has radius 0.01 Mpc and a center
0.1 Mpc along the declared source axis. Both paths use the constant-rate
algorithm and a 0.2 Mpc flight budget. Rotated cone angles are
(theta,phi)=(0.4,0.7) radians; the zero-opening pencil uses +z.

Seventeen historical file hashes are checked against public snapshot
`e1fc71f613efe7fecf3210e15339c34814245cd4`. The reference uses the original
ordinary-cone block and RNG, free motion, observer and continuation functions.
Its driver is the reduced interaction-free loop, not the complete historical
application. Modern sampling uses SourceGeometry and propagation uses the
actual Result/PropagationEngine. Weights are retained verbatim; modern final
energy, identity and weight are checked. No source normalization correction,
resampling or replacement history is applied.

## Descriptive results

All 6000 histories completed in macOS arm64 Release, emulated Linux amd64
Release and the same Linux environment with ASan+UBSan/leak detection.
Each environment produced the following counts:

| Source | Fork hits / primaries | Modern hits / primaries |
| --- | ---: | ---: |
| +z pencil | 1000 / 1000 | 1000 / 1000 |
| Axial 20-degree cone | 76 / 1000 | 84 / 1000 |
| Rotated 20-degree cone | 82 / 1000 | 84 / 1000 |

Every emitted weight was 1, so hit weights equal hit counts. No failed
histories, nonfinite records, negative/zero weights or final-time overshoots
occurred. The largest historical per-history step count was 2110. Linux
Release and sanitizer history records are identical after parsing; timing
records are excluded from that comparison. Replaying these same pilot seeds
across toolchains does not create additional independent statistical evidence.
These counts have no acceptance intervals and are not used to select final
observables or bins.

No sanitizer diagnostic occurred. Execution logs retain a startup warning
about absent optional interaction-table directories; this interaction-free
model does not load those tables. It does not establish readiness for cascades.
Twelve self-contained Python guards passed, including incomplete/duplicate/
nonfinite pilot rejection, and use the existing Release CI entry point.
The helper refactoring also retained byte-identical output for all 210 native
[observer controls](fork_observer_comparison.md). The earlier complete C++
suites passed 244/244 in all three environments after the time-progress fix;
the pilot adds no production-code change.

## Cost measurements and limits

| Environment | Combined sampling/propagation time, 6000 histories | Process wall time |
| --- | ---: | ---: |
| macOS arm64 Release | 0.0323 s | 0.130 s |
| Emulated Linux amd64 Release | 1.410 s | 2.215 s |
| Emulated Linux amd64 ASan+UBSan | 1.576 s | 2.474 s |

Compute timers exclude compilation, startup and JSON output. Wall times
include startup/output and container launch where applicable. The historical
reduced driver and the complete modern engine have different overheads, so
these are planning measurements, not a performance comparison of equivalent
engines. Linux ran under amd64 emulation on an arm64 workstation using the
frozen offline image recorded in the earlier observer evidence. Its cost is
not a native Linux server benchmark. Millisecond totals and process launch
overheads do not justify precise extrapolations.

Raw histories, inputs, extracted sources, executable, manifests and logs remain
under ignored `build/checks/source-observer-pilot-{native,linux,sanitizers}`.
The CLI refuses to overwrite a nonempty output directory. The report and
hashes above preserve what was run, while the public protocol/adapter support
independent reproduction without an enclosing private repository.

## Next statistical step

The accepted target is simultaneous 95% intervals and absolute hit-fraction
half-width 0.001. Before collecting final data, freeze the observable family,
interval method/confidence allocation, counts and disjoint source seeds.
Do not pool this cost pilot into that ensemble or treat zero-containing
intervals as proof of arbitrarily small systematic error. Longer than
30-minute checks require a separate external Linux amd64 server script.
Field coupling and interacting cascades remain later stages.
