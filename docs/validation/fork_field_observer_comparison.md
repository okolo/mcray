# Historical field/observer/engine coupling diagnostic

On 2026-11-01, all **729 fixed-state cases** completed on native macOS arm64,
emulated Ubuntu 22.04 Linux amd64 Release and the same Linux build with
ASan+UBSan/leak detection. Hit flags agree in every old/new pair, including
216 tangent controls per environment. This is a deterministic coupling check,
not a sampled-source ensemble or a physical-equivalence verdict.

[The protocol](../../tests/fork_compat/COUPLED.md) owns the frozen matrix and
reproduction command. [Machine-readable evidence](fork_field_observer_comparison.json)
records exact source/reference/input/library hashes, compiler versions,
execution prefixes, raw-artifact hashes, per-field errors and refinements.
The adapter revision is clean `f53b53351fdbdbc38d2280142c881afac978b50a`, based on
merged PR #17 (`e82f8493fbe54cc21f852ad6b8ee80efb1765db3`). Production code was
not changed for this diagnostic. The linked cores have the same production
sources as that base; their distinct configuration hashes are retained.

## What executes

The read-only reference is snapshot
`e1fc71f613efe7fecf3210e15339c34814245cd4`. Seventeen original-file hashes
are verified before extraction. The existing actual old field, Deflection3D,
particle and observer functions are combined with the verbatim original
PropagationEngine::move body. The reduced outer driver follows the old
constant-rate path with empty interaction lists, no filters and no abort.
It is not a build of the complete historical application. Modern runs use
the actual constant-rate Result/PropagationEngine, SphereObserver and
Deflection3D, capturing one final particle and rejecting failures.

The 27 common initial directions cover the +z pencil and axial/rotated
20-degree cone domains, with axis, interior, tangent and near-rim rays.
They are geometric states, without source RNG or vMF weights. Cross them
with electrons, positrons and protons at 1 EeV, no/constant/tabulated fields
and accuracy 10/20/40. The sphere has radius 10 kpc and its center is 100 kpc
along the profile axis; initial z=0.01, duration=0.2 Mpc, H0=70 and
OmegaLambda=0.7. There are no random/CEL interactions or energy losses.
The grid is exactly the existing isolated transport fixture.

## Old/new differences

| Environment | Hit differences | Maximum position difference, Mpc | Maximum direction difference | Maximum flight difference, Mpc |
| --- | ---: | ---: | ---: | ---: |
| macOS arm64 Release | 0/729 | 5.3919e-16 | 7.2218e-16 | 0 |
| Linux amd64 Release, emulated | 0/729 | 1.0000e-5 | 2.2725e-17 | 1.0000e-5 |
| Linux ASan+UBSan, emulated | 0/729 | 1.0000e-5 | 2.2725e-17 | 1.0000e-5 |

Each implementation records 360 hits. There are no end-time overruns or old
zero-scale calls in this matrix; the largest engine move count is 2097.
No sanitizer diagnostics occurred. Startup missing-interaction-table warnings
are retained: this interaction-free probe does not read those tables.

### Linux one-step boundary amplification

The Linux maximum is **10 pc**, confined to 18 no-field axial cases:
+z pencil and the axial-cone axis, crossed with three species and three
accuracy values. The old path stops after one engine move; modern needs two.
Both detect the same sphere. All other engine-move counts agree. Among those
other 711 cases the largest position difference is 2.2976e-17 Mpc and flight
difference 1.4746e-17 Mpc. A total of 176 Linux rows have nonzero time
subtractions, including the 18 one-step cases and much smaller residuals;
the evidence retains these without applying a pass/fail cutoff.

The first step is the distance to the axial surface. The old observer computes
it in double; the modern observer retains coord_type (long double on Linux):

```cpp
double center = 0.1 * units.Mpc, radius = 0.01 * units.Mpc;
double oldStep = center - radius;
coord_type modernStep = coord_type(center) - radius;
```

A separate bounded Linux first-step control confirms:

- oldStep / Mpc = 0.090000000000000016517;
- modernStep / Mpc = 0.0900000000000000109876;
- their difference is 5.5296429542e-18 Mpc;
- the modern first point has distance-radius exactly zero: both membership
  functions reject it under the strict inside rule;
- the old first point is inside by 5.5296429542e-18 Mpc: both membership
  functions accept that point.

Thus, for these cases the difference originates in step arithmetic rather
than a disagreement between membership functions at the same point. Starting
exactly on the surface, modern takes the shared minimum inside-observer step
0.001*R=10 pc. This amplifies a last-bit step difference into a different
stopping endpoint/time. macOS long double has double precision and does not
expose this particular distinction. The first-step control source, output and
compiler-command hashes are retained in the evidence.

The earlier [observer report](fork_observer_comparison.md) documents the wider
coordinate precision and overflow-safe membership differences separately.
This diagnostic does not lower modern precision, change strict membership,
remove the historical step floor or add a new crossing solver to force matching
endpoints. No new production correction is justified by these rows alone.

## Refinement within each implementation

All 486 accuracy-pair refinements per environment are retained, separately for
old and modern implementations. No hit flag changes between 10/20 or 20/40 in
this matrix. Maximum position changes in the modern path are:

| Field | Accuracy pair | macOS, Mpc | Linux Release / sanitizers, Mpc |
| --- | --- | ---: | ---: |
| No deflector | 10→20 and 20→40 | 0 | 0 |
| Constant | 10→20 | 1.0222e-5 | 1.0222e-5 |
| Constant | 20→40 | 3.6290e-6 | 3.6290e-6 |
| Tabulated | 10→20 | 3.3711e-5 | 3.3711e-5 |
| Tabulated | 20→40 | 1.7986e-5 | 1.7514e-5 |

Old refinements follow the same scale; exact values, direction/time changes
and hit-change flags are in the evidence and raw reports. Integration and
observer stopping are coupled: unchanged hit flags do not imply identical
stopping times. These maxima are diagnostics, not estimated continuum errors
or a universal convergence rate. Boundary and platform sensitivities remain
visible instead of selecting only favorable trajectories.

## Reproduction, checks and limits

Use the protocol command with the pinned public snapshot and a configured
current build. Raw 21-digit states, extracted fragments, inputs and logs remain
under ignored build/checks/coupled-final-{native,linux-release,linux-sanitizers}.
The supplemental first-step check is under build/checks/coupled-boundary-linux.
All three final probes ran from the clean adapter revision above; results/docs
were added afterward. Development runs are retained separately.

Nineteen self-contained Python tests passed, including fixed input/grid
coverage, duplicate/nonfinite/incomplete rejection, retention of hit differences
and decimal refinement below double precision. They run through the existing
Release CI entry point. Each final C++ probe is separately compiled and
executes all 729 cases. No production source changed, so the already completed
253-test PR #17 suites were not rerun for documentation/test-tool additions.
Linux runs use the pinned offline image in the evidence under amd64 emulation;
they are not native Linux performance measurements.

The reduced historical driver does not validate the old engine's interaction,
filter, retry, threading or variable-rate behavior. No field exits or zero-scale
shortcut cases are reached here; those remain covered by the isolated
transport diagnostic. Energy/identity are checked in both coupled paths;
weight is checked in modern output only, and is absent from the minimal old
particle shell because this fixed-state path has no weight-changing operation.

A real-source field/observer pilot and statistical ensemble need separate
frozen histories, seeds, observables and precision/confidence decisions. The
completed photon ensemble and restored vMF point controls keep their original
scope; these trajectories do not expand either into field/cascade acceptance.
