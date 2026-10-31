# Historical observer and ballistic stopping comparison

This deterministic prerequisite follows the merged field/transport comparison.
It is not the proposed source ensemble or an interacting cascade comparison.
[The protocol](../../tests/fork_compat/OBSERVER.md) defines extraction and the
210 fixed controls. [Machine-readable evidence](fork_observer_comparison.json)
records source/reference/input/adapter/library hashes and execution artifacts.

## Scope and implementation

The historical reference is the public snapshot
`e1fc71f613efe7fecf3210e15339c34814245cd4`. Seventeen file hashes are checked
before extraction. Original Observer membership/step methods, free motion and
Result continuation/time-left functions execute verbatim in the reference
translation unit. Its interaction-free constant-rate driver disables filters,
abort and deflection; it is not a build of the complete old engine.

Modern ballistic controls run the actual Result/PropagationEngine in both
constant- and variable-rate modes, capturing the one final particle and
checking energy, weight and identity. Common photons start from initialized
states; H0=70 km/s/Mpc and OmegaLambda=0.7 are explicit. The audited shared
CosmoTime class and units are reused; no cosmology-equivalence claim is made.

There are 174 point/step queries and 36 ballistic cases. Radii 0.01/0.1 Mpc,
origin/translated centers, strict inside/outside, birth hits, surface directions,
misses, short end-time remainders and no observer are covered. Adjacent-value,
diagonal near-surface and huge-radius cases are numeric controls, not source
samples. No statistical confidence interval or equivalence margin is chosen.

## Results

On 2026-10-31 all 210 controls completed in each environment:

| Environment | Membership differences in point queries | Old step > remaining | Ballistic hit-flag differences | Modern rate-path bitwise differences |
| --- | ---: | ---: | ---: | ---: |
| macOS arm64 Release | 8 | 56 | 0/36 | 11/36 |
| Emulated Ubuntu 22.04 Linux amd64 Release | 32 | 56 | 0/36 | 16/36 |
| Same Linux build with ASan+UBSan and leak detection | 32 | 56 | 0/36 | 16/36 |

No sanitizer diagnostics occurred. Linux runs used the frozen offline image
recorded in the evidence under amd64 emulation on an arm64 workstation; they
are not native Linux performance measurements. Ten self-contained Python
adapter tests passed and are included by the existing Release CI entry point.
The reported final adapter revision is `5d5006c`, with clean working-tree state.
Raw results, inputs, extracted code and logs remain under ignored
`build/checks/observer-fixed-*`; pre-guard artifacts remain separately retained.

### Membership precision and range

Normal center/interior/exterior and axial-surface controls agree. Differences
are confined to the numerical boundary/range controls:

- The fork narrows coordinate offsets to double and compares a sum of squares
  against double R squared. Modern code keeps the wider offsets and uses hypot.
  Linux therefore distinguishes adjacent long-double points that the old code
  rounds back onto the surface. Those account for 24 differences; macOS long
  double has double precision and does not expose this distinction.
- Six diagonal cases differ near the represented sphere surface. Decimal
  direction components and translated-coordinate arithmetic do not place all
  represented points exactly on the mathematical surface. A different result
  in this set is a rounding sensitivity, not a change from strict inequality.
- Two deliberately huge-radius controls expose overflow in the old squared
  distance/radius arithmetic. Modern hypot remains finite and classifies these
  points correctly without reproducing the old overflow.

These controls do not justify restoring lower precision or squared overflow.
The previously accepted strict endpoint-membership rule is preserved.

### Step floor and stopping time

The old observer can return its radius-based minimum step even when remaining
time is smaller. There are 56 such raw point-kernel queries; some points would
already stop at birth, so that count is not a count of realized bad trajectories.
Modern LimitStep never exceeds remaining time in the matrix.

The old ballistic driver actually overshoots the specified end time in the
short-remainder controls (three macOS cases, four Linux cases), by up to
9.9e-5 Mpc of elapsed c*time, about 99 pc. This is the known step-cap correction
already made in the modern observer, not a new migration mismatch to restore.
All 36 final hit flags agree despite these possible timing differences.

The two modern rate paths are not bitwise identical: variable-rate propagation
uses two half-moves and double step selection. Its final time is one ulp above
the represented endpoint in one macOS case and two Linux cases. The excess is
7.55e-15 Mpc on macOS or 3.69e-18 Mpc on Linux. Constant-rate time does not exceed
the represented endpoint in these controls. This finite-precision residual is
recorded explicitly and is distinct from a minimum-step overshoot. No blanket
exact-time or whole-engine equivalence is asserted.

## Confirmed no-progress bug and correction

A separate bounded probe found a real pre-existing modern engine defect:
when a nonzero end-time remainder is one ulp of the current absolute time,
both variable-rate half-moves can round back to that same time. Position changes
but the outer loop never finishes. A two-second watchdog reproduced the stall
on main `473b99e`; an adjacent initial time completed, demonstrating the
tie-to-even dependency.

Production commit `5d5006c` records time before the complete variable-rate move
and rejects a history if time has not increased after both halves. It throws
the existing fatal numerical-error type; the existing Run handler increments
the failure count and discards the entire failed history. It does not change
the half-step algorithm, invent a time increment, resample the history or
silently emit a successful partial result.

The new C++ regression constructs the one-ulp remainder and uses a bounded
logger watchdog so old code fails promptly instead of hanging CI. Against
the previous engine it reaches that watchdog and fails. Corrected code rejects
the history before the watchdog and captures no successful particle. A complete
native Release suite passed 244/244 after the fix.

Complete rebuilt suites passed in all three environments:

| Environment | CTest passed | Elapsed |
| --- | ---: | ---: |
| macOS arm64 Release | 244/244 | 16.78 s |
| Emulated Linux amd64 Release | 244/244 | 119.89 s |
| Emulated Linux amd64 ASan+UBSan, leak detection | 244/244 | 174.92 s |

All 210 original controls were repeated after adding the guard with identical
input hashes. Their raw output is byte-identical before/after in each environment;
the guard changes only the unresolved-progress failure path. The evidence
retains both the timeout-control and bounded-regression log hashes.

## Source comparison prerequisite and next scope

Source inspection also found a setup distinction to respect in a future
ensemble: the old CRbeam applies direction customization only when jet angle
is greater than 1e-14. With zero opening it leaves the default +z direction,
even if rotation parameters are supplied. Modern pencil explicitly accepts
an arbitrary axis. The shared pencil baseline must therefore use +z; rotated
comparisons should use an ordinary nonzero cone rather than assume rotated
zero-angle behavior is shared. This was a source audit, not a source-sampling
execution in this diagnostic.

The next statistical stage still needs a frozen observable family, desired
resolution, confidence allocation, seeds and sample budget. These deterministic
controls do not provide that agreement. Source/field coupling, conditional
arrival distributions, observer ensembles and matched EM cascades remain
subsequent work.
