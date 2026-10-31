# Source/observer cost pilot

This is a separate, frozen **cost and failure diagnostic**, approved on
2026-10-31. It is not the final statistical comparison, provides no acceptance
verdict and must not be pooled into a later acceptance ensemble.

## Frozen design

Use photons at 1e12 eV, initial weight 1, source redshift 0.01,
H0=70 km/s/Mpc and OmegaLambda=0.7, no interactions or magnetic field.
The source is at the origin. An inside sphere with radius 0.01 Mpc is
centered 0.1 Mpc along the declared source axis. Both implementations use
the same represented center, strict endpoint membership, a 0.2 Mpc flight
budget and the constant-rate path. Stop at the first endpoint satisfying
the respective observer, or the final time; this preserves the previously
documented historical step-floor overshoots.

| Profile | Half-angle | Axis angles (theta, phi), radians | Fork seed | Modern seed | Primaries per implementation |
| --- | --- | --- | --- | --- | --- |
| pencil_z | 0 degrees | 0, 0 | 2003101 | 2003102 | 1000 |
| cone_z | 20 degrees | 0, 0 | 2003103 | 2003104 | 1000 |
| cone_rotated | 20 degrees | 0.4, 0.7 | 2003105 | 2003106 | 1000 |

Reserve these six source seeds for this pilot only. Each case maintains one
RNG stream per implementation; do not reseed each primary or require matching
trajectories. The modern propagation engine uses seed 17 and consumes no
interaction random draws in this model. No retries, replacements or adaptive
history counts are allowed. A failed process, failed modern history, reference
step-budget failure or missing/duplicate/nonfinite record invalidates the
entire pilot. Preserve partial logs, without a successful report.

The old cone initializer, RNG, free photon movement, observer membership and
step patching, continuation and time-left methods are extracted verbatim from
the hash-pinned historical snapshot. The historical driver is the reduced
interaction-free loop described in [OBSERVER.md](OBSERVER.md), rather than
the entire old executable. The modern path uses SourceGeometry and the actual
Result/PropagationEngine. Source weights are retained; no renormalization or
correction of historical directions is applied. Historical ballistic motion
does not modify weights; modern final energy/type/weight are checked against
their inputs. Rotated zero-angle pencils, translations, outside observers,
vMF/Gaussian sources and interacting cascades are outside this pilot.

## Run and retained output

Use the dependencies, public snapshot retrieval and configured build described
in [OBSERVER.md](OBSERVER.md). Normal builds and CI do not need the snapshot.
The output directory must be new or empty. From the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 tests/fork_compat/source_observer_pilot.py \
  --fork-root /path/to/pinned/reference --build-dir build/Release \
  --output build/checks/source-observer-pilot
```

`--command-prefix` accepts a shell-free compiler/program argv prefix, for
example a prepared offline Linux amd64 container. The program has a 180-second
process timeout; the reference driver also has a 100000-step per-history cap.
Keep temporary/cache directories inside the permitted workspace.

Before executing, the tool records the fixed design, reference hashes, source
and adapter hashes, Git revision/dirty status, linked core hash and compiler
command/version. Raw JSONL records every emitted direction/weight, primary ID,
hit/miss, elapsed time and endpoint, signed end-time error and historical step
count. The report requires all 6000 primary records and six timing records.
It contains descriptive totals and failure diagnostics, with no intervals,
significance tests, acceptance decision or pilot-driven bin selection.

Compute timers include source sampling and propagation, exclude JSON output,
compilation and startup. A separate process wall clock includes output/startup
and, when used, container launch. Historical and modern paths have different
driver overheads; these timings do not benchmark equivalent engine complexity.
Emulated amd64 timings on an arm64 workstation are not native Linux estimates.

The accepted target for the final comparison is simultaneous 95% intervals
and absolute hit-fraction half-width 0.001. The observable family, interval
method, counts, disjoint seeds and stopping rules still need to be frozen
before collecting final data. No final ensemble is authorized by running
this pilot; checks expected to exceed 30 minutes need an external server script.
