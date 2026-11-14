# Electron/photon source in a tabulated field: cost pilot

The user approved proceeding with a cost pilot on 2026-11-01 and revised the
species/energy to **electrons and photons at 1 PeV**. Other proposed settings
remain unchanged. This separate pilot measures cost and failures only; it
provides no confidence intervals or statistical-equivalence decision and must
not be pooled into a later acceptance ensemble.

## Frozen pilot

Reuse the three source profiles, geometry, cosmology and tabulated field from
[the coupled diagnostic](COUPLED.md). Source position is the origin, initial
weight 1 and redshift 0.01, H0=70 km/s/Mpc and OmegaLambda=0.7. The strict
inside sphere has radius 0.01 Mpc, its center is 0.1 Mpc along the declared
axis, and the flight horizon is 0.2 Mpc. Use the existing divergence-free
3x4x5 transport grid with zero extension and Deflection3D accuracy=20.
There are no random/CEL interactions, energy losses, thinning or filters;
photon deflection is absent because its electric charge is zero. Run the
constant-rate path, preserving each implementation's endpoint/step arithmetic.

| Profile | Species | Half-angle, degrees | Axis theta/phi, radians | Fork seed | Modern seed | Count per implementation |
| --- | --- | ---: | --- | ---: | ---: | ---: |
| pencil_z | electron | 0 | 0, 0 | 2005101 | 2005102 | 1000 |
| pencil_z | photon | 0 | 0, 0 | 2005103 | 2005104 | 1000 |
| cone_z | electron | 20 | 0, 0 | 2005105 | 2005106 | 1000 |
| cone_z | photon | 20 | 0, 0 | 2005107 | 2005108 | 1000 |
| cone_rotated | electron | 20 | 0.4, 0.7 | 2005109 | 2005110 | 1000 |
| cone_rotated | photon | 20 | 0.4, 0.7 | 2005111 | 2005112 | 1000 |

There are **12000 complete primary histories**. These twelve reservations
are separate from the spent photon pilot/final ensemble seeds. One uninterrupted
source stream serves each case/implementation; do not reseed every primary or
interpret equal indices as matched histories. No vMF/Gaussian source, source
energy law or cascade is tested. The modern propagator uses seed 17 and
consumes no random interaction draws in this model.

A separate two-primary-per-stream software fixture (24 histories total) uses
seeds 2006101–2006112. It is the default mode and is never pilot evidence.
Fixtures may be replayed for software checks; do not rerun spent pilot seeds,
replace failed histories or adapt the count after viewing results.

## Implementation and failure contract

Reuse the hash-pinned original cone initializer and RNG, plus actual historical
field/deflector/move/observer functions from the coupled adapter. Modern paths
use actual SourceGeometry and the shared modern coupled engine implementation.
No production function or physical model is changed. The new common adapter
files split the existing modern driver from its fixed-state comparison main;
this avoids a second independently maintained implementation.

Load one immutable grid per implementation before trajectory timing, then reuse
it sequentially across source streams. A fresh counting wrapper, deflector,
result and engine are built per primary. Field values, source directions and
weights are not renormalized to manufacture agreement. Energy and identity are
checked on both propagation paths; modern output weight is checked too. The
minimal old particle shell has no weight field: the actual source weight is
retained separately and there is no weight-changing propagation operation.

A 100000-engine-move per-history guard, compiler timeout 240 seconds and
execution timeout at most 600 seconds bound the work. Compilation/execution/
reduction share a **900-second whole-run budget**. A failed process, engine
failure, stalled/over-budget trajectory, missing/duplicate/nonfinite record or
deadline overrun invalidates the entire pilot. Preserve partial raw data/logs
and failure.json; publish report.json only after complete validation. There
are no retries or replacement histories. Pilot collection requires a committed
clean source revision and a new/empty output directory.

## Run and artifacts

Use the pinned public snapshot and configured build dependencies described in
AGENTS.md and the coupled protocol. Python standard library is sufficient.
From this public repository root:

```sh
mkdir -p build/checks/field-source-tmp
TMPDIR="$PWD/build/checks/field-source-tmp" PYTHONDONTWRITEBYTECODE=1 \
  python3 tests/fork_compat/field_source_pilot.py \
  --fork-root /path/to/pinned/reference --build-dir build/Release \
  --output build/checks/field-source-fixture
```

Use explicit **--mode pilot** only for the authorized new pilot reservation;
its output directory must differ from the fixture. --command-prefix supplies
a shell-free prepared offline compiler/runtime prefix for Linux checks. Any
container timeout must cover its selected compiler/probe phase. Keep all caches,
outputs and temporary files in the permitted workspace.

Before execution, preserve source/adapter/reference/core hashes, exact Git
revision and dirty state, generated design/grid hashes and compiler command/
version/logs. Raw 21-digit JSONL retains primary IDs, source weights/directions,
final position/direction, actual flight parameter, hit/miss, signed end-time
error, engine-move and field/scale counters. A successful report requires all
12000 pilot histories and twelve source-stream timing records.

Per-history compute timers include source sampling, engine setup, propagation
and diagnostic counters; they exclude grid parsing, source-stream setup, JSON,
compilation and startup. Whole-process timing includes output/startup. Historical
propagation uses the reduced old driver; modern uses the full current engine.
These costs do not benchmark identical driver overhead and must not be quoted
as native Linux performance when run under emulation.

The reducer reuses complete-primary accounting from the photon cost pilot and
adds field/engine counters. Descriptive hit counts and cost are not acceptance
criteria. The next statistical protocol still needs agreed observables,
confidence/precision choices, final counts, new seeds and stopping rules.
Checks expected to exceed 30 minutes belong on an external Linux amd64 server.

See the [confirmed pilot results](../../docs/validation/fork_field_source_pilot.md)
for the completed reservation, measured costs, known stopping differences
and limits. Those seeds are spent; do not turn the run example into a rerun.

The separately approved [1 EeV three-species pilot](EEV_FIELD_SOURCE_PILOT.md)
uses explicit --design eev with new reservations. The default PeV design
and its completed evidence remain separate.
