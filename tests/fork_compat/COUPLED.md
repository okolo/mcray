# Fixed-state field/observer coupling diagnostic

This follows merged PR #17. It implements the requested next coupling step
using previously approved physical models and numerical controls, without
choosing new statistical acceptance criteria. It reports every numerical
and hit/stopping difference; it does not declare physical equivalence.

## Scope and frozen inputs

Reuse the isolated transport energy 1 EeV for electrons, positrons and protons,
accuracy values 10/20/40 and its exact divergence-free 3x4x5 tabulated field.
The old transport input generator supplies the grid rather than maintaining
a second prescription. Auxiliary generated transport cases are not executed.
The grid spans x=[-0.2,0.2], y=[-0.3,0.3], z=[-0.4,0.4] Mpc. Other fields
are no deflector and the existing oblique (0.2,0.3,0.5) nG constant field
with positive scale 0.1 Mpc. Grid extension is zero; the baseline starts at
the origin inside the grid. Old zero-scale calls remain visible if encountered.

Reuse the source/observer geometry: origin source, +z pencil, axial and
rotated 20-degree cone axes (theta,phi)=(0.4,0.7), strict inside sphere with
R=0.01 Mpc at distance 0.1 Mpc along the declared axis, horizon 0.2 Mpc.
Common cosmology is H0=70 km/s/Mpc, OmegaLambda=0.7 and initial z=0.01.
There are no random/CEL interactions, adiabatic losses, filters or aborts.
Use the constant-rate engine path only.

The driver supplies identical fixed initial particle states to both paths.
No source sampler, random ensemble or vMF weighting is invoked. Directions
are geometric probes of the cone domain, not a substitute source distribution:

- One +z pencil direction.
- Thirteen directions per cone: its axis, plus four azimuths at each of
  half the geometric tangent angle asin(R/d), that tangent angle, and a
  near-rim offset with (1-cos(alpha))/(1-cos(20 degrees))=0.99.
- Initial directions use the fork's two rotation formulas and are normalized
  once before serialization. Observer centers are shared represented inputs.
- The 27 states crossed with three species, three fields and three accuracy
  values give **729 cases**, frozen and hashed before trajectory evaluation.

Tangency states are floating-point/grazing controls. Keep their outcomes
separate from interior/near-rim states. Strict endpoint sampling, its minimum
step and the already documented old remaining-time overrun are retained;
this does not replace them with a crossing solver. No numerical tolerance
for an equivalence verdict is introduced.

## Original functions and modern execution

The same seventeen reference pins as observer_reference.json are verified
against immutable snapshot e1fc71f613efe7fecf3210e15339c34814245cd4. Reuse the
existing particle/time/unit/field/Deflection3D/observer/continuation extraction.
Add the **verbatim original PropagationEngine::move body** and a small storage/
logger shell with one actual deflector. The outer driver is the exact
constant-rate reduction when interaction/filter/abort lists are disabled:
continue, get time left, patch the observer step, then call original move().
No old library or production physics function is replaced by a lookalike.

This is a reduced reference driver, not a complete old-engine build. Trace
pointers are null. Counting loggers observe engine moves and impose a
100000-move diagnostic budget. Source RNG is absent; reference/modern
propagator RNG objects are initialized but deflection consumes no draws.
A stalled/over-budget/incomplete/failed trajectory aborts the report without
replacement. Modern execution uses the actual Result/PropagationEngine,
SphereObserver and Deflection3D with captured final output. Energy and identity invariants are checked in both paths; modern output
weight is checked too. The minimal historical shell has no weight field
because this fixed-state path has no weight-changing operation. Current APIs and ownership are preserved;
no production implementation is changed for this diagnostic.

## Output, refinement and reproduction

```sh
mkdir -p build/checks/coupled-tmp
TMPDIR="$PWD/build/checks/coupled-tmp" PYTHONDONTWRITEBYTECODE=1 \
  python3 tests/fork_compat/coupled.py --fork-root /path/to/pinned/reference \
  --build-dir build/Release --output build/checks/coupled-run
```

Use the configured build/link dependencies and public reference retrieval in
AGENTS.md. Python standard library is sufficient; --command-prefix supports
the prepared offline Linux container and sanitizer build. Compiler timeout is
240 seconds and execution timeout is 180 seconds. Output must be new/empty,
outside the read-only reference and within the permitted writable workspace.

Preserve raw 21-digit JSONL states, case/grid design and hashes, original
fragments, compiler command/version/logs, source/adapter/core hashes and exact
Git revision/dirty status. Raw old/new state records include hit/past-end,
position, direction, actual flight parameter, beta/mass/charge/energy,
field/scale calls, zero-scale calls and engine-move counts. Old/new position,
direction and time differences are computed in C++ before serialization.

Summary groups retain profile, ray kind, field, species and accuracy. Pairwise
10-to-20 and 20-to-40 refinements are reported separately for each implementation,
with a hit-changed flag: endpoint changes across different stopping events
cannot be interpreted purely as integration error. Python parses decimals
before subtraction so extended-precision refinement differences do not vanish
through intermediate conversion to double. Serialized summary diagnostics use
ordinary floats; raw states retain their original precision.

Self-contained guards use the existing tests/test_fork_transport.py CI entry
point. They check the unchanged reused grid, 729-case coverage, every
refinement, path independence, duplicate/nonfinite rejection and retention
of disagreement and sub-double decimal differences. Normal builds/CI do not
need the external snapshot or run this optional coupled probe.

See the [confirmed results](../../docs/validation/fork_field_observer_comparison.md). A later source-RNG cost pilot
and statistical field/observer ensemble need their own frozen reservations,
counts and confidence/precision decisions. Do not rerun the spent final seeds
from the photon baseline or infer interacting-cascade equivalence here.

The subsequent [1 PeV electron/photon cost pilot](FIELD_SOURCE_PILOT.md)
reuses this adapter with source sampling. The modern propagation/output
helpers now live in coupled_modern.cpp/.h; fixed-state case inputs and the
frozen results above remain unchanged.
