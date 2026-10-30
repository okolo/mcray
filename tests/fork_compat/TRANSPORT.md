# Deterministic historical transport comparison

This optional diagnostic extends the [field/source adapter](README.md) with
the actual historical `Deflection3D::GetRotationRate`, `Rate` and `Propagate`
bodies. It reports errors and numerical convergence, without a physical or
statistical equivalence verdict. Ordinary builds and CI do not need the fork.

## Inputs and extraction

Use a configured compatible Unix Makefiles build with the grid test link
recipe and static mcray library. Python >=3.8 standard library orchestrates
the compiler and probe. Compiler, GSL, Boost and Fortran runtime dependencies
are those of the configured modern build. Keep that build current with the
checked-out production source; its static-library hash is recorded.

```sh
mkdir -p build/fork-transport/tmp
TMPDIR="$PWD/build/fork-transport/tmp" PYTHONDONTWRITEBYTECODE=1 \
  python3 tests/fork_compat/transport.py \
  --fork-root /read-only/path/to/CRbeam_v0 \
  --build-dir build/Release --output build/fork-transport/run
```

The existing `--command-prefix` mechanism supports an offline Linux container;
see [the adapter instructions](README.md). Both compilation and execution
use that prefix. Mount only the current repository and its output/temp paths;
the compiler uses copied/extracted reference sources, not the external fork.
An ASan/UBSan build is supported by selecting its build directory and exporting
the same sanitizer runtime options for the probe. This instruments the extracted
historical functions too. Each probe invocation has a 120 s execution timeout.

`transport_reference.json` pins twelve reference files. Extraction preserves
the original particle enum, light-particle mass/charge tables and functions,
double-precision `beta`, tracer function, SavedMF methods and transport bodies.
The only enabled species are named Electron, Positron and Proton; nuclear
branches remain in the extracted functions but reject if reached. The tracer
stream is null. Reference RNG sources are compiled but Propagate consumes no
random numbers. Driver shells supply only fields/lifetimes and the particle
members accessed by those original functions.

The complete historical and modern CosmoTime classes are compared after
comment/whitespace/namespace normalization before sharing that type. Unit
implementation identity is checked by the existing adapter. Both paths receive
the same raw initial time and energy. The modern cosmology is initialized at
H0=70 km/s/Mpc, OmegaLambda=0.7, maximum z=10.01 and default accuracy 100000;
the starting time corresponds to z=0.01. This is not an independent comparison
of cosmological integration. No interaction, adiabatic loss, observer, source
ensemble or engine scheduling is involved.

## Frozen diagnostic design

`inputs.json` and the case/grid hashes are written before compilation and
trajectory evaluation. There are 756 cases:

- electron, positron and proton;
- accuracy 10, 20 and 40, crossed with 1, 4 and 16 equal outer calls;
- elapsed c*time 0.25 Mpc, common normalized direction (0.6,0,0.8);
- zero field, axial 1 nG field and an oblique (0.2,0.3,0.5) nG field;
- a mildly relativistic constant-field control with common E=10 times the
  modern species mass, and amplitude selected for angular rate 0.5/Mpc;
- a complete centered 3x4x5 grid spanning x=[-0.2,0.2], y=[-0.3,0.3],
  z=[-0.4,0.4] Mpc, with six inside/face/corner/near-outside/far-outside starts
  and all four extension modes.

Except the mildly relativistic control, energy is 1 EeV. The grid comes from
a nonuniform divergence-free analytic field (each component is independent of
its own coordinate); the trilinear interpolant retains that property. The
constant-field scale is explicitly positive at 0.1 Mpc. Field energy/charge
identity is checked; unexpected energy changes or incomplete/nonfinite results
abort the report. The fork's zero-scale shortcut is retained and counted.

In constant fields, compare against the analytic helix with the correct sign
for dn/dt=(q e/E) n cross B. Position uses each implementation's own beta,
so beta rounding is not confused with integrator error. The modern and old
beta values and their difference are reported separately. Before processing
cases, the analytic formula is checked against free motion and quarter-turns
for both charge signs, with a tolerance derived from double roundoff.

## Output and checks

- `trajectories.jsonl`: every case, raw coordinates/directions, actual elapsed
  time, beta, energy, mass/charge, scales, field/scale call counts, norm/time
  errors, analytic errors where available, and old/new differences. Original
  21-digit C++ decimal output is retained without Python float rewriting.
- `report.json`: grouped error maxima for every accuracy/outer-call combination,
  and total counts. It never selects a best run or assigns acceptance.
- `manifest.json`, `fragment-sha256.json`: exact reference/adapter/modern-source
  and linked-library hashes, Git revision/dirty status, environment/prefix and
  verbatim extracted-fragment hashes.
- `inputs.json`, `transport.cases`, `transport.grid`: frozen design and data.
- `compile-command.json`, `compiler-version.txt`, `compile.log`, `execution.log`:
  actual commands, toolchain and diagnostics. Failed reruns remove stale reports.

Self-contained Python guards do not need the fork or a compiled library:

```sh
TMPDIR="$PWD/build/fork-transport/tmp" PYTHONDONTWRITEBYTECODE=1 \
  python3 tests/test_fork_transport.py -v
```

They check complete/fixed case coverage, path-independent inputs, rejection of
changed reference hashes, duplicate/nonfinite trajectories and retention of
all refinements, plus preservation of the original compiler-helper defaults.
Confirmed execution results belong in the public validation documentation.
