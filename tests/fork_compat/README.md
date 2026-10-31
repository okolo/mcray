# Historical fork kernel comparison

For the subsequent interaction-free charged-particle trajectory diagnostic,
see [deterministic transport comparison](TRANSPORT.md). The field/source
comparison below remains independently runnable and unchanged in scope.
For strict sphere membership, step limits and interaction-free stopping,
see the [observer diagnostic](OBSERVER.md).
The separately approved [source/observer cost pilot](SOURCE_OBSERVER_PILOT.md)
adds the common pencil/ordinary-cone source samplers without an acceptance verdict.

This optional diagnostic compiles the actual historical `SavedMF` implementation
and uniform-cone initialization block (`fJetType=0`), rather than a transcription
of their formulas. It also compiles the fork's complete Randomizer and Jet
sources. The eight input SHA-256 hashes in `reference.json` pin the read-only
reference snapshot. Changed hashes abort before compilation.

Dependencies: Python >=3.8 standard library; the compiler, GSL and other libraries
from a configured modern **Unix Makefiles** Release build with tests enabled.
Build `grid_magnetic_field_tests` first. The adapter reuses its `link.txt` recipe
and static mcray core. It compiles the current GridMagneticField.cpp and
SpatialGeometry.cpp explicitly, overriding those archive members. Use a build
from the same source revision and configuration: this is not an ABI migration
test. It does not require a private enclosing workspace or a sibling checkout.

From the modern repository root, with paths chosen for the local environment:

```sh
mkdir -p build/fork-comparison/tmp
TMPDIR="$PWD/build/fork-comparison/tmp" PYTHONDONTWRITEBYTECODE=1 \
  python3 tests/fork_compat/run.py \
  --fork-root /read-only/path/to/CRbeam_v0 \
  --build-dir build/Release --output build/fork-comparison/run
```

For a compiler/runtime in an existing offline container, `--command-prefix`
accepts its command as an argv prefix (parsed with `shlex`, never a shell).
For example, prefix with `docker run --rm --network none --read-only` and
explicit project/temp mounts, working directory and image, using an entrypoint
that executes its arguments. Python can then run on the host; it need not be
installed in the compiler image. The prefix executes both compiler and probe,
so their absolute project paths must be accessible unchanged inside it.
The manifest distinguishes the orchestrator environment from this execution
prefix and retains the compiler version.

Keep output and temporary directories inside the project when required by local
workspace instructions. Never place output in the reference directory.

The adapter extracts exact pinned source fragments into the output directory.
Only the namespace, unused interface shells and logging are adapted:

- An empty Utils.h satisfies the isolated old Randomizer include.
- The old field's base and time shells have no behavior used by its methods.
- A shared Mpc constant is allowed only after verifying that old/new Units.cpp
  are identical after removing comments, whitespace and the namespace rename.
- The original cone block includes all branches but executes only type 0. The
  original mutable phi and source weights are retained across histories.
- Reference progress and zero-weight warnings are suppressed by stream buffers;
  the arithmetic and RNG source are unchanged.

`manifest.json` records reference, adapter, current kernel and linked-core
hashes; `compile-command.json` and `compile.log` retain the build invocation.
`inputs.json` hashes generated grid/query files. `field-queries.json` retains
individual field errors; `report.json` retains aggregates and cone statistics.
Missing-table warnings from the modern library startup are harmless here: these
kernels do not read interaction tables.

## Fixed diagnostic design

Two complete regular grids, each containing the origin, have 3x4x5 and 4x3x5
nodes, unequal axes and non-affine field values of order one Gauss. For each of
four boundary modes, compare the same 1000 interior points, 1000 points spanning
inside/outside, and 24 face/seam/outside controls: 16192 queries total. Python's
query seed is 918273. Absolute error is the largest component difference in
Gauss; scaled error divides that difference by `max(1, |B_fork|_inf,
|B_modern|_inf)`. Bitwise counts are diagnostics, not acceptance criteria.

Cone tests use 100000 independent primary draws per implementation and case,
with seeds 17311 (fork ranlux) and 91827 (modern ranlxd2). The seeds are
reused across scenarios; samples are never pooled across those scenarios.
Half-angles are 20,
90 and 180 degrees, plus 1e-6 and 1e-10 radians, each with axis angles
(theta,phi)=(0,0) and (0.4,0.7) radians. Compare weight, weighted first/second
moments of `(1-axis.dot(n))/(1-cos(angle))`, three weighted direction means,
and 20 weighted polar CDF values. The polar offset uses squared chord length
to avoid subtractive cancellation. Means/variances use Welford accumulation;
no ensemble of daughter rows is treated as independent histories.

The 260 differences have simultaneous approximate 95% normal intervals with
Bonferroni allocation over the entire family, including the narrow-cone bug
controls. The interval uses `sqrt(var_fork/N + var_modern/N)`. Zero variance
is reported explicitly; no artificial error floor is added. These are
**diagnostic intervals**, not a user-approved physical equivalence margin or
proof of exact equality. The script reports discrepancies without assigning a
whole-program acceptance verdict. It does not retune counts, bins or intervals
after observing outcomes.

See [results and limitations](../../docs/validation/fork_kernel_comparison.md).

The approved [source/observer ensemble](SOURCE_OBSERVER_ENSEMBLE.md) freezes
the 39-interval family and 48-million-history count. Its compressed streaming
runner defaults to a fixture replay; final collection requires explicit mode.
