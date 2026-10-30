# Grid magnetic fields and field composition

The C++ library provides `Interactions::SavedMF` and `Interactions::ComboMF`
in `GridMagneticField.h`. They extend the existing `MagneticField` interface
without changing the wave fields, particle identity, random generator or
simulation defaults. The CRbeam CLI and YAML version 2 can select them as
specified below; version 1 remains compatible.

## Input and units

Construct `SavedMF(filename, boundary_mode)` with a text file containing six
columns per node:

```text
# x/Mpc y/Mpc z/Mpc Bx/G By/G Bz/G
0 0 0 0 1e-9 0
0 0 1 0 1e-9 0
0 1 0 0 1e-9 0
0 1 1 0 1e-9 0
1 0 0 0 1e-9 0
1 0 1 0 1e-9 0
1 1 0 0 1e-9 0
1 1 1 0 1e-9 0
```

Coordinates are comoving Mpc; components are Gauss. File paths are used as
supplied, relative to the caller's working directory when not absolute.
Rows may occur in any order. Blank lines, full-line comments and trailing
`#` comments are accepted. Other trailing columns are errors. Numeric
parsing uses the classic locale with a decimal point.

Every axis must have at least two distinct, regularly spaced coordinates;
axes may have different spacings, and the box need not contain the origin.
The file must contain exactly one node at every Cartesian combination.
Missing, duplicate, malformed or nonfinite entries are fatal errors. Axis
spacing is checked against its mean with tolerance
`1e-10 * step + 64 * epsilon(coord_type) * max(abs(axis endpoints))` in Mpc.
This allows decimal rounding; nodes are not snapped or synthesized. Axis
extents, internal coordinates and the finite positive minimum scale must
remain representable. Node-count products are checked before allocation.

Queries use the existing interface: `x` has three coordinates in mcray
internal length units, and the supplied output vector must already contain
three components. The field returns Gauss. Invalid queries are fatal and
cannot be bypassed by retrying a cascade. Interior values use trilinear
interpolation over the actual enclosing nodes, with local cell fractions
rather than a polynomial in absolute coordinates.

This is a **static field**: `CosmoTime` does not change its values or scales.
No `(1+z)` rescaling, unit conversion of field components or automatic
physical/comoving conversion is applied. Interpolation does not generally
preserve a discrete divergence constraint; no solenoidal projection is
performed and no physical field model is inferred from the file.

## Boundary modes

To reproduce the historical `SavedMF` model, trilinear interpolation applies
strictly inside the box. Touching any face invokes the selected extension
for **all three coordinates**:

| Mode | On or outside the original box |
| --- | --- |
| `SavedMF::FILLING_ZERO` (default) | Return the zero vector. |
| `SavedMF::FILLING_NEAREST` | Round every coordinate to the nearest Cartesian node, clamping outside indices. Halfway ties round toward the higher index. |
| `SavedMF::FILLING_REPEAT` | Wrap every coordinate modulo its box width. |
| `SavedMF::FILLING_MIRROR` | Fold every coordinate with a period twice its box width. |

Repeat and mirror retain the fork's seam rule: after mapping, a coordinate
at/below the lower face moves inward by `width/node_count`; a coordinate
at/above the upper face moves inward by the same amount. Mirror folds
positions without reversing component signs. These displacements can be a
substantial fraction of a cell. They are compatibility behavior, not an
exact periodic interpolation across matching endpoints.

The inset applies to the **mapped** endpoint, not the original query's face.
For example, on a zero-based axis `[0,L]`, repeat maps `x=L` to zero before
applying the lower inset, giving `L/N`. Mirror maps the same point to `L`
before applying the upper inset, giving `L-L/N`. This distinction is part
of the historical seam behavior; floating-point remainder rounding near a
seam can also select its side.

Box limits and periods use the fork's double-precision arithmetic in internal
units, including the order of subtraction before the remainder. This matters
at a seam: a rounding difference in the remainder can select a finite inset.
Interpolation retains the numerically stable weighted-corner formulation.
Nearest extension rounds the two other axes too, even when they are inside.
All modes may therefore be discontinuous. Choose the boundary mode with this
field model in mind; the propagator does not locate exact face crossings.

`MinVariabilityScale` returns `min(width_x/Nx, width_y/Ny, width_z/Nz)` in
internal length units at every position and time, including zero-field
regions. This restores the historical extent/node-count convention, with a
correct minimum when axes tie. The fork incorrectly selected the z scale
for some ties. A zero sentinel outside the box is deliberately not restored:
the modern deflector requires a positive scale for subdivision, including
trajectories entering the field. This scale alone does not establish
trajectory convergence near a discontinuity.

The robust loader still uses actual extrema and rejects malformed/incomplete
grids. It does not reproduce the fork's zero-initialized extrema or undefined
uninitialized-node behavior. See the [direct comparison and exceptions](
../validation/fork_kernel_comparison.md).

## Composition and ownership

```cpp
#include "GridMagneticField.h"

Utils::SmartPtr<Interactions::MagneticField> grid =
    new Interactions::SavedMF("field.txt", Interactions::SavedMF::FILLING_ZERO);
Utils::SmartPtr<Interactions::MagneticField> waves =
    new Interactions::TurbulentMF(randomizer, 0.05, 5.0, 1e-15, 300);
Utils::SmartPtr<Interactions::MagneticField> sum =
    new Interactions::ComboMF({grid, waves});
Interactions::Deflection3D deflection(sum, 10);
```

`ComboMF` retains each component through the existing intrusive `SmartPtr`.
Do not give it pointers to stack objects. It rejects empty and null component
lists, sums field vectors and returns the smallest component scale at the
requested time. Each scale must be finite and strictly positive; an invalid
scale is a fatal error rather than an instruction to skip that component.
Malformed/nonfinite component vectors and nonfinite final sums are also
fatal. Composition does not normalize the result or alter component models.

`SavedMF` and `ComboMF` are assembled once and have no mutation API. Their
queries use local scratch storage and support concurrent reads. Other
component classes must also support concurrent const queries. The existing
intrusive reference counts are not atomic: assemble fields, create clones
and release references in serial lifetime-management code; concurrent reads
do not copy component smart pointers. PropagationEngine serializes both worker
construction and destruction under the same OpenMP critical section, including
references retained by Deflection3D clones. Other callers must provide the
same lifetime synchronization when cloning or releasing shared fields. No owning pointer is added to
`Particle` and no RNG is used during grid loading or field queries.

## Verification

`grid_magnetic_field_tests` uses temporary synthetic files and the production
library. Run the bounded checks from the repository root:

```sh
cmake --build build/Release --target grid_magnetic_field_tests --parallel 2
ctest --test-dir build/Release -R '^unit.magnetic\.' --output-on-failure
```

The checks cover affine interpolation, box faces and translated bounds, all
four extensions, malformed grids/queries, component lifetime, finite positive
scales, concurrent reads and shared deflector clones. A uniform-field control
checks an electron and an explicit Fe56 identity without changing their
energy or weight. These are implementation controls, not physical acceptance
of an external field realization or a long propagation calculation.

## CRbeam CLI and YAML

CRbeam now accepts a grid through the following long options. Relative CLI
file paths resolve against the working directory:

```sh
CRbeam --particle proton --z 0.01 --emax 1e19 \
  --magnetic-field grid --magnetic-field-file ./field.txt \
  --magnetic-field-boundary mirror --accEGMF 10 --seed 42 --output results
```

`--magnetic-field` accepts `none`, `monochromatic`, `kolmogorov` and `grid`.
The two other new options require `grid`. Its file is mandatory and boundary
choices are `zero` (default), `nearest`, `periodic` and `mirror`. No short
aliases are introduced. Compositions are configured through YAML only.

Existing CLI invocations keep their previous settings and RNG path. For
explicit wave models, existing `--EGMF`, `--lmaxEGMF`, `--lminEGMF`,
`--nmodesEGMF`, `--randomEGMF` and `--accEGMF` supply the applicable settings.
`monochromatic` rejects an explicit `--turbulentEGMF`, minimum scale or mode
count. `grid` and `none` reject explicitly supplied strength/spectral/random
flags; a field file is never silently rescaled. Duplicate new options and
invalid numeric arguments are errors. `--accEGMF` remains the accuracy option.

YAML version 1 retains its existing magnetic-field schema, zero-strength
switch and behavior. Version 2 uses a model-specific mapping:

```yaml
schema_version: 2
source:
  redshift: 0.01
  particle: {type: proton}
  spectrum: {type: monoenergetic, energy_eV: 1.0e19}
magnetic_field:
  model: sum
  realization: shared
  deflection_accuracy: 10
  components:
    - model: grid
      file: ./field.txt
      boundary: zero
    - model: kolmogorov
      strength_G: 1.0e-15
      minimum_scale_Mpc: 0.05
      maximum_scale_Mpc: 5.0
      modes: 300
output:
  directory: ./results
```

Relative grid paths resolve against the input YAML directory. Other sections
retain the version 1 format. The following model keys are supported:

| Model | Model-specific keys |
| --- | --- |
| `none` | None; explicitly disables deflection. |
| `grid` | Required `file`; optional `boundary` defaults to `zero`. |
| `monochromatic` | `strength_G`, `maximum_scale_Mpc` (wavelength). |
| `kolmogorov` | `strength_G`, `minimum_scale_Mpc`, `maximum_scale_Mpc`, `modes`. |
| `sum` | Required nonempty `components` sequence of grids and wave fields. |

The top-level mapping also accepts `deflection_accuracy` and `realization`.
Accuracy defaults to 10. Realization choices are `shared` (default) and
`per_primary`. Pure `grid`/`none` require `shared`. For a stochastic sum,
`per_primary` regenerates only the wave components for each primary history;
its grids stay fixed. The same field realization is shared by that primary's
secondary cascade, as in the existing per-primary wave-field mechanism.

Components do not accept their own realization/accuracy, nested sums or
`none`. Inapplicable/unknown/duplicate keys are errors. Wave defaults remain
those of the legacy command model: strength 1e-15 Gauss, maximum scale 5 Mpc,
minimum scale 0.05 Mpc and 300 modes where applicable. Kolmogorov requires at
least two modes and minimum scale strictly below maximum. Explicit zero
strength is allowed for a wave component and does not disable its neighbors.
If the entire field mapping is omitted in version 2, the legacy default is
preserved; there is no hidden change to the default model or amplitude.

Every grid is loaded and validated before creating simulation output. The
application freezes its file bytes, hashes and parses those same bytes once,
and retains the immutable grid. Repeated components using the same file and
boundary reuse the loaded field; a change detected between separate boundary
loads is an error. Neither propagation nor per-primary construction reopens
the file. The loader uses temporary memory for the file bytes in addition to
the grid's parsing/storage buffers.

The effective `config.yaml` records version 2, resolved absolute file paths,
model defaults and applicable settings. `run.yaml` records the SHA-256 of each
loaded grid under `magnetic_field.files_sha256`, without duplicating settings.
The hashes identify the bytes parsed for that execution; snapshots do not
embed grids or require a prior run's hash for a later intentional rerun.
