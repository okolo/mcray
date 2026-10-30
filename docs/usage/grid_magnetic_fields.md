# Grid magnetic fields and field composition

The C++ library provides `Interactions::SavedMF` and `Interactions::ComboMF`
in `GridMagneticField.h`. They extend the existing `MagneticField` interface
without changing the wave fields, particle identity, random generator or
simulation defaults. The current CRbeam CLI/YAML interface does not select
these new classes yet.

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

The original box is closed, including its faces and corners. An exact node
at the maximum endpoint retains that endpoint's value. Extensions apply
outside the box independently on each axis:

| Mode | Outside the original box |
| --- | --- |
| `SavedMF::FILLING_ZERO` (default) | Return the zero vector if any coordinate is outside. |
| `SavedMF::FILLING_NEAREST` | Clamp each coordinate to the nearest box face, then interpolate. |
| `SavedMF::FILLING_REPEAT` | Wrap each coordinate modulo the box width. |
| `SavedMF::FILLING_MIRROR` | Fold coordinates with a period twice the box width. |

Nearest extension retains interpolation in the remaining axes; it does not
round to the nearest node. Mirror extension folds positions and does not
reverse magnetic component signs. Periodic endpoints are not automatically
made equal: incompatible opposite faces produce a seam. Zero extension can
also be discontinuous at the boundary. Choose data and boundary modes
appropriate to the intended physical model.

`MinVariabilityScale` returns the smallest node spacing in internal length
units at every position and time, including zero-field regions. It never
returns a zero sentinel. The existing deflector can therefore use its
normal step subdivision, including for trajectories entering the box. This
scale alone does not locate an exact boundary crossing or establish
trajectory convergence near a discontinuity.

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
do not copy component smart pointers. Existing `Deflection3D::Clone` shares
the field through that ownership mechanism. No owning pointer is added to
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
