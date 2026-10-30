# mc_rotate.py

[Postprocessing utility index](README.md)

mc_rotate.py rotates the original +z-pencil output and extrapolates each
recorded trajectory along its final momentum to the sphere R=cT. The source
must be at the origin, with the same birth redshift >0 and end redshift zero
for all histories; T must be the actual common flight time. This routine is
not for arbitrary source/observer output or fallback normalization scales.
The extrapolation assumes speed c and no further interactions or deflections.
The first two columns, the named rotated columns and trailing input columns
retain the original order. Angles are in degrees; delays are Delta_t/T.

```python
import sys
sys.path.insert(0, 'tools/postprocess')
from mc_rotate import mc_rotate
output = mc_rotate('output/run/z0/photon', 5, 1, psf_deg=0.1,
                   block_size=32768)
```

## Algorithm

The calculation uses c=1. Let T be the common elapsed flight time at the
recorded endpoint. A straight photon emitted at the same instant reaches a
sphere of radius T in time T. Each simulated trajectory is extended from
its recorded endpoint to that sphere along its final momentum direction.
The additional flight time is therefore the delay relative to that straight
photon. This is a geometric continuation; the script does not simulate new
interactions, magnetic deflections or cosmological evolution on the extension.

**1. Decode position and momentum.** The first eight numeric input columns are:

| Column | Quantity |
| --- | --- |
| 1 | Energy in eV |
| 2 | Statistical weight |
| 3, 4 | x/T, y/T |
| 5 | d = 1-z/T |
| 6, 7 | n_x, n_y |
| 8 | 1-n_z |

Remaining columns are carried into the output in their original order.
The direction n is reconstructed and normalized. With x and y now denoting
normalized coordinates, the recorded position is r=(x,y,1-d). All later
positions are also normalized by T.

**2. Find the extra flight time.** A forward continuation has position
r + tau*n, where tau=Delta_t/T. The sphere condition gives

\[
|r+\tau n|^2=1,
\qquad \tau^2+2q\tau-\delta=0,
\]

where

\[
q=r\cdot n,
\qquad \delta=1-|r|^2=d(2-d)-x^2-y^2.
\]

The deficit delta is evaluated directly from the saved small complement d.
Computing 1-|r|^2 after reconstructing z/T would lose small delays when the
radius rounds to one. Short sums use Neumaier compensation independently
for every particle and retain np.longdouble throughout.

For h=sqrt(q^2+delta), evaluated with hypot(q,sqrt(delta)), the root is

\[
\tau=
\begin{cases}
\displaystyle\frac{\delta}{h+q}, & q>0,\\[6pt]
h-q, & q\le0.
\end{cases}
\]

The first branch avoids subtracting nearly equal numbers for outward motion.
The second handles inward and tangent motion without that cancellation.
The endpoint u=r+tau*n is normalized to remove residual rounding error.
For the axial case r=(0,0,1-d), n=(0,0,1), the result is exactly tau=d in
real arithmetic, even for delays much smaller than machine epsilon.
For an inward direction starting exactly on the sphere, the positive root
selects the later exit on the opposite side; a zero tangent root remains zero.

**3. Rotate the endpoint to the observer axis.** A rotation R maps u to +z.
The complete history is interpreted in this rotated frame: the original
initial direction becomes p_initial=R*(0,0,1)=(-u_x,-u_y,u_z), and the final
momentum becomes R*n. This is a geometric rotation of the recorded history,
not a new propagation calculation in an independently fixed external field.
At the north pole R is the identity; at the south pole the convention is a
pi rotation about x. Rodrigues' formula is evaluated by components, without
allocating a 3x3 matrix for each row. Near the north pole, 1-cos(theta) is
computed as sin(theta)^2/(1+cos(theta)) to retain tiny transverse rotations.

**4. Apply the jet and PSF selections.** For the supplied jet inclination psi,
the jet axis in the rotated frame is j=(0,sin(psi),cos(psi)). Retain particles
whose initial direction lies strictly inside the specified jet half-angle.
The legacy convention treats nonpositive half-angles and angles >=180 degrees
as full-sky selection. The PSF cut retains particles with arrival angle
Theta <= psf_deg, where

\[
\Theta=\operatorname{atan2}(|u\times n|,u\cdot n).
\]

Angles between initial direction and jet axis use the same atan2 construction.
It retains small angles that arccos of a nearly unit dot product would lose.
Azimuth Phi is calculated from the x/y components of R*n with atan2 and mapped
to [0,360) degrees; Phi is set to zero when Theta is zero. Jet and PSF cuts
are applied before computing the rotated final components and azimuth.

**5. Write selected rows.** The output columns are:

```text
E/eV weight Theta Phi dT_raw/T dT_calc/T [original columns 9 onward]
```

Theta and Phi are in degrees. dT_raw/T is the stored d, with tolerated negative
roundoff corrected to zero. It is the original longitudinal deficit, not a
general sphere-arrival delay. dT_calc/T is the calculated tau. Multiplying tau
by T gives the estimated additional flight time; if the header provides
T/Mpc in c=1 units, the corresponding time is tau*(T/Mpc)*Mpc/c.

Selected rows retain input order and weights. Comments and column headers
are copied or updated in their original position. NumPy's scalar formatter
preserves the actual longdouble precision instead of limiting the output to
six significant digits. The script does not use a summation or formatting
helper that narrows extended scalars to float64.

## Numerical limits and validation

Negative d or radial deficit is clamped only within the specified roundoff
bound. Larger violations, non-finite rows and zero momentum directions raise
ValueError with the input line number. Positive deficits are never clamped to
zero. By default the bound uses eight longdouble epsilons, scaled by the
magnitudes of the deficit terms. For older low-precision text files, pass an
explicit input_rtol, for example `input_rtol='5e-6'` when justified by the
file's precision. The parameter represents uncertainty in normalized input
coordinates and widens the bounds; it cannot recover information lost when
the input was written. Cancellation between nearly equal deficit terms still
limits the relative accuracy of their difference.

## Tests

Run the deterministic numerical and file-format checks (NumPy only):

```sh
PYTHONDONTWRITEBYTECODE=1 python3 tests/test_mc_rotate.py
```

## Block processing

The file reader buffers at most block_size data rows and uses NumPy to parse
and process a complete block. Workspace scales with block size and column
count, not total file length. The default is 32768 rows; smaller blocks reduce
memory at the cost of more calls and parsing overhead. Formatting streams
accepted rows without accumulating the complete output in memory.

The numerical operation is also exposed separately for arrays:

```python
from mc_rotate import rotate_block
rotated_rows = rotate_block(particles, 5, 1, psf_deg=0.1)
```

particles must be a two-dimensional numeric array in the original column
order. The function does not modify it. The returned longdouble array contains
only selected rows in input order; an empty selection retains its column
shape. Calculations are vectorized across rows, including compensated sums.
Rotation components are computed directly rather than allocating a 3x3 matrix
per particle. Jet/PSF cuts precede the azimuth calculation. Numerical error
objects carry a zero-based block row index; file processing maps this index
to the actual input line, including blank/comment lines. Column counts must
remain consistent across blocks.

Tests create temporary files under build/mc-rotate-tests, inside the project.
See [verification results](../../docs/validation/mc_rotate.md) for the
high-precision reference and platform coverage.
