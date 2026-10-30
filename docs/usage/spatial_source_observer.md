# Source geometry, sphere observers and spatial output

YAML schema version 2 and the long CLI options below configure source geometry.
Existing CLI invocations and version 1 keep their defaults and output format.
A version 2 configuration containing only default spatial values also keeps
that format, including when replayed from config.yaml.

| CLI option | YAML field | Default / meaning |
| --- | --- | --- |
| `--source-position-Mpc x,y,z` | `source.position_Mpc` | `[0, 0, 0]`, absolute engine coordinates in Mpc |
| `--source-angular-model pencil\|cone` | `source.angular_distribution.model` | `pencil` |
| `--source-axis x,y,z` | `source.angular_distribution.axis` | `[0, 0, 1]`, finite nonzero direction, normalized |
| `--source-half-angle-deg angle` | `source.angular_distribution.half_angle_deg` | `0`; cone supports 0 through 180 degrees |
| `--end-redshift z` | `propagation.end_redshift` | `0`, nonnegative and below source redshift |
| `--observer none\|sphere` | `observer.model` | `none` |
| `--observer-center-Mpc x,y,z` | `observer.center_Mpc` | `[0, 0, 0]`, sphere only |
| `--observer-radius-Mpc radius` | `observer.radius_Mpc` | Required positive finite radius for sphere |
| `--observer-stop-when inside\|outside` | `observer.stop_when` | `inside`, sphere only |
| `--output-primary-direction` | `output.primary_direction` | `false`, record primary direction |

YAML vectors are sequences of three numbers. CLI vectors are comma-separated;
quote the complete value if it contains spaces. Geometry and time settings are
validated before output creation. Source coordinates, fields and observer
coordinates use the same existing engine convention; no new cosmological
coordinate transformation is applied.

A cone samples uniformly in solid angle about its axis. Pencil and zero-width
cone sources consume no extra random values; geometry does not change primary
energy, weight or identity. Nonuniform jet distributions are not implemented.

## Sphere semantics

The first version preserves the fork's strict endpoint membership rule:
inside detects `distance < radius`; outside detects `distance > radius`.
Equality detects neither. A primary born in the selected region is recorded
immediately, subject to the existing output energy filter. Detection stops the
trajectory. Particles reaching end redshift without detection are also written,
with `in_target = 0`, if they pass the output filter.

This is sampled membership, not an exact crossing solver. Step limiting uses
the fork's radius fractions (0.001 inside, 1e-7 outside), capped by remaining
simulation time. A step below representable time resolution is an error.
One shared immutable observer controls termination, step limiting and the
recorded `in_target` flag. The API permits a later observer implementation;
crossing events and multiple-hit logging are not implemented.

See [the runnable configuration](../../examples/spatial_observer.yaml):

```sh
mkdir -p output
(cd bin && ../build/Release/CRbeam ../examples/spatial_observer.yaml)
```

## Spatial table contract

Nondefault spatial settings select `spatial_format = 1`, identified in the
first header line and run.yaml. This layout preserves the fork's column order:

1. `E/eV`, `weight`, `x/T`, `y/T`, `1-(z/T)`, `Px/P`, `Py/P`, `1-(Pz/P)`, `z_production`.
2. `in_target` when a sphere is enabled.
3. `source_x`, `source_y`, `source_z` when primary-direction output is enabled.
   These are dimensionless primary momentum-direction components.
4. `z_src`, `E_src/eV`; `Id` when trajectory logging is enabled.
5. `PathMpc`, followed by the new last column `FlightMpc`.

`PathMpc` is elapsed time since this particle's production. `FlightMpc` is
elapsed time since the primary's birth. Both are existing simulation time
intervals expressed as lengths at c=1, divided by units.Mpc; neither is
inferred from displacement. Source redshift/energy remain available for each
row. Optional columns keep their relative order. Ordinary runs retain the
modern RawOutput3D layout, including its interaction columns.

Each output filename, including each individual isotope filename, has an
independent immutable scale. The first written particle selects T:

- If its primary flight time is positive, T is exactly that interval, as in
  the fork.
- If it is zero, T is `max(1 Mpc, abs(x0), abs(y0), abs(z0))`, using the first
  recorded position. Later positive flight times do not replace this scale.
- Negative/non-finite flight times or unrepresentable normalized coordinates
  are errors.

The first line contains `T/Mpc = ...` and either
`normalization = first_particle_flight` or `normalization = fallback`.
T is printed as an ordinary decimal for old regex-based readers. Values use
max_digits10 precision for the coordinate type. The same line records primary
birth position as `source_x_Mpc`, `source_y_Mpc`, `source_z_Mpc` for standalone
coordinate reconstruction. The tab-separated column header remains the last
comment line before numeric data. No additional columns precede FlightMpc.

For S from `T/Mpc`, reconstruct coordinates in Mpc as
`x = column(x/T)*S`, `y = column(y/T)*S`,
`z = (1-column(1-(z/T)))*S`. The stored complement can still lose precision
when z/T is near zero; increased text precision does not eliminate cancellation.
