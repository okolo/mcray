# Truncated von Mises-Fisher source

The optional `vmf` source samples an axisymmetric von Mises-Fisher angular
profile conditioned on a cone. It preserves the incoming primary weight,
energy and identity. The angular density integrates to one inside the chosen
cone, so changing kappa or the opening does not rescale the source's total
angular normalization. The existing energy/redshift sampling weights retain
their meanings. This describes one jet; no second opposite jet is added.

## Input

```sh
CRbeam --particle electron --z 0.01 --emax 1e13 \
  --source-angular-model vmf --source-axis 0,0,1 \
  --source-half-angle-deg 10 --source-vmf-kappa 100 \
  --seed 123 --output vmf-run
```

The equivalent angular settings in YAML schema version 2 are:

```yaml
source:
  angular_distribution:
    model: vmf
    axis: [0, 0, 1]
    half_angle_deg: 10
    kappa: 100
```

[examples/vmf_source.yaml](../../examples/vmf_source.yaml) is a complete
configuration. It uses the existing CMB processes with EBL and magnetic field
disabled and includes primary-direction output.

| Setting | Contract |
| --- | --- |
| `model` / `--source-angular-model` | Adds `vmf` alongside unchanged `pencil` and `cone`. |
| `kappa` / `--source-vmf-kappa` | Required explicitly for vmf, finite and >=0, dimensionless. Forbidden with other CLI/YAML models, even if zero. |
| `half_angle_deg` / `--source-half-angle-deg` | Required explicitly for vmf, finite in [0,180]. Hard outer cutoff, not a standard deviation. |
| `axis` / `--source-axis` | Existing normalized nonzero finite vector; default [0,0,1]. |

An explicitly supplied zero opening gives the axis direction for any finite
kappa. Kappa=0 gives the uniform-in-solid-angle cone. Positive kappa
concentrates emission toward the axis, and a 180-degree opening covers the
full sphere. Negative kappa is unsupported. Missing, duplicated, malformed,
nonfinite or inapplicable settings fail before output creation. Version 1
retains its original schema and does not accept this angular model/parameter.

The effective model, normalized axis, opening and kappa are saved in
config.yaml. run.yaml retains the existing configuration hash. The model
uses the existing spatial-output path, including zero opening, with no new
columns. `output.primary_direction` / `--output-primary-direction` retains
its current optional meaning; see [the spatial contract](spatial_source_observer.md).
Readers can inspect the complete angular settings through config.yaml.
The legacy +z-pencil assumptions of mc_rotate do not extend to these sources.

## Distribution and sampler

For unit axis a and direction n, write u=1-a dot n and
U=1-cos(theta_max)=2*sin(theta_max/2)^2. For U>0 and kappa>0:

```text
f_Omega(n) = kappa*exp(-kappa*u) / (2*pi*(1-exp(-kappa*U))),  0<=u<=U.
H(u)       = (1-exp(-kappa*u)) / (1-exp(-kappa*U)).
```

The azimuth is uniform. For kappa=0 the density is 1/(2*pi*U); for U=0 the
source is the deterministic axis. The density is with respect to solid angle,
dOmega=dphi*du. Its normalization is conditional on the chosen cone; it is
not an unnormalized full-sphere intensity restricted to that cone.

For a uniform draw r, the implementation samples

```text
u = -log1p(-r*(-expm1(-kappa*U))) / kappa.
```

It uses the existing orthonormal frame and direction normalization. Sine is
computed as sqrt(u*(2-u)), preserving transverse components when 1-u rounds
to one. For a positive opening there are exactly two RNG draws: offset and
azimuth. Kappa=0 reuses the existing cone arithmetic, giving identical
directions and next RNG value for the same modern inputs. Zero opening uses
no draws and matches the pencil direction. Legacy pencil/cone arithmetic,
RNG and defaults are preserved.

The constructor caches the truncated mass. When kappa*U<=epsilon(coord_type),
the relative uniform-limit correction is below the coordinate precision and
the sampler uses r*U. Large kappa*U is handled without overflowing the product
or computing a positive exponential. Near zero, log1p is evaluated with the
rapidly convergent identity log(1+x)=2*atanh(x/(2+x)), stopping when the sum
no longer changes. This avoids observed small-argument x87-emulation loss
of precision; larger arguments use the standard function. Possible edge
rounding is bounded by U.
No rejection loop, PDF epsilon, probability floor, weight correction or
arbitrary physical kappa cap is introduced.

For the new model, a positive opening whose width underflows to zero is an
input error. Finite coordinate/RNG precision still limits extreme tails and
very narrow structures. For a rotated axis, angular perturbations smaller
than its representable component spacing can disappear. Large finite kappa
therefore does not promise an exactly resolved angular profile at arbitrary
scales. Equal seeds also do not promise identical histories across platforms.

The C++ SourceGeometry API adds AngularModel::Vmf and an optional kappa
constructor argument (zero is the unused default for other C++ models).
Substantial implementation stays in SpatialGeometry.cpp; the immutable
settings are safe to query with worker-local RNGs and particles.

The fork's vMF branch sampled uniform cone directions and multiplied their
weights by a full-sphere density. This implementation follows the approved
normalized source contract above; it does not reproduce those absolute
weights, Gaussian-in-1-cos(theta) sampling, or the fork's RNG streams.

For short-check scope and results, see [vMF verification](../validation/vmf_source.md).
