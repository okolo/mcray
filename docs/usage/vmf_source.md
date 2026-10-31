# Fork-compatible weighted vMF source

The `vmf` model follows CRbeam_v0: sample a direction uniformly in the
configured cone, then multiply the incoming primary weight by `vMF_weight`.
It does not directly sample a cone-conditioned vMF probability distribution
and does not preserve total angular weight when kappa or opening changes.
This restores the original contract at the user's request on 2026-10-31,
superseding the normalized conditional-vMF implementation from PR #12.

## CLI and YAML

The existing names and schema version 2 are retained:

```sh
CRbeam --particle electron --z 0.01 --emax 1e13 \
  --source-angular-model vmf --source-axis 0,0,1 \
  --source-half-angle-deg 10 --source-vmf-kappa 100 \
  --seed 123 --output vmf-run
```

```yaml
source:
  angular_distribution:
    model: vmf
    axis: [0, 0, 1]
    half_angle_deg: 10
    kappa: 100
```

The [complete example](../../examples/vmf_source.yaml) enables primary-direction
output. Kappa must be explicitly supplied and finite, and may be positive,
negative or zero, as allowed by the original weight function. A negative
value favors the opposite pole. Half-angle must be explicitly supplied in
[0,180] degrees. A finite nonzero axis is normalized. Kappa remains forbidden
with pencil/cone models. Schema version 1 does not accept vmf/kappa.
Effective settings are saved in config.yaml and run.yaml identifies code and
execution provenance. This changes results for existing nonzero-kappa vMF
configurations from the superseded normalized implementation; replay must
use the original recorded code revision when that old behavior is intended.
No second normalized-vMF mode or new selector is introduced.

## Original sampling and weight contract

For unit axis a, direction n, u=1-a dot n and U=1-cos(theta_max), the proposal
is uniform in solid angle on the cone: q=1/(2*pi*U). Each active draw uses two
random values, one for u and one for azimuth. The weight update is:

```text
if abs(kappa) < 1e-14:
    multiplier = 1
else:
    multiplier = kappa*exp(-kappa*u) / (2*pi*(1-exp(-2*kappa)))
weight *= multiplier
```

This is the original full-sphere density multiplier. Do not multiply by
2*pi*U, divide by q, replace the denominator by 1-exp(-kappa*U), or otherwise
renormalize the total source. The special small-kappa branch is preserved
exactly, including its jump from multiplier 1 to a limit of 1/(4*pi).
Outside the small-kappa branch, for positive kappa and an active opening
with U>0, the expected angular multiplier is:

```text
(1-exp(-kappa*U)) / (2*pi*U*(1-exp(-2*kappa))).
```

Thus changing source parameters can change absolute flux. Angularly
normalized observables and absolute weighted output must be distinguished.
Large concentration with a broad uniform proposal can have very high weight
variance; it does not acquire the variance properties of direct vMF sampling.
Initialize multiplies the current weight once; each primary starts with the
weight from its energy/redshift sampler, as in the fork source loop.

The fork activates angular customization only for opening >1e-14 radians.
Vmf retains this activation threshold: zero/smaller openings consume no
angular random draws and leave the weight unchanged, regardless of kappa.
The already accepted explicit-axis interface still orients the pencil along
its axis; the fork's ignored zero-opening rotation is not restored. Ordinary
pencil/cone behavior and their RNG consumption are unchanged. For active
opening and small kappa, vMF and cone use exactly the same modern directions
and next RNG value for the same inputs.

## Code structure and numerical exceptions

`src/lib/Jet.h` and `Jet.cpp` restore the original `vMF_weight` function name
and angle-based entry point. A small axis overload adapts the existing
vector-axis CLI. SourceGeometry keeps the recognizable fork sequence:
uniform cone direction, then vMF weight multiplication. There is no inverse
vMF CDF, conditional-mass cache or direct-vMF sampler. Current particle
identity, ownership and RNG interfaces are retained. Equal seeds do not
produce the fork's random sequences. Finite RNG/direction precision still
limits extreme concentrations; pointwise weight checks do not establish
statistical agreement for unresolved rare directions.

The original formula is evaluated with documented numeric corrections:

- `expm1` prevents cancellation in 1-exp(-2*kappa) near the threshold.
- Reflecting the cosine for negative kappa avoids overflowing exp(+2*abs(kappa)).
- Unit-vector dot products may be clipped only for small rounding beyond a
  pole; substantial range violations are errors. Stable cone-width/sine
  arithmetic and the existing source frame remain.
- Log-space fallback retains a representable density when exp(-kappa*u)
  alone underflows. A truly unrepresentable tail remains zero.
- Nonfinite inputs, unrepresentable multipliers or overflow of the resulting
  primary weight fail explicitly; they are not resampled or renormalized.

These changes repair arithmetic, not angular normalization. Energy, identity,
ancestry and time are preserved. Gaussian jets remain outside this port.
See [the original-contract verification](../validation/vmf_fork_contract.md)
for direct reference comparisons and regression scope. The earlier
[normalized-model checks](../validation/vmf_source.md) remain a historical
record and do not validate the restored weighted contract.
