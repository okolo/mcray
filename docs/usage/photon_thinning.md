# Photon-threshold thinning

CRbeam keeps its existing energy-based thinning by default. The optional
`photon` policy replaces it; the two policies are never applied sequentially.
It acts on secondary particles produced by interactions, including lists with
only one secondary. It does not thin source primaries or change enabled
physical interactions, propagation filters or output columns.

For a photon with energy E below the threshold Et, retain it with probability
p=(E/Et)^alpha and replace its weight w by w/p. Other species and photons with
E>=Et retain their original weight. Alpha=0 disables this policy without RNG
consumption. No random numbers are consumed for nonphotons or photons at or
above threshold. Secondary order and all metadata except retained weights are
preserved. This is the fork's photon-thinning probability rule, using the
modern RNG and processes; it does not promise fork-identical trajectories.

## CLI and YAML

```sh
CRbeam --particle electron --z 0.01 --emax 1e13 \
  --thinning-model photon --thinning 1 --thinning-threshold-eV 1e12 \
  --seed 123 --output photon-run
```

The new long options are `--thinning-model energy|photon` and
`--thinning-threshold-eV VALUE`. Existing `--thinning` / `-t` selects alpha;
its default remains 0.9. The threshold is mandatory for `photon` and forbidden
for `energy`, even when alpha=0. Alpha must be finite and in [0,10] for photon
thinning or [0,1] for energy thinning. The threshold must be finite, positive
and representable after conversion from eV to internal units. Explicit legacy
CLI alpha values are now checked for numeric syntax and the existing energy
policy range in Release as well; valid legacy settings retain their behavior.

YAML uses schema version 2:

```yaml
schema_version: 2
source:
  redshift: 0.01
  particle: {type: electron}
  spectrum: {type: monoenergetic, energy_eV: 1.0e13}
simulation: {primaries: 100, seed: 123}
sampling:
  thinning:
    model: photon
    alpha: 1
    threshold_eV: 1.0e12
output: {directory: photon-run}
```

Omitted model means `energy`; a supplied threshold alone does not enable
photon thinning. Schema version 1 keeps its existing alpha-only mapping and
does not accept model/threshold keys. The effective version-2 config.yaml
records the model, alpha and applicable threshold. run.yaml records its hash
and the usual execution status; settings are not duplicated there.

## Weights and failure behavior

For one independent contribution, the thinned estimator is I*w/p, where I
is Bernoulli(p). Its conditional mean is w and variance is w^2*(1/p-1).
Apply the same relation to w times any linear observable. Thinning preserves
these expectations, not event-by-event multiplicity or energy. Large weight
variance can make a simulation statistically inefficient. Estimate errors
from complete independent primary histories, not correlated daughter rows.
Finite RNG resolution and finite sample size limit extremely rare retention
probabilities; no finite ensemble establishes a general physical acceptance.

The implementation evaluates the power through a logarithm of the energy
ratio, with a difference of logs when division underflows. It does not
multiply separate powers of dimensional energies. Active photon thinning
rejects nonpositive/nonfinite photon energies and negative/nonfinite photon
weights. Zero weights are allowed. An unrepresentable positive probability
or corrected weight is a fatal, non-retryable error, checked before the
Bernoulli draw. No probability floor or silent weight clipping is used.
A failed history increments the engine failure count; CRbeam fails the run
after the current batch and marks run.yaml failed. Any partial output from
that run is invalid. Alpha=0 returns immediately without validating particles.

The policy is immutable and can be shared by worker engines using their own
RNGs. Stable in-place compaction has O(n) time and O(1) auxiliary storage.
Changing policy changes RNG consumption and histories; retaining the existing
policy preserves its original algorithm and RNG behavior.

## Short checks

`unit.thinning.*` covers constructor/input errors, RNG bypass cases,
single-photon selection, mixed-list compaction and metadata, wide numeric
ranges, and engine failure handling for both propagation-rate paths.
The fixed software moment check uses 100000 independent one-photon histories
for each p=0.1,0.5,0.9, incoming weight 3 and seed 73102. It compares the first
and second raw moments with their analytic values using six analytic standard
errors for each, without adaptive reruns or tolerance changes. For X=I*w/p,
Var(X^2)=w^4*(1-p)/p^3. This checks the sampler/weight rule, not an astrophysical
acceptance threshold. CLI/YAML/snapshot and legacy-regression checks run in
the configuration suite. All are included in both Release and ASan+UBSan CI.
