# Nuclear GZK and residual policies

The existing GZK process supports opt-in bound-nucleon photopion interactions
for A<=56. Free proton/neutron sampling uses the existing path. PSB masses
are excluded. Nuclear GZK is Off by default.

## Controls

Enable `--nuclear-gzk On` and select
`--nuclear-gzk-residuals Strict|CRPropa` (default Strict).
Both policies use the same prepared bundle:

| Policy | Pure p/n spectator | Other missing masses/decays |
| --- | --- | --- |
| Strict (mode 1) | Immediately split into free nucleons | Reject incomplete reachable data |
| CRPropa (mode 2) | Immediately split into free nucleons | Use explicitly tagged approximate masses and artificial dripping |

CRPropa requires `--nuclear-decay On` and a bundle marked
`gzk_compatibility=true`. Its missing-mass estimate is A*amu-Z*me, using
pinned source constants. An unknown decay gets mean lifetime 1e-99 s and
proton emission when Z>N, neutron emission otherwise. These are completion
approximations, not measured masses or lifetimes. Known corrected decay
channels and explicit stable records take precedence. Unsupported source
modes still fail. Strict never loads approximate mass rows or artificial
fallback channels.

Pure spectators, including (2,0) and (2,2), split in both policies without
creating an artificial bound particle or requiring a tabulated mass. Their
energy is divided equally among their free nucleons and their weights are
inherited. This compresses the reference's effectively instantaneous pure
spectator decay sequence; it is not finite-lifetime transport of a bound state.

The selected policy, mass model and provenance are recorded in `params`.
The internal TALYS mass enum identifies bundle-supplied numeric masses;
`CRPropaApprox` source tags and the effective `CRPropa-compatible-bundle`
metadata distinguish approximations from tabulated measurements.

## Data and commands

First prepare and verify a corrected decay bundle as described in
[nuclear_decay.md](nuclear_decay.md). Then create a new switchable bundle:

```sh
python3 tools/prepare_gzk_data.py --repo /path/to/pinned/CRPropa3-data \
  --decay-data /path/to/corrected-decay-bundle --output /path/to/gzk-bundle
```

The source revision must be `34b10e27375c0e9894d7c984ac15fb96e12975ef`.
The tool preserves measured mass precedence, the corrected exclusive decay
inventory, raw PD cross sections and the original PD terminal whitelist.
It records input/payload hashes and keeps artificial decays in a separate
file. It never overwrites an existing directory.

Run these commands from `bin` after a Release build, choosing a new output
path for each run:

```sh
../build/Release/CRbeam --particle nucleus --A 4 --Z 2 --z 0.00001 \
  --emax 4e20 --nparticles 100 --background 0 --noEMcascade \
  --photodisintegration Off --nuclear-decay On --nuclear-gzk On \
  --nuclear-gzk-residuals Strict --talys-data /path/to/gzk-bundle \
  --output /path/to/strict-output
```

Repeat with `--nuclear-gzk-residuals CRPropa` and a different output path.
Keep source, energy, background, seed and numerical controls identical when
assessing the policy effect. Identical seeds do not imply identical
trajectories after different decay paths become available.

Reachability is checked before creating output. The pinned bundle does not
close the Strict Fe56 GZK/decay network: Li4 (A4Z3) has no verified decay
record. That run fails explicitly; completing it requires additional verified
data. CRPropa closes that network with photodisintegration Off. Enabling
TALYS additionally requires actual cross sections for every reachable state;
missing PD data remain errors in both policies. The terminal whitelist is
not expanded by selecting CRPropa.

For reproducible isolated controls with an explicit seed, run from the
repository root (the output files must be new):

```sh
for policy in Strict CRPropa; do
  build/Release/tests/nuclear_scenario --process nuclear-gzk \
    --a 4 --charge 2 --gzk-residuals "$policy" \
    --talys-data /path/to/gzk-bundle --energy-per-nucleon 1e20 \
    --distance-mpc 0.1 --primaries 100 --batches 1 --seed 20261009 \
    --output "/tmp/gzk-${policy}.tsv"
done
```

This adapter isolates GZK plus nuclear decays. The Release CRbeam CLI retains
its existing time-derived RNG seed; use the adapter for fixed-seed diagnostics.
For roots where Strict cannot close, the useful result is a coverage error,
not a second propagated spectrum.

## Rates and kinematics

At equal gamma, proton/neutron rates are multiplied by 0.85*X^(2/3) for
2<=A<=8 and by 0.85*X for A>8, with X=Z or A-Z. Select the target by its
weighted rate. Existing SOPHIA sampling receives energy E/A; the residual
receives E*(A-1)/A. Products keep their parent weight and ancestry. Clones
share immutable nuclear data and clone the background integral.

This reproduces the selected bound-nucleon prescription, not every CRPropa
numerical detail: backgrounds and convolution remain native mcray, free
nucleon masses are unchanged, and existing SOPHIA and decay approximations
remain. Exact beta-decay event energy conservation is not imposed.

See [implementation diagnostics](../validation/nuclear_transport_baseline.md)
and [test tools](../../tests/README.md). PD and decay photons remain pending.
