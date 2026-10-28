# Photons from nuclear processes

Two independent opt-in controls enable photon creation:
`--photodisintegration-photons On|Off` and `--decay-photons On|Off`.
Both default to Off. PD photons require a TALYS nuclear primary; decay photons
require `--nuclear-decay On`. PSB photon emission is not implemented.

Photon creation is independent of transport filtering. With `--noEMcascade`,
requested photons are created and then discarded by the existing filter.
Without that flag they enter the existing EM cascade. The effective creation
switches, yield policy and transport/filter state are recorded in `config.yaml`;
data provenance is recorded in `run.yaml`. Photon-free runs do not need photon payloads and retain their
existing sampling path. Enabling creation consumes random numbers even when
photons are later filtered, so individual nuclear trajectories can differ.

## Energy and lines

The initial approved approximation matches CRPropa: sample an isotropic
rest-frame line with energy E0 and boost using E=gamma*E0*(1-cos(theta)), with
the ultrarelativistic beta=1 approximation. PD uses the residual gamma for
the boost; decay uses the parent gamma before the compound decay sequence.
Photons inherit weight and ancestry. Their direction follows the existing
one-dimensional secondary transport convention.

PD residual energy and beta Q are not reduced by photon energy. Nuclear and
beta kinematics remain unchanged. Exact event energy conservation is therefore
not claimed; the photons add an explicit energy excess under this baseline.
Energy-conserving refinement is separate work. The beta sampler is unchanged.

Decay lines come from the pinned corrected generator's ground-state,
ionization-filtered gamma data. They are attached once to the selected
exclusive decay channel, before its beta/fragment sequence, rather than once
per emitted fragment or an inclusive parent branch. Each retained line has
an independent Bernoulli draw. An explicitly reviewed channel with no lines
is distinct from absent photon coverage; the latter fails loading. Invalid
line energies/probabilities are errors.

## Native PD yields and switchable policy

The same native background integral convolves the raw retained line cross
section and its parent/residual transition cross section. Conditional yield
is R_line/R_transition at the incoming nuclear gamma. The denominator is the
pinned `xs_photon_sum.txt`; exclusive nuclear channel selection retains the
existing independently normalized native channel rates. These separate
normalizations are preserved explicitly, not asserted to be interchangeable.
Interpolation and the 513-node zero-padding policy match existing TALYS data.

Choose `--pd-photon-yields Strict|CRPropa`:

| Policy | Conditional ratio above one | Positive line rate, zero transition rate |
| --- | --- | --- |
| Strict (mode 1) | Error | Error |
| CRPropa (mode 2, default) | One guaranteed photon of that line | One guaranteed photon of that line |

For valid ratios, both use independent Bernoulli emission per line. Zero
numerator and denominator give no photon. This is reference saturation, not
an inferred physical multiplicity model. The raw ratio remains available in
`PhotonYields` and the diagnostic probe. Strict errors occur when an affected
transition is evaluated; a small transport run can finish if it never visits
an inconsistent transition. Strict data errors bypass the engine's generic event-retry mechanism; they
fail the affected history and make the run invalid. Failed runs and partial
output are invalid.

## Explicit missing-source exceptions

The following ten approved transitions suppress photon creation, matching the
reference's absent-entry behavior. The list is stored in
`pd_photon_exceptions.tsv` and copied into effective metadata:

| Parent | Residual |
| --- | --- |
| D (2,1) | p (1,1) |
| T (3,1) | p (1,1) |
| T (3,1) | D (2,1) |
| He3 (3,2) | p (1,1) |
| He3 (3,2) | D (2,1) |
| He4 (4,2) | D (2,1) |
| He4 (4,2) | He3 (3,2) |
| Li7 (7,3) | He6 (6,2) |
| Li8 (8,3) | Li7 (7,3) |
| B8 (8,5) | Be7 (7,4) |

The first eight belong to the original Fe56 PD+decay network; the last two
are additionally reachable when nuclear GZK is enabled. Existing bundles
containing only the first eight remain valid for their supported networks.
Prepare a new photon bundle for the expanded network; do not edit an immutable
bundle or its hashes in place.

The residual label does not list emitted fragments. Nuclear interactions and
products remain enabled. These exceptions do not establish physical zero
photon production; their omitted secondary photons require investigation.
Other missing reachable photon transitions are errors. An explicit source
transition with no retained lines instead represents the supplied thinned
photon model. Neither policy changes the nuclear PD terminal whitelist or
completes missing GZK/PD masses and cross sections.

## Prepare and verify data

Run from the repository root, using a new output directory:

```sh
python3 tools/prepare_nuclear_photon_data.py \
  --repo /path/to/pinned/CRPropa3-data --data /path/to/verified-decay-bundle \
  --output /path/to/photon-bundle
python3 tools/prepare_nuclear_photon_data.py \
  --verify --output /path/to/photon-bundle
```

The source revision is `34b10e27375c0e9894d7c984ac15fb96e12975ef`.
The tool preserves existing nuclear tables and original manifests. It adds
raw PD transitions/lines, explicit exceptions, reviewed decay-channel
coverage and lines, photon provenance and SHA-256 hashes. Numerical payloads
remain outside Git. See [decay preparation](nuclear_decay.md) first.

For a switchable GZK bundle, pass `--decay-reference` pointing to the verified
corrected decay bundle containing `decay_sources/corrected/data/nuclear_decay.txt`.
GZK network and PD data limitations still apply. A PD-only bundle can be
extended without decay files; decay photons then remain unavailable.

After a Release build, run from `bin`:

```sh
../build/Release/CRbeam --particle nucleus --A 56 --Z 26 --z 0.00001 \
  --emax 5.6e20 --nparticles 100 --photodisintegration TALYS \
  --nuclear-decay On --photodisintegration-photons On --decay-photons On \
  --pd-photon-yields CRPropa --talys-data /path/to/photon-bundle \
  --output /path/to/new-output
```

## Short diagnostic tools

`nuclear_scenario` accepts the same three photon controls with processes
`talys`, `decay` and `talys-decay`; `nuclear-gzk` also accepts decay photons.
It retains photons as free output particles with EM transport disabled and
supports fixed seeds, primary histories, `--particle-output` and `--preflight`.

The isolated probe does not require the full descendant network to close:

```sh
build/Release/tests/nuclear_photon_probe /path/to/photon-bundle \
  56 26 1e19 CRPropa > /tmp/fe56-photon-yields.tsv
```

It reports rest-frame line energies, raw conditional ratios and effective
probabilities for each distinct transition, using native CMB at z=0. Values
are conditional; summing across different transitions is not a mean yield per
interaction. It does not run transport or weight ratios by transition frequency.
Nonzero exit status marks incomplete diagnostic output.

See [confirmed short checks](../validation/nuclear_photon_results.md).
Large-ensemble spectral agreement, saturation effects and the contribution of
missing light-transition photons remain unestablished.
