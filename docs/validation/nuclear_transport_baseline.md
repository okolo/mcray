# Nuclear transport regression baseline

Date: 2026-10-09. Source: `3cc383564dfd79cb1c79134d1f64854be1ea62ec`,
branch `talys`, with no tracked source changes at baseline capture.

## Native macOS arm64

An incremental Release rebuild with `cmake --build build/Release --parallel 2`
completed. `ctest --test-dir build/Release --output-on-failure` passed all 104
entries in 29.95 seconds: 72 unit, 28 integration and four E2E entries (two
also labelled smoke). Tests used their existing seeds, data and tolerances.

## Linux amd64

The standard `bash bin/test-docker.sh` Release build and short test run passed
all 104 entries in 98.47 seconds with GCC/G++/gfortran 13.3.0. It explicitly
selects linux/amd64 on an arm64 workstation. No physical Linux hardware result
is claimed. These results precede the nuclear PPP additions below.

## Scope

These checks establish the existing implementation baseline, including
proton PPP, free-nucleon GZK, nuclear identity, TALYS and photon-free nuclear
decays. They do not implement or validate nuclear PPP, nuclear GZK, PD photons
or decay photons, and do not establish numerical agreement with CRPropa.
No acceptance ensembles or SOPHIA table generation were launched.

## Nuclear PPP implementation checks

The subsequent native Release rebuild passes all 108 entries in 29.14 seconds.
Four new cases cover nuclear rate scaling at equal gamma, non-A*mp masses,
clones, unchanged parent identity/energy/weight, secondary ancestry and energy,
thresholds, invalid mass/energy/thinning and CLI selection errors. The existing
independent proton differential-kernel quadrature checks remain passing.
The nuclear implementation uses those unchanged kernels; these tests do not
establish absolute agreement of CEL and secondary first moments over all
energies and backgrounds.

Short real-data isolated Fe56 CEL and He4 secondary controls completed with
the pinned TALYS mass bundle. CRbeam Losses with `--noEMcascade` and Full with
two threads also completed. These are execution diagnostics, not acceptance
ensembles. The changed implementation also passes all 108 entries in Docker
Linux amd64 emulation in 99.59 seconds, using the standard Release build.

A larger short He4 secondary diagnostic used 100 complete primary histories,
E/A=1e19 eV, L=0.1 Mpc, seed 20261001, frozen native CMB and no EM transport
or redshift losses. It produced 62 retained lepton rows. Mean nuclear energy
loss was 2.9776e15 eV/history; mean weighted lepton energy was 2.5921e15
eV/history with standard error 0.6078e15 eV from independent primary histories.
This pilot does not establish equivalence or a precision target.


## Nuclear GZK and switchable residual checks

The subsequent macOS Release build passes all 115 CTest entries in 28.80 s.
The standard Docker Linux amd64 build also passes all 115 entries in 99.01 s
under arm64 host emulation. Six GZK cases and one CLI case were added to the
108-entry PPP suite. They
cover A=8/9 modifiers, equal-gamma rates, free p/n regression, clone behavior,
both policies on a common synthetic bundle, artificial decay activation,
pure p/n breakup, preflight and unsupported/missing-data errors. Two hundred
SOPHIA events check baryons, charge including leptons, weights and energy
(relative tolerance 1e-6); target frequencies use five binomial standard
errors. This does not establish full CRPropa spectral agreement.

A real switchable bundle was prepared from pinned data
34b10e27375c0e9894d7c984ac15fb96e12975ef and the hash-verified corrected decay
inventory. It contains 605 known channels, 1145 separate artificial channels
and 1096 separately tagged approximate masses (not all are reachable).
Known source modes with unspecified branch weights but valid lifetimes are
preserved. Zero-weight unsupported modes do not create channels; reachable
unsupported modes are errors in both policies.

He3 preflight succeeds in both modes (five reachable identities). Fe56
Strict preflight rejects missing Li4 decay data; CRPropa succeeds with PD Off
(836 identities, including pure spectators and free nucleons). Adding TALYS
to that CRPropa network rejects missing reachable cross sections. Real CRbeam
checks confirm both errors occur before output creation; the CRPropa PD-Off
case completes with two worker threads.

Short isolated transport uses native CMB, E/A=1e20 eV, seed 20261009, one
batch, no redshift loss, EM transport or free-neutron decay. For He4,
100 histories over 0.1 Mpc produce identical 112 physics rows in Strict and
CRPropa. Fe56 CRPropa completes 20 histories over 0.01 Mpc. A separate Li4
CRPropa control (four histories, E/A=1e19 eV, L=1e-6 Mpc, same seed) produces
four He3 and four protons, demonstrating artificial dripping. These are
execution controls, not precision estimates of astrophysical policy effects.

Source patch/file hashes, bundle manifests, commands, seeds, logs and outputs
are retained outside Git. See [usage](../usage/nuclear_gzk.md) for public
preparation and reproduction commands. PD/decay photons and physical
acceptance of the new processes remain pending.


## Photon implementation follow-up (2026-10-11)

PD and decay photon controls, native conditional convolution, selectable
Strict/CRPropa yields and explicit light-transition exclusions are now
implemented. All 125 entries pass on macOS and Docker Linux amd64 emulation.
See [photon diagnostics](nuclear_photon_results.md) for the independent dated
scope, short controls and remaining physical-validation limitations.

## Missing-PD compatibility and expanded photon coverage, 2026-10-20

Implemented on talys over source `7fc04ae6ce2e76887287ba85562b721d6b31ddf8`;
these checks used the uncommitted implementation patch SHA-256
`35ab10daddcb388fdac2ebd5a99007be0dba5bb819007b5349c9c3835fd88011`.
`--pd-missing-data Strict|CRPropa` now independently controls wholly absent
PD tables. Strict remains the default. Compatibility skips only PD for those
states; partial/malformed tables and missing masses/required decays still fail.
The original 13 explicit terminals are unchanged. Clones share the selected
immutable policy. The CLI records the policy and each reachable missing-table
identity in `params` in these pre-YAML controls. Current runs instead record
settings in `config.yaml` and provenance/exceptions in `run.yaml`, as documented
in the [YAML contract](../usage/yaml_configuration.md).

The combined Fe56 GZK/PD/decay graph contains 836 identities; 582 reachable
bound states lack PD tables and are explicitly recorded as skipped. This
includes states formerly rejected by preflight, not evidence that 582 states
were visited in these short histories. An exhaustive photon-coverage inventory
of that reachable graph finds exactly two additional absent transitions beyond
the original eight: Li8 -> Li7 and B8 -> Be7. Both are explicit approved
no-photon exceptions, matching CRPropa; their nuclear interactions/products
remain enabled. Other absent photon transitions still fail. Fresh photon-data
preparation includes all ten explicit exceptions and preserves their hashes.
Original eight-exception bundles remain supported for their original network.

Synthetic tests verify Strict errors, compatibility zero PD rates, unchanged
rates for covered states, clone behavior, rejection of partial/unloaded data,
continued GZK and artificial decays, metadata, independent CLI policies and
retained nuclear products/energy/weights for both new photon exceptions.
All 134 CTest checks pass on macOS arm64 (14.90 s) and Docker Linux amd64
emulation (47.68 s). These timings describe regression runs, not a performance
benchmark. All 21 deterministic Python controls also pass in the existing
NumPy/SciPy environment.

Short real-data CLI controls use Fe56, seed 20261020, total energy 5.6e19 eV,
source z=1e-8, CMB, two primary histories and two worker threads. Strict with
GZK residual compatibility rejects missing PD data before creating output.
Missing-PD compatibility completes a TALYS+decay+GZK run. With the newly
prepared ten-exception photon bundle, native and Linux controls also complete
TALYS + decays + GZK + PPP Full + both photon sources + the EM cascade.
Effective params confirm all requested switches and the 582-state exception
inventory. These tiny controls verify functional wiring and data preflight,
not interaction statistics, propagated spectra or physical acceptance.

Reproduce the full-process control using the command in
[the GZK contract](../usage/nuclear_gzk.md), omit `--noEMcascade`, use a fresh
photon-enriched GZK bundle and add `--nuclear-ppp Full
--photodisintegration-photons On --decay-photons On`. The prepared control
bundle's manifest SHA-256 values are:

- nuclear: `f1eb796059cc66e5351cdeb1eff025d12e8fa42757a893f543d4267f6249c608`;
- decay: `a32f343739fee7b01bd85d5944e8b90cc882c0143c0f1eae211b2f778e23beae`;
- photon: `01e4517adb191b11d567039d1e58485f4e0268f3603a9d4d007aae94de29a9ab`.

The [frozen initial-point acceptance workflow](../../tests/crpropa/ACCEPTANCE.md)
is separate: it retains the original piloted nuclear/photon input hashes,
uses independent deterministic M/C targets, and excludes GZK/PPP/EM to isolate
PD+decays and photon creation. Its 272-job plan is prepared without launching
new final histories on the workstation. Neither this functional control nor a
successful workload archive substitutes for its statistical verdict.
