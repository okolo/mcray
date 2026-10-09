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
