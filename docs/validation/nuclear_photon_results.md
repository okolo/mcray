# Nuclear photon implementation diagnostics

Date: 2026-10-11. Source baseline: `f8b0c2aab79127466ddc16a942a59b146a03028c`,
branch `talys`, plus the recorded uncommitted photon implementation patch.
Nuclear data: CRPropa3-data `34b10e27375c0e9894d7c984ac15fb96e12975ef`;
runtime prescription inspected at CRPropa3
`73fa865b3304113da2059755a83f5939c1039d68`.

## Regression and deterministic checks

All 125 CTest entries pass in the native macOS arm64 Release build (29.26 s)
and the standard Docker Linux amd64 build under emulation (96.56 s).
Ten entries were added to the 115-entry suite: nine photon unit cases and
one CLI case. Coverage includes conditional native convolution, default
saturation versus strict errors, zero denominators, independent Bernoulli
multiplicities, boosts, weights, ancestry, clones, missing coverage, approved
exceptions, compound-decay association and the photon-free RNG path.

Synthetic line sampling uses five-standard-error checks for counts and
isotropic boosted energy means. Nuclear energy remains unchanged by photon
creation under the approved approximation; exact total event conservation
including photons is deliberately not asserted. Existing beta sampling and
nuclear process regression checks remain passing. The Docker result is not
a claim of execution on physical Linux amd64 hardware.

## Prepared real-data controls

The immutable `nuclear-photons-fe56-v2-20261011` bundle extends the verified
corrected Fe56 decay bundle. Original tables/manifests are retained; added
raw photon payloads and corrected-reference outputs have SHA-256 manifests.
The eight approved light-transition exclusions are explicit data and metadata,
not silent missing-data substitutions. The earlier preliminary photon bundle
without those exceptions is not the final control dataset.

Fe56 PD photon preflight succeeds with these exclusions. A native CMB control
with E/A=1e19 eV, L=0.1 Mpc, 20 histories, one batch and seed 20261011 completes
in CRPropa mode and produces 19 photon rows. Strict reports a fatal data
error on an inconsistent selected transition; its partial output is invalid.
An initial diagnostic exposed that the legacy engine retried generic sampling
exceptions by resampling channels. Strict data errors now explicitly bypass
that recovery path in both constant-rate and variable-rate stepping, and an
engine regression test checks that only one sampling attempt occurs. The
earlier apparently successful Strict output is discarded as biased.
TALYS plus nuclear decay with both photon sources enabled also completes all
20 histories and produces 19 photon rows. The same combined control completes
on Linux amd64 with the same species counts. These small samples do not
estimate the astrophysical difference between policies.

A photon-enabled B12 decay control uses E/A=1e13 eV, L=1e-9 Mpc,
10,000 complete histories, one batch and seed 20261011, without EM transport
or redshift losses. It produces 132 photons on both platforms. The corrected
branch-weighted expectation is 117.514 with count standard deviation 10.776;
the observed count passes a five-standard-deviation diagnostic check. Mean
photon energy is 4.7980e10 eV. This is a selected-line implementation control,
not broad acceptance of all nuclear radiation data.

Two-thread CRbeam B12 controls complete with EM transport enabled and with
`--noEMcascade`. The first retains a photon output file; the second creates
requested decay photons and filters them out. Creation/transport settings
are recorded separately. A further two-thread Fe56 TALYS run with photon creation and filtering also
completes. These are cascade-wiring smoke checks, not a
propagated EM spectral comparison.

## Conditional-yield probe

The isolated Fe56 probe convolves all 60 retained lines across its distinct
transitions on native CMB at z=0. It does not weight lines by the probability
of selecting their transition in transport:

| E/A (eV) | Lines with raw ratio >1 | Maximum raw ratio |
| --- | --- | --- |
| 1e18 | 28 | 3.2244e22 |
| 1e19 | 16 | 152.2643 |
| 1e20 | 1 | 1.3354 |

Strict correctly exits with an error on each scan; its partial diagnostic
output is incomplete. CRPropa completes and saturates the reported ratios.
Large conditional ratios can involve very small denominators. These numbers
are not unconditional photon yields, event frequencies or propagated effects,
and do not prove that saturation is physically justified. Source support,
thinning, denominator consistency and interpolation require investigation.

## Reproduction and limits

Public preparation commands, both policies, the exception list and probe
usage are in [nuclear_photons.md](../usage/nuclear_photons.md). The scenario
adapter supports fixed seeds and particle rows; numerical data are outside
Git. Source patch/file hashes, input/payload manifests, commands, logs,
container identity and output summaries are retained in the ignored external
validation directory under `nuclear-photons-implementation-20261011.json`.

No large acceptance ensembles or new SOPHIA tables were generated. Full
CRPropa photon-spectrum agreement, energy-conserving kinematics, a replacement
for saturation and the effect of secondary photons omitted from the eight
light transitions remain unestablished. Changes to their physical treatment
require separate approval.
