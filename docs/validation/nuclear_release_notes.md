# Nuclear transport release notes: pending final eight-hour check

The release scope covers nuclei through Fe56, free nucleons and the existing
electron/photon cascade. Nuclear functionality includes TALYS and PSB
photodisintegration, tabulated-mass decays, nuclear pair and photopion
processes, secondary nuclear photons, explicit missing-data compatibility
policies and YAML configuration with separate configuration/provenance files.
The merged correction to ICS and pair-secondary sampling is retained.

## Release decision and evidence

On 2026-10-26 the user limited remaining validation to one eight-hour server
check and selected release with documented limitations rather than further
rare-method development. The check is prepared but has not yet run.
Record its exact source revision, archive hash and observed outcomes here
before publishing. An uncertain comparison must remain uncertain.

Confirmed evidence is indexed in [the validation results](../README.md).
The latest [weighted pilot report](rare_transport_pilot_results_20261025.md)
does not establish final rare-channel acceptance. The last server build
passed all 155 CTest entries. These isolated PD/decay ensembles do not by
themselves validate all GZK, PPP, EM-cascade or full-production combinations.

## Known limitations retained at release

- No claim of complete simultaneous 99% acceptance for all nuclear secondary
  outputs at the 1% integral / 5% spectral tolerances. Rare-electron statistics,
  weighted normalization alarms and CRPropa step diagnostics remain unresolved.
- The 51-node beta-spectrum sampler retains its known discretization bias.
- Nuclear photons use the agreed CRPropa-compatible initial energetics;
  emitting nuclear/beta energies are not debited by photon energies, so exact
  event energy conservation is not claimed.
- CRPropa-compatible photon saturation and explicit transitions without
  photon tables are approximations requiring separate research.
- The finite TALYS grid/high-energy tail, data gaps reached by combined
  nuclear processes and extension beyond Fe56 remain separate work. Runtime
  Strict/CRPropa policies diagnose or record the approved data limitations.
- The official convolution endpoint omission is documented in the rate
  validation results. It is not corrected by retuning mcray. Preparation of
  an English upstream report remains a publication task; no message is sent
  automatically by the release tools.

See [the eight-hour check contract](../../tests/crpropa/RELEASE_CHECK_8H.md)
for resource limits and interpretation. The release tag/version and final
result citation are to be recorded after the check; no release is published
by this preparation document.
