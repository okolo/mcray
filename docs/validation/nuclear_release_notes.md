# v2.0.0: nuclear transport, YAML configuration and numerical corrections

The release scope covers nuclei through Fe56, free nucleons and the existing
electron/photon cascade. Nuclear functionality includes TALYS and PSB
photodisintegration, tabulated-mass decays, nuclear pair and photopion
processes, secondary nuclear photons, explicit missing-data compatibility
policies and YAML configuration with separate configuration/provenance files.
The merged correction to ICS and pair-secondary sampling is retained.
Nuclear processes and nuclear photon creation remain opt-in. The complete
YAML example explicitly enables its Fe56 settings; it is not a statement of
legacy defaults. Full-process runs require a separately prepared, compatible
data package; unsupported Strict-mode networks fail with diagnostics.

## Configuration and dependencies

Use `CRbeam simulation.yaml`, or keep using the compatible simulation CLI.
The complete commented schema is [config.full.yaml](../../examples/config.full.yaml).
Both interfaces write `config.yaml` and `run.yaml` instead of the old `params`
file. Consumers of that metadata file must migrate; particle output formats
are unchanged. Nuclear tables are installed in the standard `tables/nuclear`
location and verified by their manifests. They are prepared separately from
the pinned source data and are not embedded in the YAML or source archive.
The old `--talys-data` override remains available in CLI mode only.

Build with C++17 and the existing Boost/GSL/Fortran dependencies. CRbeam also
requires yaml-cpp >=0.6 and OpenSSL Crypto >=1.1.1 for SHA-256. CMake >=3.16 is
supported, including Ubuntu 22.04's CMake 3.22; the build script has a fallback
for versions lacking `--fresh`. See [YAML usage and data installation](../usage/yaml_configuration.md)
and [Linux build instructions](../../docker/README.md).

## Release decision and evidence

On 2026-10-26 the user limited remaining validation to one eight-hour server
check and selected release with documented limitations rather than further
rare-method development. That check ran on source
`18b28b51ed24945f3a9638fb55b804897159f1b6`. The release preparation adds
documentation and an offline audit tool; production C++/Fortran sources
are unchanged from that tested revision.

Run `20261026T210541Z-03385791` reached its wall-time limit (exit 124):
1017/1024 jobs and 7860369 complete primary histories were retained. Seven
jobs remain incomplete, with unknown consumption within their reservations;
they were not rerun. Build and all 155 Linux amd64 CTest entries passed.
All ten production coefficient controls matched. The archive SHA-256 is
`6d892bd176de81d8e0c367050c76aa0bbf9ecc15235039e1b3f36e54c7e6a892`;
all 4123 payload hashes and completed raw-history statistics were verified.

A subsequent offline audit verified all 1024 ensemble logs and found no
exception-driven retries or recovery failures. The exact tested constant-rate
path logs every retry directly to captured stderr; this mechanism can be
excluded as a cause of the retained mcray statistics (505 complete jobs,
3892433 primary histories). This does not remove the other acceptance limits
or establish the same diagnostic coverage for variable-rate production runs.

For the 24 ordinary quantities, 108 comparisons in complete groups meet the
existing agreement rule, 132 remain inconclusive and none meets its
disagreement rule. The 48 ordinary comparisons from incomplete groups remain
inconclusive. All 96 rare comparisons remain diagnostic. Three exploratory
single-antineutrino energy intervals nominally flag disagreement; their
interval coverage is unestablished and they are retained as unresolved
signals. No tolerance was relaxed and no uncertain result was relabelled
as successful validation. See [the complete retained-results report](release_check_8h_results_20261027.md).

Confirmed evidence is indexed in [the validation results](../README.md).
The latest [weighted pilot report](rare_transport_pilot_results_20261025.md)
does not establish final rare-channel acceptance. These isolated PD/decay ensembles do not by
themselves validate all GZK, PPP, EM-cascade or full-production combinations.

## Known limitations retained at release

- No claim of complete simultaneous 99% acceptance for all nuclear secondary
  outputs at the 1% integral / 5% spectral tolerances. The eight-hour matrix
  is incomplete; rare statistics and a simultaneous numerical-step bound
  remain unresolved. Weighted normalization alarms belong to experimental
  validation adapters, not to the ordinary production transport in this run.
- The 51-node beta-spectrum sampler retains its known discretization bias.
- Nuclear photons use the agreed CRPropa-compatible initial energetics;
  emitting nuclear/beta energies are not debited by photon energies, so exact
  event energy conservation is not claimed.
- CRPropa-compatible photon saturation and explicit transitions without
  photon tables are approximations requiring separate research. Such
  exceptions mean photons were not modeled, not that their physical yield
  is established to be zero.
- The finite TALYS grid/high-energy tail, data gaps reached by combined
  nuclear processes and extension beyond Fe56 remain separate work. Runtime
  Strict/CRPropa policies diagnose or record the approved data limitations.
  The GZK compatibility policy includes approximate residual completion and
  the PD compatibility policy skips missing PD tables only; see
  [GZK residuals and missing-data policies](../usage/nuclear_gzk.md).
- The official convolution endpoint omission is documented in the rate
  validation results; resolved discrepancies reach approximately 8.28% under
  the documented rate threshold. This is not a universal bound for weak
  rates. mcray is not retuned to reproduce that numerical omission. An
  English upstream report is prepared separately and has not been sent.

See [the eight-hour check contract](../../tests/crpropa/RELEASE_CHECK_8H.md)
for resource limits and interpretation. This is a release with documented
validation limits, not a claim that every planned physical acceptance gate
has passed. Release v2.0.0 is published from main, which incorporates the
prepared talys release through a pull request without changing production sources.
The exact release commit is identified by the v2.0.0 tag and the attached
validation-summary.json; the tested production-source revision remains
18b28b51ed24945f3a9638fb55b804897159f1b6.
