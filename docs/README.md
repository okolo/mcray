# mcray documentation

This repository owns implemented contracts, usage instructions, runnable test
adapters, and confirmed validation results and limitations. Shell examples
run from the mcray repository root unless explicitly stated otherwise.

## Usage and tests

- [Automated Linux build and short tests](usage/ci.md)
- [Source geometry, sphere observers and spatial output](usage/spatial_source_observer.md)
- [Truncated vMF angular source](usage/vmf_source.md)
- [Grid magnetic fields, composition and CLI/YAML input](usage/grid_magnetic_fields.md)
- [YAML configuration, CLI compatibility and run records](usage/yaml_configuration.md)
- [Optional photon-threshold thinning](usage/photon_thinning.md)
- [Complete commented configuration](../examples/config.full.yaml)
- [Unified nuclei and TALYS](usage/talys_usage.md)
- [Nuclear decays and data preparation](usage/nuclear_decay.md)
- [Nuclear pair production](usage/nuclear_pair_production.md)
- [Nuclear GZK and residual policies](usage/nuclear_gzk.md)
- [Nuclear process photons](usage/nuclear_photons.md)
- [Tests and validation tools](../tests/README.md)
- [CRPropa adapters and external Linux runs](../tests/crpropa/README.md)
- [Frozen PD/decay/photon acceptance and multicore server commands](../tests/crpropa/ACCEPTANCE.md)
- [Transport/photon v2 preparation and rare-channel tools](../tests/crpropa/ACCEPTANCE_V2.md)
- [Full weighted Fe56 adapters and independent server pilots](../tests/crpropa/RARE_TRANSPORT_PILOTS.md)
- [One-shot eight-hour release check](../tests/crpropa/RELEASE_CHECK_8H.md)
- [Retained eight-hour release-check results, 2026-10-27](validation/release_check_8h_results_20261027.md)
- [Nuclear transport release notes and known limitations](validation/nuclear_release_notes.md)
- [Docker Linux amd64](../docker/README.md)

## Validation results

- [Direct historical-fork field and cone comparison](validation/fork_kernel_comparison.md)

- [OpenMP worker lifetime regression](validation/worker_lifetime.md)
- [Photon-thinning verification](validation/photon_thinning.md)
- [Truncated vMF verification](validation/vmf_source.md)
- [Fork integration review against current main](validation/fork_integration_review.md)
- [Stable legacy mc_rotate delay checks](validation/mc_rotate.md)
- [Spatial source, observer and output checks](validation/spatial_source_observer.md)
- [Grid magnetic field implementation checks](validation/grid_magnetic_fields.md)
- [Grid magnetic field CLI and YAML checks](validation/grid_magnetic_field_input.md)
- [Corrected ICS/GammaPP spectra with fast invariant selection](validation/em_secondary_sampling.md)

- [YAML configuration regression checks](validation/yaml_configuration.md)
- [PSB and nucleon results](validation/psb_validation_results.md)
- [SOPHIA table generation](validation/sophia_table_generation.md)
- [Nuclear decay / TALYS diagnostics](validation/nuclear_decay_comparison_results.md)
- [Nuclear transport regression baseline](validation/nuclear_transport_baseline.md)
- [Nuclear photon diagnostics](validation/nuclear_photon_results.md)
- [Photodisintegration physical-validation diagnostics](validation/photodisintegration_physical_validation.md)
- [Native Linux transport and secondary-photon pilots](validation/nuclear_transport_photon_pilot_results.md)
- [Deterministic photodisintegration reference](validation/deterministic_pd_reference.md)
- [Deterministic PD+decay energy-state reference](validation/deterministic_coupled_reference.md)
- [Native Linux full PD+decay numerical matrix](validation/deterministic_coupled_full_matrix.md)
- [Weighted transport adapter preflight](validation/weighted_transport_preflight.md)
- [Independent weighted transport pilot results, 2026-10-25](validation/rare_transport_pilot_results_20261025.md)
- [Completed frozen transport/photon acceptance](validation/transport_photon_acceptance_results.md)
- [Frozen acceptance publication checks](validation/transport_photon_acceptance_preflight.md)

Results distinguish confirmed checks from remaining limitations and pilot
observations. Read their version and configuration scope before reuse.
