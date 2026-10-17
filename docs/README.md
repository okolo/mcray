# mcray documentation

This repository owns implemented contracts, usage instructions, runnable test
adapters, and confirmed validation results and limitations. Shell examples
run from the mcray repository root unless explicitly stated otherwise.

## Usage and tests

- [Unified nuclei and TALYS](usage/talys_usage.md)
- [Nuclear decays and data preparation](usage/nuclear_decay.md)
- [Nuclear pair production](usage/nuclear_pair_production.md)
- [Nuclear GZK and residual policies](usage/nuclear_gzk.md)
- [Nuclear process photons](usage/nuclear_photons.md)
- [Tests and validation tools](../tests/README.md)
- [CRPropa adapters and external Linux runs](../tests/crpropa/README.md)
- [Docker Linux amd64](../docker/README.md)

## Validation results

- [PSB and nucleon results](validation/psb_validation_results.md)
- [SOPHIA table generation](validation/sophia_table_generation.md)
- [Nuclear decay / TALYS diagnostics](validation/nuclear_decay_comparison_results.md)
- [Nuclear transport regression baseline](validation/nuclear_transport_baseline.md)
- [Nuclear photon diagnostics](validation/nuclear_photon_results.md)
- [Photodisintegration physical-validation diagnostics](validation/photodisintegration_physical_validation.md)
- [Native Linux transport and secondary-photon pilots](validation/nuclear_transport_photon_pilot_results.md)
- [Deterministic photodisintegration reference](validation/deterministic_pd_reference.md)
- [Deterministic PD+decay energy-state reference](validation/deterministic_coupled_reference.md)

Results distinguish confirmed checks from remaining limitations and pilot
observations. Read their version and configuration scope before reuse.

## Private development workspace

Design proposals, research, decision history, validation plans and the future
task list are maintained in the private
[mcray_dev workspace](https://github.com/okolo/mcray_dev). Authorized users
can read its [documentation index](https://github.com/okolo/mcray_dev/blob/main/docs/README.md)
and [future tasks](https://github.com/okolo/mcray_dev/blob/main/docs/tasks/FUTURE_TASKS.md).
These links provide development context; public builds, tests and usage do
not require access to that repository. Do not duplicate its maintained plans
here. Keep README.md and AGENTS.md as the standard repository entry points.
