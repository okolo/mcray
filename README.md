# mcray
Framework for Monte Carlo simulation of ultra-high energy cosmic rays and electromagnetic cascade propagation.

For v2.0.0-rc1 nuclear transport, YAML configuration, migration notes and known
validation limits, see the [release notes](docs/validation/nuclear_release_notes.md).

### Authors:
   Oleg Kalashev and Alexander Korochkin

### Features
 - propagation of protons, neutrons, nuclei, electon-photon cascade and neutrino can be simulated
 - support for interactions with arbitray photon background including trajectory simulation for secondary particles, produced in the interactions
 - trajectories of individual particles in presence of magnetic field can be calculated

### Applications

[CRbeam](src/app/crbeam) - cosmic ray beam simulation

### Installation

The code requires C++17. CRbeam additionally requires yaml-cpp >=0.6 and
OpenSSL Crypto >=1.1.1 for SHA-256. See [YAML configuration and table installation](docs/usage/yaml_configuration.md). Test builds require CMake 3.16 or newer and use
GoogleTest 1.18.0. See [tests/README.md](tests/README.md).

- Install required libraries to [external](src/external) folder
<pre><code>
cd bin
cmake -S ../src/app/crbeam -B .
make
</code></pre>

### Linux amd64 checks

The primary target platform is Linux amd64. Build and run short tests in Docker:

```bash
bash bin/test-docker.sh
```

For the environment, test coverage, and rerun instructions, see
[docker/README.md](docker/README.md).

### Code review

[Independent OpenCode PR review](docs/usage/opencode_review.md) uses GLM-5.3
through OpenRouter, automatically and on request with `/opencode`.

## Documentation

See the [documentation index](docs/README.md) for usage, reference tools and
validation results. The [release notes](docs/validation/nuclear_release_notes.md)
describe the supported scope and known limitations.
