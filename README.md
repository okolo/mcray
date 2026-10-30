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

On Ubuntu 22.04/24.04 amd64, build the application, run the short tests and
install the complete pinned nuclear table package with one command from a
source checkout:

```bash
bash bin/install.sh --install-deps
```

`--install-deps` uses `sudo apt-get` (or apt-get directly when already root)
to install the system build and Python dependencies. Build products, generated
data, Python environment and logs remain in the checkout; no system-wide
CRbeam executable is installed. Use `--jobs 4` to change build/test parallelism
from its default of two. When dependencies are already installed, omit
`--install-deps`; `bash bin/install.sh --help` lists the options.

For a clean Ubuntu machine without Git, installation including source download
can be started with this single shell command:

```bash
sudo apt-get update && sudo apt-get install -y git ca-certificates && \
  git clone https://github.com/okolo/mcray.git && \
  bash mcray/bin/install.sh --install-deps
```

The installer builds the checked-out application revision. It fetches
CRPropa3-data at `34b10e27375c0e9894d7c984ac15fb96e12975ef`, then prepares
TALYS cross sections, corrected decays, the switchable GZK data and both nuclear
photon sources. No TALYS executable or CRPropa application build is needed.
Prepared data and their SHA-256 manifests are installed in `bin/tables/nuclear`,
alongside the other application tables. Python preparation uses an isolated
environment with NumPy 2.2.6, SciPy 1.15.3 and GitPython 3.1.45. Internet access
is needed for Ubuntu packages, GoogleTest, Python wheels and reference data.

Every attempt retains its log and preparation files under
`.crbeam-install/run.*/`. Failures stop installation before a new nuclear
package is published. Repeating the command verifies and preserves an existing
package with matching GZK/photon preparer fingerprints, without downloading or
regenerating its data. A partial, corrupt or incompatible existing package causes an error;
the installer never deletes or replaces it. Use a fresh checkout when changing
an installation whose data you need to preserve. Other table directories,
including `sophia2`, are unchanged.
Stage progress is printed to the terminal; detailed build, test and preparation
output is in the displayed `install.log` file.

After installation, run from `bin` so the common table-location mechanism finds
all installed tables:

```bash
cd bin
../build/Release/CRbeam /absolute/path/to/simulation.yaml
```

See [the complete commented YAML example](examples/config.full.yaml) and
[configuration and run records](docs/usage/yaml_configuration.md). Nuclear
processes remain opt-in; installing their tables does not enable them. The
data retain the documented compatibility approximations and validation limits.

For a manual build, install C++17, CMake >=3.16, GFortran, GSL, Boost.Filesystem/
System, yaml-cpp >=0.6 and OpenSSL Crypto >=1.1.1, then run
`bash bin/build.sh Release`. Test builds use pinned GoogleTest 1.18.0; see
[tests/README.md](tests/README.md).

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
