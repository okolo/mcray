# mcray
Framework for Monte Carlo simulation of ultra-high energy cosmic rays and electromagnetic cascade propagation.

### Authors:
   Oleg Kalashev and Alexander Korochkin

### Features
 - propagation of protons, neutrons, nuclei, electon-photon cascade and neutrino can be simulated
 - support for interactions with arbitray photon background including trajectory simulation for secondary particles, produced in the interactions
 - trajectories of individual particles in presence of magnetic field can be calculated

### Applications

[CRbeam](src/app/crbeam) - cosmic ray beam simulation

### Installation

The code requires C++17. Test builds require CMake 3.16 or newer and use
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

## Documentation

See the [documentation index](docs/README.md) for design, reference, usage and
validation documents, and the [future task list](https://github.com/okolo/mcray_dev/blob/main/docs/tasks/FUTURE_TASKS.md)
for deferred work.
