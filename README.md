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
- Install [GSL](https://packages.debian.org/sid/libgsl-dev) with headers: `sudo apt-get install libgsl-dev`
- Install required libraries to [external](src/external) folder
- Build code:
<pre><code>cd bin
cmake -S ../src/app/crbeam -B .
make
</code></pre>
