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

- Install required libraries to [external](src/external) folder
<pre><code>
cd bin
cmake -S ../src/app/crbeam -B .
make
</code></pre>
### Experimental isotope interactions

`Interactions::PD` is a library interaction for the explicit isotopes in
`ParticleType` (`H2` through `Fe61`). Add it to a `PropagationEngine` with
`engine.AddInteraction(new Interactions::PD(backgroundIntegral));`.
The CRbeam CLI currently accepts only light primary particles; this branch does
not yet provide a validated nuclear-propagation CLI.

PD reads `PD/eps.txt`, `PD/xs_pd_sum.txt` and `PD/xs_pd_thin.txt` below
`Utils::tables_dir`, using the same installation-relative lookup as other tables.
These TALYS data are **not distributed in this repository**. Supply a validated
export and retain its TALYS version, calculation settings and provenance with it;
the regression fixtures are synthetic and must not be used as physical data.
The loader expects whitespace-separated text (blank lines and `#` comments are
allowed):

- `eps.txt`: one strictly increasing positive rest-frame photon energy in MeV
  per row, at least two rows.
- `xs_pd_sum.txt`: `Z N sigma_1 ... sigma_K`, cross sections in millibarns on
  the `eps.txt` grid; one row per isotope.
- `xs_pd_thin.txt`: `Z N channel sigma_1 ... sigma_K`. The six decimal digits
  of `channel` count emitted neutrons, protons, deuterons, tritons, He3 and He4,
  respectively (leading zeroes may be omitted).

Supply totals for every propagated isotope, including light fragments; an
explicit all-zero row represents an intentionally disabled interaction. Missing
isotopes are errors rather than implicit zero rates. Positive total cross
sections require positive channel support. Thinned channel weights are
normalized at the sampled energy; their sum need not equal the total table.
Cross sections are interpolated linearly on a logarithmic `s` axis so zeroes
remain exact. Tables and interpolators are cached and shared by interaction
clones.

Unsupported residual isotopes are reported as model errors. PD validates them
at load time; a GZK residual outside the catalogue stops propagation rather than
silently dropping statistical weight or resampling the event. Model errors
propagate out of worker threads and the executable exits unsuccessfully;
partial output from a failed run must not be treated as a completed simulation.

The nuclear implementation remains an approximation: masses are sums of free
nucleon masses, nuclear decays are not included, and PPP retains the branch's
independent-proton (`Z` times the proton rate at the same Lorentz factor) model.
`ProtonPPcel` applies the matching nuclear fractional energy loss and should
be registered alongside PPP. This is not a validated coherent nuclear
pair-production model. The numerical checks below do not validate TALYS cross sections or these physical assumptions.

### Nuclear regression checks

With the normal build dependencies available, configure the library with
`cmake -S src/lib -B build -DBUILD_TESTING=ON`, build with
`cmake --build build`, then run `ctest --test-dir build --output-on-failure`.
Repeat with `-DCMAKE_BUILD_TYPE=Debug` and `Release` in separate build directories.
The tests create temporary synthetic PD tables and check more than 50 channels,
channel probabilities, residual nucleons, charge/energy/nucleon conservation,
invalid data, cached clones, GZK rate-weighted sampling, PPP energy scaling and
propagation of fatal model errors from worker threads.
