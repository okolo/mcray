# Nuclear pair production

Nuclear PPP is opt-in for A=2..56 with tabulated nuclear masses. It uses the
existing native proton continuous-loss and differential secondary kernels
at equal Lorentz factor. It does not use CRPropa's secondary spectrum tables.

```sh
CRbeam --particle nucleus --A 56 --Z 26 --talys-data /path/to/data \
  --photodisintegration Off --nuclear-ppp Losses --noEMcascade --z 0.01
```

`--nuclear-ppp Off|Losses|Full` defaults to Off. Losses enables only nuclear
continuous energy loss; Full additionally samples e-/e+ products and requires
EM cascade transport. Full with `--noEMcascade` is rejected. Both active modes
require `--talys-data`, even for PSB-chain source isotopes, and reject PSB
photodisintegration. TALYS and nuclear decays can be selected independently
with complete reachable data. Heavy nuclear PPP is outside the current scope.

For gamma=E/M, evaluate the proton kernels at Ep=gamma*mp. The nuclear
relative loss rate is Z²*mp/M times the proton relative loss rate. The
secondary event rate is Z² times the proton rate; absolute lepton energies
at equal gamma follow the proton kernel. Only CEL reduces primary energy.
Secondary generation preserves A,Z, mass, energy and weight of the surviving
nucleus. Proton and free-neutron behavior is unchanged.

`--ppp-thinning` divides the secondary event rate and multiplies secondary
weights, as for proton PPP. One electron or positron is sampled per accepted
event with an additional factor of two in weight. This represents pairs in
expectation, rather than producing two particles in each event. The existing
native high-rest-energy approximation above k=20000 electron masses remains.
Exact event energy conservation is not asserted. Energy-loss/secondary first
moments need independent checks before physical acceptance.

The module consists of NuclearPPPcel and NuclearPPP, using cloned background
integrals and unchanged proton kernels. Particle ancestry and production time
are inherited through a proton surrogate; the engine assigns child identities.
Effective mode and scaling/energy ownership appear in the params file.

## Isolated diagnostics

The test-only scenario runner disables redshift energy loss and subsequent
EM transport, keeping complete primary histories. From the repository root:

```sh
build/Release/tests/nuclear_scenario --process nuclear-pair-cel \
  --talys-data /path/to/data --a 56 --charge 26 --distance-mpc 1 \
  --energy-per-nucleon 1e19 --primaries 20 --batches 1 --output /tmp/losses.tsv
build/Release/tests/nuclear_scenario --process nuclear-pair \
  --talys-data /path/to/data --a 4 --charge 2 --distance-mpc 0.01 \
  --energy-per-nucleon 1e19 --primaries 20 --batches 1 --output /tmp/pairs.tsv
```

`--scan` writes rates across the existing energy scan; `--photon-table` allows
the existing shared tabulated background. A mass-only bundle suffices for
these isolated PPP modes. Short checks do not establish CRPropa agreement,
light-nucleus validity, first-moment normalization or full transport acceptance.
