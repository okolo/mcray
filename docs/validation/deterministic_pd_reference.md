# Deterministic PD-only reference: initial controls, 2026-10-17

The [PD+decay energy-state extension](deterministic_coupled_reference.md) is
now implemented and checked separately. This report retains the PD-only scope
and original controls below.

The first implementation provides independent branching transport and photon
first moments, plus I/M/C coefficient adapters described in the
[tool contract and commands](../../tests/deterministic/README.md). This is a
PD-only diagnostic, not completed physical acceptance or a combined decay
reference. Production source code and the physical model were not changed.

## Frozen scope

Production baseline: mcray 82df04adc32420fd68c12656e23873d7977d17fa on talys,
with the new uncommitted validation tools captured by per-run hashes. Reference:
CRPropa3 73fa865b3304113da2059755a83f5939c1039d68 and CRPropa3-data
34b10e27375c0e9894d7c984ac15fb96e12975ef. Inputs and official photon tables
are those of retained server pilot 20261016T102021Z-03eb6c57; see its
[audit and archive identity](nuclear_transport_photon_pilot_results.md).

Monoenergetic Fe56, E/A = 3.9810717055349857e18 eV, L = 1 physical Mpc,
frozen z=0 CMB or CMB+Inoue12. TALYS PD and PD photons enabled; all decays,
losses, thinning and EM propagation disabled. The PD-only closure has 160
states including free p/n. Fixed E/A uses the approved PD product rule exactly.
Photons retain saturated Bernoulli probabilities, residual-gamma boost, eight
explicit no-photon exceptions and no debit from nuclear energy.

## Initial deterministic means

Means per injected Fe56; I uses order 64 per raw cross-section interval, M uses
background step .0025 and eps-rel 1e-6, C uses unchanged official tables and
runtime rules. Photon energy is the total expected emitted energy per history.

| Background | Provider | Fe56 survival | Photon count | Photon energy, eV |
|---|---|---:|---:|---:|
| CMB | I | 0.567727467712 | 0.416874111584 | 1.716941908884e15 |
| CMB | M | 0.567734363944 | 0.416864183648 | 1.716900398459e15 |
| CMB | C | 0.558131124237 | 0.431165206947 | 1.775823758182e15 |
| CMB+Inoue12 | I | 0.524134041063 | 0.487268824661 | 1.985527522210e15 |
| CMB+Inoue12 | M | 0.524136665935 | 0.487264529175 | 1.985508590160e15 |
| CMB+Inoue12 | C | 0.515250396309 | 0.501280757159 | 2.044537110862e15 |

Relative M/I differences in these three means are at most 0.00242% for CMB and
0.000954% for CMB+Inoue12. C/I survival shifts are -1.6903% and -1.6949%;
photon-count shifts are +3.4282% and +2.8756%. These are deterministic
coefficient-representation differences, without Monte Carlo noise or an
empirical correction. They must be retained separately from transport sampling.
The previously identified high-gamma endpoint discrepancy is a separate effect;
these injection-point results do not measure that high-gamma limit.

## Numerical controls and limits

Thirteen short analytic/rule tests pass: exponential survival/rewards, stable
and zero-length controls, sequential equal/unequal rates, distinct branches,
fragment multiplicities, line/bin moments, conservation and malformed inputs,
two path solvers, full/exclusive normalization, C++ nearest-node rounding,
serialized branch-tail handling and saturated/zero-denominator line yields.

Native macOS arm64 exporter built through CMake with Homebrew Clang 22.1.5
and gfortran. The same adapter compiled against verified matching production
sources on Docker Linux amd64. This does not replace a native-server campaign.
Both backgrounds were evaluated by the Linux exporter; resulting mean rewards
agree with the native exporter to at most 5.6e-16 relatively.
Additional real-data controls give exactly zero photon rewards for Off with
each I/M/C provider, and unit survival/zero rewards for C outside its open
gamma grid.

Independent order 32 to 64 changes the three reported means by less than
2e-11 relatively on both backgrounds (CMB order 16 is also retained).
Sparse matrix exponential versus independent DOP853 path integration changes
the photon count by less than 3.1e-15 absolutely and total photon energy by less
than 15 eV. Both checks preserve expected A=56 and Z=26 to rounding precision.
Diagnostic photon bins [0,1e14,1e15,1e16,infinity) eV are integrated analytically;
they are not final acceptance bins. Refinement differences are observed
numerical controls at this injection point, not certified global error bounds.

The observed elapsed time of a coefficient build and solve was about 1.6 s
for C, 9 s for I/order32, 15 s for I/order64, and 10 s for M on the workstation.
No new primary histories were generated. All artifacts, tool/input hashes,
settings and logs are retained under the ignored deterministic-reference run
directory; public workflows do not depend on its local path.

## Retained pilot sampling versus each program's coefficients

The comparison adapter matched the PD-only, photons-On, finest-level pilot
batches (20 batches of 1000 histories per background and program). Payload
hashes were verified. Below are (pilot mean - deterministic mean) divided by
the exploratory batch standard error; mcray uses M and CRPropa uses C.

| Background | Program/control | Survival | Photon count | Photon energy |
|---|---|---:|---:|---:|
| CMB | mcray / M | +1.670 | -0.955 | -1.841 |
| CMB | CRPropa / C | -1.495 | -0.805 | -0.263 |
| CMB+Inoue12 | mcray / M | -0.545 | -0.863 | -0.766 |
| CMB+Inoue12 | CRPropa / C | +0.024 | -0.189 | +0.546 |

These diagnostics do not distinguish a 1% transport error at the available
pilot precision and are not simultaneous confidence intervals. They support
keeping coefficient differences separate from sampling/step checks, but do
not establish transport acceptance or certify rare isotope/photon yields.

Combined PD+decays requires energy-state remapping and its own convergence
checks. The current solver supplies means, not higher moments or history
probabilities. Rare decay photons, final spectral intervals, simultaneous
confidence construction and numerical envelopes remain outside these results.
The agreed simultaneous 99% confidence, 1% integrated / 5% spectral tolerances
and 200000-history limit have not changed. No final acceptance ensemble or
acceptance verdict follows from these controls.
