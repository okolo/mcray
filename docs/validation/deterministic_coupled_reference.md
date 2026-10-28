# Deterministic PD+decay initial controls, 2026-10-17

The [full native Linux numerical matrix](deterministic_coupled_full_matrix.md)
completed on 2026-10-20. Its audited results supersede the earlier pending
server status for this injection point; final Monte Carlo acceptance remains
open. Historical implementation checks and the first failed attempt follow.

The joint energy-state reference is implemented on talys. Its
[public contract, commands and external runner](../../tests/deterministic/COUPLED.md)
describe the actual reduction, numerical controls and limits. Production
physics sources are unchanged. This extends the
[PD-only controls](deterministic_pd_reference.md); it does not establish final
physical acceptance or freeze a statistical protocol.

## Versions and scope

The production baseline remains mcray 82df04adc32420fd68c12656e23873d7977d17fa,
with new validation tools captured by hashes. CRPropa3 is pinned to
73fa865b3304113da2059755a83f5939c1039d68; CRPropa3-data is pinned to
34b10e27375c0e9894d7c984ac15fb96e12975ef. Prepared raw data and official
tables are those used by retained server pilot 20261016T102021Z-03eb6c57;
its [audit](nuclear_transport_photon_pilot_results.md) records the archive identity.

Fe56, E/A=3.9810717055349857e18 eV, L=1 physical Mpc, frozen z=0 CMB or
CMB+Inoue12, TALYS PD and nuclear decays enabled. Both photon sources are On,
with saturated PD yields, the eight approved omissions and no energy debit.
Free-neutron decay, cosmology, losses, GZK, PPP, thinning and EM propagation
are Off, matching the isolated pilot fixtures. Corrected decay channels and
the 51-node beta approximation are retained, including the two programs'
existing constants/mass conventions.

## Energy-state and stiff-path checks

The 160-isotope closure is unchanged. All allowed nuclear paths imply
0.994279554757 <= (E/A)/(injected E/A) <= 1 in this network. Each isotope
has its own propagated bounds; four states have point intervals and 156
have finite intervals. N=5 gives 784 energy states; N=9 gives 1408.
Beta daughters retain gamma, compound channels execute in production order,
and positive linear-energy remapping preserves number and mean energy.
Rates are evaluated directly at every isotope/node energy, without a second
coefficient interpolation grid. Targets outside the bounds are errors.

All finite decay times are retained, including the very short-lived states.
BDF solves the sparse stiff system; independent Radau checks use tenfold
tighter numerical tolerances. No instantaneous-decay approximation is used.
At this pilot point, BDF versus Radau differences are below 2.4e-11 relative
for total photon count/energy and below 4.3e-10 for decay-photon count/energy.
These are observed numerical differences, not certified global error bounds.

I and M were compared at N=3 and 5; C at N=5 and 9, on both backgrounds.
Among Fe56 survival, total/decay photon number and energy, and the diagnostic
photon-bin moments, the largest observed grid-refinement change is below
8e-11 relatively. The bins [0,1e14,1e15,1e16,infinity) eV are analytic reward
integrals and remain diagnostic choices. These controls cover one injection
point and distance; they do not certify rare isotope yields, energy-state
distribution shapes, other injection energies, or long paths.

Expected baryon number remains 56, total charge including beta leptons remains
26, and (nuclear plus lepton mean energy)/(injected energy) remains one to
rounding precision in the retained runs. Photon energy is additional under
the accepted compatibility model and is excluded from that invariant. Tiny
negative energy-state populations from solver roundoff (around 1e-33 or
smaller in these runs) remain visible; they are not clipped into resolved
rare yields. Refinement reports omit relative precision claims for components
at/below the solver's coordinate tolerance.

## Initial means

Means per injected Fe56. I uses order16 per cross-section interval and N=5;
M uses native background step .0025, eps-rel 1e-6 and N=5; C uses the
unchanged official tables and N=9. All use BDF rtol=1e-10, atol=1e-16.

| Background | Provider | Fe56 survival | PD photon count | Decay photon count | Total photon energy, eV |
|---|---|---:|---:|---:|---:|
| CMB | I | 0.567727467712 | 0.416905757836 | 0.00144762700807 | 1.72083964342e15 |
| CMB | M | 0.567734363944 | 0.416895827448 | 0.00144754810551 | 1.72079789707e15 |
| CMB | C | 0.558131124238 | 0.431201354742 | 0.00155346417386 | 1.78002745975e15 |
| CMB+Inoue12 | I | 0.524134041058 | 0.487331951948 | 0.00224903514001 | 1.99184527314e15 |
| CMB+Inoue12 | M | 0.524136665936 | 0.487327656155 | 0.00224900327788 | 1.99182623814e15 |
| CMB+Inoue12 | C | 0.515250396311 | 0.501349963506 | 0.00235810084488 | 2.05118562411e15 |

M/I total photon number/energy differences are below 0.00243%; decay-photon
count differs by -0.00545% on CMB and -0.00142% on CMB+Inoue12. C/I total
photon-count shifts are +3.4424% and +2.8855%, while decay-photon counts
shift by +7.3111% and +4.8494%. These deterministic differences retain the
distinct official coefficient representation and conventions. They are not
Monte Carlo transport errors, and no empirical correction is applied.
They should not be identified with the earlier high-gamma endpoint limit.

## Retained pilots against their own coefficients

The comparison adapter selects matching complete `talys-decay` cases with
both photon flags On, finest pilot level, 20 independent batches of 1000
histories per program/background. Output hashes are verified. The table gives
(pilot mean - deterministic mean)/exploratory batch standard error, using M
for mcray and C for CRPropa.

| Background | Program/control | Survival | Photon count | Photon energy | Positron count |
|---|---|---:|---:|---:|---:|
| CMB | mcray / M | +0.576 | -0.357 | -0.582 | -0.657 |
| CMB | CRPropa / C | -0.181 | +0.465 | +1.665 | +0.258 |
| CMB+Inoue12 | mcray / M | +2.520 | -1.055 | +0.332 | -0.998 |
| CMB+Inoue12 | CRPropa / C | -0.000 | -0.720 | -0.712 | +0.042 |

The +2.52 survival residual is retained for investigation with the frozen
protocol; it is not discarded or treated as a standalone acceptance failure.
These batch residuals are not simultaneous 99% intervals and do not resolve
a 1% transport budget. Sparse decay photons remain statistically inconclusive
despite their precise deterministic mean. No additional histories were used.

Separate local coefficient controls keep the energy grids fixed: doubling I
order16 to order32 changes the reported main means by less than 1.2e-11
relatively. Halving the M background step from .0025 to .00125 changes total
photon number/energy by at most 1.82e-5 relative and decay-photon number/energy
by at most 4.13e-5 relative. The remaining M/I difference at this point is
therefore larger than the observed energy-grid effect and responds to native
background-integration refinement. These are diagnostics; production settings
and the physical model were not changed.

## Verification and remaining limits

Twenty-one Python controls pass. New controls cover the exact mean of the
51-node inverse-CDF approximation, compound beta/emission ordering, positive
number/energy-preserving remapping with overflow rejection, DAG energy bounds,
finite decay survival and parent-boosted lines, competing PD/decay rates,
a finite 1e24/Mpc rate at both resolved short and long distances, stiff BDF/Radau
agreement, and reduction to the independent PD-only matrix-exponential result.

The extended exporter built on macOS arm64 and Docker Linux amd64. On the
same CMB isotope/energy slice, exported PD coefficients agree to below 6.4e-16
relatively; total/partial decay rates and decay lines agree exactly. All 129
existing Linux CTest checks pass on the verified matching production library.
The Python joint solver has been checked locally; the full native Linux
refinement campaign remains to be run with the prepared wrapper.

Observed coefficient/solve/path-check times at N=5 were about 21.5 s for I
and 49-51 s for M on the workstation; C at N=9 took 7 s. These are diagnostic
elapsed times with concurrent background cases, not performance benchmarks.
Manifests, coefficient slices, kernels, grids, sparse matrices, populations,
refinements, pilot comparisons and logs are retained in the ignored
coupled-reference artifact directory. Public workflows do not require its
local path or additional documentation.

The external initial/full matrices separate energy-grid, independent
quadrature, native-background and stiff-path refinements. Wider injection/path
coverage, rare-output numerical envelopes, lepton spectra/correlations and
final statistical protocol freezing remain outside these initial results.
The physical model, simultaneous 99% confidence, 1% integrated / 5% spectral
budgets and 200000-history cap remain unchanged. No acceptance verdict follows
from these checks.

## Native-server launch failure audited, 2026-10-19

Run 20261017T155033Z-90561e31 fetched and built the requested public commit
7af65eda1afcd2470989ff7c139f199efe6108a1 successfully. Its archive SHA-256 is
`63300fc007c8b016307253dc004a864114007c251080af6c95da47eaf378c405`;
the archive and all ten recorded payload checksums were verified.

CTest passed 125 of 129 checks. Four existing
`ProtonPppSampling.MatchesIndependentRestFrameQuadrature` parameter cases
(k=3,5,20,50 electron masses) reached the per-test 30-second wall-clock limit.
The other two cases passed in 27.06 and 27.49 seconds. There were no reported
assertion failures in the timed-out cases; their numerical outcomes remain
unknown. Missing-table warnings are not the cause established by this log:
these cases construct their own fixed-rest-frame photon field. The archived
logs do not establish why the server timing varied relative to earlier runs.

The runner stopped at CTest, before starting any deterministic matrix case.
The archive therefore contains no full-matrix results and cannot establish
grid convergence or physical agreement. No acceptance histories were used.

The six PPP sampling cases now have a 180-second test timeout; all other
timeouts remain 30 seconds. Their quadrature orders, sample counts, fixed
seeds and numerical/statistical assertions are unchanged. The physical model,
transport acceptance budgets and history cap are unchanged. A new server run
must use the published fix and pass the complete CTest gate before the full
matrix starts; failed tests are not excluded or treated as passes.

Verification of this correction on Docker Linux amd64: all six PPP cases
passed in 80.16 seconds total (12.96-14.32 seconds each). CTest JSON discovery
contains exactly 129 unique names, six with TIMEOUT=180 and the remaining
123 with TIMEOUT=30. Production-library source checksums matched the previous
verified baseline. These local checks do not replace the required new native
server run or provide any full-matrix result.
