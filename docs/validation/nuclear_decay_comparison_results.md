# Nuclear decay and TALYS comparison diagnostics

Date: 2026-10-08. Production baseline: mcray dba3132 (physics introduced in
14f253d and 448b48b). Comparison adapters are working-tree additions. No
production physics source changed during these checks. These results are
pilot diagnostics, not a completed statistical acceptance campaign.

## Frozen references and common background

CRPropa code: 73fa865b3304113da2059755a83f5939c1039d68.
CRPropa3-data: 34b10e27375c0e9894d7c984ac15fb96e12975ef.
CRPropa was built from a Git archive in Linux amd64 Docker; its checkout was
left unchanged. Main decay comparisons load the approved corrected generator
output in CRPropa and the corresponding full-precision mcray bundle.

The user approved an identical tabulated photon density in the comparison
adapters. Export dn/dE in cm^-3 eV^-1 at z=0 from mcray, with 4097 logarithmic
nodes. Both convolution adapters use the agreed log-log interpolation and
linear density against log(E) on zero-ended intervals. CRPropa transport
uses rates preconvolved with this density. Its runtime TabularPhotonField
interpolates differently, but only supplies static redshift scaling here.
Production temperatures and constants remain unchanged. Their effect on
native-background results has not yet been measured. A common CMB+Inoue12
background was exported, but the TALYS rate diagnostics below use CMB only.

## Deterministic decay checks

The audit covers 158 bound isotopes, 92 active parents, 66 inactive/stable
states and 103 channels. The combined reachable network has 160 states with
free p/n. The normalization changes affect Be11, B12, N16, Na20 and Al24;
Be11 B-P is included as previously approved.

At gamma=10^6, 10^8, 10^10 and 10^12, all 632 mass/rate rows agree with the
corrected CRPropa reference. Maximum relative decay-rate difference is
3.637e-7, consistent with the official lifetime writer's serialization;
maximum relative bound-mass difference is 2.22e-16. All 264 inactive rows have
zero rate in both codes. The three external GoogleTest assertions also pass
using the Linux amd64 mcray probe.

## Isolated decay transport pilot

Thirteen scenarios cover B12, Na20 and Be11 at optical depths 0.1, 1 and 3,
stable Fe56, and prompt Be8, He5 and Li5. Each code uses 20 independently
seeded batches of 1000 primaries per scenario. All 520 jobs completed in
342 seconds in Docker Linux amd64 emulation. Every history passed baryon
number and charge bookkeeping, including charged leptons. Free-neutron decay
is disabled in these isolated fixtures.

Compare batch means and combine the two standard errors in quadrature; do
not treat secondary particle rows as independent samples. Maximum absolute
diagnostic z score for primary survival is 1.35. The largest score across
reported integrated yield/energy observables is 2.73. These scores do not
establish equivalence: the simultaneous observable family, precision budgets
and final validation seeds are not frozen. Rare Be11 B-P needs more events.
Energy is reported under the approved approximate beta kinematics, without
requiring exact event energy balance.

## Conditional beta spectra and a shared numerical limitation

Conditioned B12 beta-minus and Na20 beta-plus samples contain 20,000 events
per code and isotope. Separately prepared benchmark bundles select one
channel; production bundles are unchanged. All three GoogleTest assertions
per isotope pass: complete events, compatible code means within the pilot
five-standard-error diagnostic, and consistency with the discrete sampler.

An independent continuum integration identifies a shared limitation of the
51-node CRPropa beta CDF adopted by mcray. Relative to that continuum model,
its predicted mean charged-lepton energy is lower by about 1.90%, and the
mean neutrino energy higher by about 1.91%. For B12, Q=13.368895514 MeV:
continuum/discrete charged means are 6.96305/6.83043 MeV after dividing lab
energy by gamma; neutrino means are 6.91685/7.04947 MeV. Na20 has the same
pattern, with Q=12.867346281 MeV. Agreement of the two programs alone would
miss this bias. No sampler correction was made; changing the approved
baseline requires a separate decision and validation.

## TALYS rate diagnostic: numerical-reference gate

For Fe56, compare all exclusive channels and the full rate at 61 Lorentz
factors, gamma=10^7..10^13, on the shared CMB table. The official generator
uses 513 rest-energy nodes and 201 gamma nodes. At gamma=3.9810717e9, a
native gamma-table node, mcray gives 0.366415/Mpc versus 0.371472/Mpc from
CRPropa: a 1.36% difference. Near a steep threshold, linear interpolation of
the tabulated rate in log(gamma) amplifies the relative discrepancy; it
reaches 34.4% at gamma=1.2589254e9, where the rate is only about 10^-6/Mpc.
These relative differences must be interpreted together with absolute rates.

A separate SciPy convolution integrates the agreed cross-section interpolant
and the identical photon density, using a 65537-node photon integral and
8193 rest-energy nodes. The CRPropa calc_rate_eps routine was also repeated
on finer rest-energy grids sampling that same interpolant. For the 30
points where the independent total rate exceeds 0.1/Mpc:

| Calculation | Maximum relative difference from independent convolution |
| --- | ---: |
| mcray direct convolution | 0.0364% |
| CRPropa original 513-node quadrature, evaluated at each gamma | 1.3434% |
| CRPropa 2x rest-grid refinement | 0.4360% |
| CRPropa 4x rest-grid refinement | 0.1206% |
| CRPropa 8x rest-grid refinement | 0.0311% |
| CRPropa 16x rest-grid refinement | 0.0078% |

This supports a significant quadrature/interpolation contribution rather than
an unexplained Monte Carlo effect. Refinement samples the agreed mcray
cross-section interpolant, so it is a diagnostic reference, not the unchanged
CRPropa data pipeline. Full all-isotope/exclusive convergence and independent
photon-integral convergence remain to be checked. Do not declare the proposed
0.1% rate budget satisfied for the complete network from this Fe56 scan.

Both combined TALYS+decay adapters additionally completed a 20-history smoke
run with baryon/charge checks. This verifies the adapters can execute; it
establishes no distributional agreement or step convergence.

## Reproducibility and remaining work

Artifacts are outside Git under ../validation:

- nuclear-decay-34b10e2-source: decay rate tables, beta samples and logs,
  photo-rates-fe56-shared-comparison.json, photo-quadrature-diagnostic-final.json.
- decay-data-audit-20261008-complete: lifetime and branching audit.
- decay-crpropa-pilot-20261008: 520 outputs, jobs.json, completion checksums
  and comparison-diagnostic.json.
- crpropa-shared-cmb-20261008-final: hashed common photons, original CRPropa
  convolution tables and pinned generator sources.
- crpropa-decay-reference-20261008-final: source manifest, base/reference
  image IDs and build log for the latest adapters.
- talys-combined-adapter-smoke-20261008: tiny combined control.

Thirteen targeted CTest entries pass in Linux amd64, including the new shared
photon-table tests and the no-interaction regression. Native targeted checks
also pass. This turn did not rerun the entire earlier 96-test suite.

Long ensembles are provided by tests/crpropa/run_external_linux.sh, with
resume checks and a local 1800-second wall limit. The extended TALYS matrix
has 480 jobs (12 cases, 20 batches, two codes, 10000 primaries per job);
only its plan was generated in talys-crpropa-external-plan-final-20261008.
The final 52-job decay smoke matrix completed in 33.7 seconds, and its
resume path verified all completion checksums without rerunning jobs. Do not launch it before choosing how to handle
the numerical-reference discrepancy and freezing acceptance observables.
The proposal to retain the official CRPropa reference and additionally use
an independently converged reference is pending user approval.

Still required: original-data transport control, rare-channel precision,
step convergence for competing processes, all-isotope/exclusive TALYS rate
convergence, native-background effect, CMB+EBL and cosmological controls,
and the frozen simultaneous statistical acceptance campaign. The external
runner currently reports integrated yields and energies, not binned spectra
or a simultaneous equivalence decision.
