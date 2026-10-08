# Unified nuclei and photodisintegration

Implementation date: 2026-10-08. Development branch: talys.

## CRbeam input and process selection

Use a named nuclear primary with explicit unsigned A and optional Z:

```sh
CRbeam --particle nucleus --A 56 --Z 26 --photodisintegration PSB --z 0.01
```

Both A and Z default to zero. For nucleus, A>=2 is mandatory and Z=0 selects
the PSB-chain representative. The first stage accepts only primary Fe56;
`--A 56` therefore resolves Z=26. Lowercase --z specifies source redshift.
For other particles A,Z are physically unused. Existing named particles and
numeric --primary remain available; do not supply both selectors.

--photodisintegration accepts exactly TALYS, PSB or Off and defaults to Off.
PSB preserves weighted branching and approximate nuclear masses. TALYS samples
one exclusive channel per event and uses pinned tabulated nuclear masses.
Off disables photodisintegration and uses the approximate PSB source mass.
The mode does not disable redshift losses or existing electromagnetic,
free-nucleon GZK/pair, magnetic-deflection or free-neutron-decay processes.
Nuclear GZK, nuclear pair production, nuclear decays and photon emission are
not implemented in this stage.

Nuclear raw output files are AXXZYY, for example A56Z26 and A04Z02. Different
isotopes at the same A receive different files. Existing non-nuclear filenames
and raw row columns are retained. SpectrumOutput appends separate named
isotope columns after the light-particle columns and writes their column
header when nuclei are present. Trajectory logs append A,Z,mass/MeV. Type
filters use Nucleus to select every bound nucleus; TotalEnergyOutput(Nucleus)
is intentionally the energy summed over all nuclear isotopes.
The params file records resolved source identity, mass model, selected mode,
data provenance, checksums and terminal-state assumptions for TALYS runs.

## Preparing the external TALYS data

The approved revision is CRPropa3-data
34b10e27375c0e9894d7c984ac15fb96e12975ef. Prepare data outside this Git
repository. Example commands from the mcray root, using the existing sibling
validation checkout:

```sh
python3 tools/prepare_talys_data.py \
  --repo ../validation/crpropa-data-34b10e2 \
  --output ../validation/talys-fe56-34b10e2 \
  --grid-policy crpropa-padding
python3 tools/prepare_talys_data.py \
  --output ../validation/talys-fe56-34b10e2 --verify
```

The preparation command requires a checkout at the exact approved revision
and reads committed blobs using git show. It refuses an existing destination.
For a fresh checkout, obtain https://github.com/CRPropa/CRPropa3-data at that
revision. A sparse/partial clone can fetch required committed blobs on demand.
No TALYS executable or background-specific rate tables are required.

The output contains masses.tsv, terminals.tsv, grid_external.tsv,
grid_talys.tsv, total.tsv, channels.tsv, provenance.txt and manifest.json.
Cross sections remain in mb and photon energies in MeV in these files.
No convolution occurs during preparation. Original nodes/values are retained;
logarithmic grids are extended with zeros to 513 nodes. Source and prepared
payload SHA-256 are recorded. Run --verify before a calculation; the C++
loader checks format, revision, dimensions, signs, identities, mass provenance
and network closure, but does not recompute cryptographic hashes itself.

The selected ordinary CRPropa set uses PD_Talys1.8_Khan full totals and thin
exclusive channels for A>=12, and PD_external for A<12. Tabulated bound-nucleus
masses follow pinned NIST then NuDat3 precedence, without approximate fallback.
The reachable Fe56 inventory has 160 states, including free p,n; all 158
bound-nucleus masses are covered by NIST. Free p,n keep mcray's existing masses.

Example from bin, with an output directory that does not already exist:

```sh
../build/Release/CRbeam --particle nucleus --A 56 \
  --photodisintegration TALYS \
  --talys-data ../../validation/talys-fe56-34b10e2 \
  --z 0.0001 --emax 1e21 --nparticles 8 --background 0 \
  --noEMcascade --output ../../validation/example-talys
```

## Particle and process contracts

Bound nuclei have one ParticleType, Nucleus. Private A,Z, numeric nuclear
mass and a mass-model tag are updated together. SetPsbNucleus(A) resolves the
PSB representative and approximate mass. TalysData::SetParticle resolves the
pinned mass for a bound nucleus, or restores the ordinary free p/n type and
mass. SetType clears nuclear state when changing to a light particle.
Use MassNumber(), AtomicNumber(), Mass(), ElectricCharge() and instance Name().
Static Mass(Nucleus) and ElectricCharge(Nucleus) cannot resolve an isotope
and raise an error. Nuclear identity must be initialized before propagation.

Copies preserve identity, selected mass, ancestry, IDs, arrays and non-owning
custom pointers. Whole-particle memset/memcpy is no longer used. Inconsistent
public Type changes are diagnosed at identity/kinematics/propagation boundaries.
PSB rejects incompatible representatives/masses; TALYS rejects the wrong mass
model or a mass inconsistent with its data set.

TalysPhotoDisintegration is separate from PhotoDisintegration (PSB). Clones
share immutable TalysData and clone their BackgroundIntegral. Full total
cross-section convolution defines the event rate. Every exclusive thin channel
is convolved separately; probabilities are R_c/sum(R_c), not R_c/R_total.
ChannelRates exposes both full total and exclusive sum for diagnostics.
A positive total with zero exclusive sum is an error.

One event removes the parent and produces the residual and all emitted
n,p,d,t,He3,He4 particles from the selected channel. Each product keeps the
parent weight and receives E_parent*A_product/A_parent. A,Z and energy are
conserved under this approved kinematic approximation. Recoil/binding-energy
corrections and nuclear-decay/photon-emission products are absent.

LinearFunc with logarithmic axes and its existing zero-handling parameter
implements interpolation. Positive intervals use log(sigma) versus log(K);
intervals containing zero use sigma versus log(K). No epsilon replacement or
extrapolation is applied. Cross sections vanish outside the padded grid.
A future task assesses the high-energy tail: https://github.com/okolo/mcray_dev/blob/main/docs/tasks/T1/README.md.

## Explicit terminal states and errors

The user approved zero photodisintegration rates for the following missing-data
states, listed as (A,Z): (5,2),(5,3),(6,2),(8,4),(9,5),(12,5),(16,9),(20,11),
(39,20),(42,21),(46,23),(47,24),(50,25). Loading prints warnings and records
this whitelist. These are terminal states of the present PD-only network,
not assertions of physical stability. CRPropa handles nuclear decays separately.

Other missing required masses/cross sections, invalid channels and invalid
rates are errors. CRbeam checks network closure before creating output.
Propagation failures are counted across worker clones; a nuclear run returns
an error if a cascade was discarded. Partial output from a failed run is
invalid. Failed-cascade secondary queues are cleared before another history.

## Current verification and remaining work

GoogleTest covers identity/mass copying, invalid updates, clearing nuclear
state, process mass-model checks, isotope output/logging, synthetic data
validation, full versus exclusive rates, channel frequencies, conservation,
terminal behavior, clone rates and failed-cascade isolation. Real CRbeam tests
exercise Off and PSB; real-data TALYS controls run separately outside Git.

A regression exposed a pre-existing MonochromaticBackgroundIntegral::GetRateK
bug: its upper cross-section bound used Xmin instead of Xmax, forcing the
rate to zero. This is corrected and checked against an analytic constant-sigma
integral and the sampling path. This fix does not affect continuous-CMB/EBL
PSB regression controls.

Before/after migration, a fixed-seed 200-primary PSB CMB+Inoue12 control at
E/A=1e18 eV and L=1 Mpc produced 208 identical rows after normalizing isotope
labels. This demonstrates preservation of that control, not full PSB validation.

Real-data isolated TALYS controls at L=1 Mpc include 1000 histories on CMB
at E/A=10^18.6 eV and 200 histories on CMB+Inoue12 at E/A=1e18 eV. Every
history conserves A=56, Z=26 and energy to rounding precision. Fe56 survival
is consistent with exp(-R*L): 540 versus 567.835 expected (1.78 binomial SE)
and 189 versus 191.872 (1.03 SE), respectively. These are internal controls,
not CRPropa E2E agreement or a rare-output precision claim.

An independent composite-Simpson integral of the full Fe56 cross section
against an analytic 2.7255 K Planck background differs from the current
ContinuousBackgroundIntegral scan by at most 0.1261% over 42 points with
R>=1e-10/Mpc. Its own 128/256-subinterval change is below 2.6e-9 relative.
The diagnostic uses mcray's existing unit constants deliberately. It verifies
convolution independently of the mcray integrator; it is not a comparison to
CRPropa's rate-generation algorithm or a proof that the omitted physical tail
is negligible. Scripts and commands are in tests/README.md.

The Fe56 scan also records the exclusive sum and its ratio to the full rate.
Over those 42 points the ratio ranges from 0.8468 to 1.2195; restricting to
R>=0.1/Mpc gives 0.8511 to 1.0162. Original-node thin/full sums range from
0.7645 to 1.000004. Separate interpolation, particularly the zero-node rule,
can therefore increase the convolved exclusive sum above the full rate near
threshold. The algorithm retains the approved separate total and normalized
exclusive probabilities; it does not force the sum to match the full rate.
Investigate this difference in the direct CRPropa comparison before assigning
numerical compatibility thresholds.

Final short-suite verification: all 88 CTest entries pass on macOS arm64
(28.87 s) and Docker Linux amd64 under emulation (94.90 s). The Linux image
also passes a real-data CRbeam TALYS run with two worker threads and the
1000-history isolated CMB E2E control with all three applicable GoogleTest
checks. Logs, image identity and source-file hashes are retained in
../validation/unified-nucleus-20261008, outside Git. The approved data are in
../validation/talys-fe56-34b10e2, outside Git. TransportCR was not changed.

Detailed ensembles, acceptance criteria, rate/grid convergence and a direct
CRPropa comparison remain to be planned/validated. This implementation is a
first-stage photodisintegration model, not complete nuclear transport.

## Optional nuclear decay extension

The separate --nuclear-decay On|Off process defaults to Off. It requires
an extended data bundle with decay channels and explicit stable-state records;
see nuclear_decay.md for preparation, approved kinematics and corrections to
the CRPropa decay generator. Photodisintegration-only data remain supported
when nuclear decays are Off. The previous statement that decays are absent
applies to the original photodisintegration-only stage.
