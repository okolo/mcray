# Nuclear decays with tabulated TALYS masses

## Approved scope (2026-10-08)

NuclearDecay is a separate stochastic process. CRbeam --nuclear-decay On|Off
selects it explicitly; the default is Off. It requires a nuclear primary and
--talys-data containing masses and decay records. Photodisintegration can be
TALYS or Off. Combining nuclear decays with PSB is rejected: PSB keeps its
original physical chain and approximate masses. Free neutrons continue to
use the existing NeutronDecay process in CRbeam.

Include beta-minus, beta-plus, alpha, proton/neutron emission and compound
channels. Create charged leptons and electron neutrinos/antineutrinos;
decay photons are disabled. Primary production runs remain limited to Fe56.
The test scenario runner also accepts daughter isotopes to isolate decays.

## Rates and event handling

The table stores a partial mean lifetime tau_i in seconds for each exclusive
channel, including its branching fraction and the CRPropa ionization
correction. The laboratory rate is sum_i(1/tau_i)/gamma, with gamma=E/M.
mcray advances physical time in natural units; units.second supplies the
conversion. Do not apply CRPropa's additional 1/(1+z) factor: that factor
converts its physical light-travel rate to comoving distance.

The propagation engine samples an event using the total rate; NuclearDecay
selects one channel with probability (1/tau_i)/sum_j(1/tau_j). The parent is
terminated. Each product retains the parent weight and inherited particle
metadata. Clones share immutable nuclear masses and decay data.

Compound channels follow the CRPropa sequence: beta-minus, beta-plus, alpha,
proton, neutron. Every intermediate identity needs an available tabulated
mass. NuclearDecayData validates the combined decay/photodisintegration
network before output creation. A missing decay entry is an error unless an
explicit stable-state record exists. Neither the photodisintegration terminal
whitelist nor absence from nuclear_decay.txt establishes physical stability.
Artificial 1e-99 s proton/neutron dripping from the original data generator
is rejected; it is not needed in the audited Fe56 network.

## Approved kinematics and limits

For beta decay, retain the residual nucleus's Lorentz factor and neglect
recoil. Sample the allowed electron/positron spectrum proportional to
E sqrt(E^2-me^2) (Q+me-E)^2, with Q=Mparent-Mdaughter-me. Use the same
51-node cumulative discretization and inverse linear interpolation as CRPropa.
The implementation scales energies by me while forming the CDF to avoid
unnecessary large/small intermediate powers. Sample an isotropic rest-frame
angle and use CRPropa's ultra-relativistic boost for the charged lepton and
neutrino. Charge conservation includes the emitted lepton.

This approximation does not conserve energy exactly in each beta event and
has no Coulomb correction. It is intended for relativistic cosmic-ray nuclei;
near-rest-frame propagation needs a separately approved treatment. Positive
beta endpoints and finite rates/energies are required.

For alpha and nucleon emission, assign E_i=Eparent*A_i/Aparent sequentially.
Keep both products even when an emitted alpha and residual are identical.
These emission steps conserve the assigned energy, baryon number and charge.

## Data preparation and the CRPropa corrections

Pin CRPropa3-data 34b10e27375c0e9894d7c984ac15fb96e12975ef. The preparation
tool copies the exact calc_mass.py, calc_decay.py, units.py, gitHelp.py and
NuDat3/NIST input tables from Git blobs into an external output bundle. It
runs the original generators for a reference and a separately patched copy
for production. It does not modify the reference checkout or the existing
photodisintegration-only bundle. NumPy, SciPy and GitPython are preparation
only dependencies; the C++ runtime does not use them.

The user approved two changes to the pinned decay generator:

1. Recompute the branching-fraction sum after converting an inclusive branch
   to exclusive branches, before normalization. The original generator divides
   by the earlier, inclusive sum, reducing the total decay rate incorrectly.
2. Map Be11 B-P to channel 10010 (beta-minus followed by proton emission).
   The original generator silently skips this known mode. Its corrected
   exclusive fraction is 8.3e-6.

Affected reachable isotopes are Be11, B12, N16, Na20 and Al24. For example,
Na20's pre-ionization fractions change from 0.6659725 and 0.1670137 (sum
0.8329863) to 0.7995 and 0.2005. Its corrected partial mean lifetimes after
the unchanged ionization correction are 0.8082378 s and 3.2228733 s, versus
0.9702895 s and 3.8690594 s in the reference. Compare production against the
corrected reference when validating these states; exact agreement with the
unmodified CRPropa table is not expected.

Keep CRPropa's fully ionized beta-plus/electron-capture correction. A channel
with a forbidden beta-plus endpoint becomes inactive; a state with no remaining
active decay is recorded explicitly as stable for this model. This can differ
from the decay of a neutral atom. The current audit finds 160 reachable states
including p/n, 158 bound isotopes, 103 active decay channels and 66 bound
states without an active decay. All required masses and intermediate identities
are already in the prepared photodisintegration network. The 13 approved
photodisintegration terminal states stay unchanged; some now decay.

From the repository root, after preparing the photodisintegration bundle:

```bash
python3 -m venv /tmp/mcray-decay-venv
/tmp/mcray-decay-venv/bin/pip install numpy scipy GitPython
/tmp/mcray-decay-venv/bin/python tools/prepare_nuclear_decay_data.py \
  --repo ../validation/crpropa-data-34b10e2 \
  --talys-data ../validation/talys-fe56-34b10e2 \
  --output ../validation/talys-decay-fe56-34b10e2-final
python3 tools/prepare_nuclear_decay_data.py \
  --verify --output ../validation/talys-decay-fe56-34b10e2-final
```

The output must not exist at preparation time. It contains the original TALYS
files plus decays.tsv, decay_stable.tsv, decay_provenance.txt,
decay_manifest.json and decay_sources/{reference,corrected}. Each source,
generated artifact, preparation script and dependency version is recorded.
SHA-256 verification is performed by the preparation/verification tool;
the C++ loader checks the schema, revision, lifetimes, channels, masses,
endpoints and network structure. Verify both manifests before production use.
The original TALYS provenance still describes its photodisintegration-only
preparation; decay_provenance and resolved_nuclear_decay record the extension
and activation separately.

CRbeam usage, from bin:

```bash
../build/Release/CRbeam --particle nucleus --A 56 --Z 26 \
  --photodisintegration TALYS --nuclear-decay On \
  --talys-data ../../validation/talys-decay-fe56-34b10e2-final \
  --z 0.01 --emax 1e21 --nparticles 1000 --output /tmp/fe56-decay-run
```

Other CRbeam interactions and output filters remain active according to their
existing settings. In particular, --noEMcascade discards propagated charged
leptons/photons; it does not prevent their creation at a decay vertex. Isolated
conservation checks must retain leptons and count their electric charge.

## Verification and remaining validation

GoogleTest covers physical lifetime conversion, inverse-gamma scaling,
absence of an extra redshift factor, stable-state handling, clones, exclusive
branch sampling, beta-minus/plus lepton signs, alpha duplicates, compound
beta-proton sequence, parent identity, weights, preflight failures and
exponential survival through PropagationEngine. The transport survival
control uses 10000 histories and a five-sigma binomial diagnostic.

The external-data scenario runner supports --process decay and talys-decay.
For decay-only daughter controls, --a and --charge specify the isotope; zero
charge selects its PSB-chain representative. Free-neutron decay is not added
by this isolated runner. Use --particle-output and nuclear_decay_comparison
to check complete histories, baryon number, electric charge with leptons and
exponential primary survival. --alpha-only additionally requires exactly two
He4 products for every Be8 history. These checks use GoogleTest and are not
rare-output precision acceptance tests.

Native controls include Be8 (100 histories, 200 He4 products), B12 (20000
histories, finite-path exponential survival), combined TALYS+decays on CMB
(1000 Fe56 histories, E/A=10^18.6 eV, L=1 Mpc), and a two-thread CRbeam
run. All histories complete and baryon number/charge checks pass. Beta energy
balance is assessed under the documented approximation rather than exact
conservation assertions. Outputs and logs are outside Git under
../validation/nuclear-decay-34b10e2-source.

The final short suite has 96 passing CTest entries on macOS arm64 (29.66 s)
and Docker Linux amd64 under emulation (95.43 s). The same Be8, B12 and
combined Fe56 external controls pass their applicable GoogleTest assertions
on both platforms. The two-thread CRbeam control also completes on both.
B12 survival is 14380/20000 versus 14388.41 expected (difference -0.13 sigma);
Be8 yields exactly 200 He4 products. The 1000-history Fe56 control retains
543 Fe56 primaries versus 567.835 expected (difference -1.59 sigma).
These counts diagnose implementation consistency, not full physical acceptance.
This implementation does not establish full CRPropa E2E agreement, precision
of rare branches, accuracy of the ionization correction, or nuclear GZK/pair
support. Dedicated rate, branching and spectral comparisons remain necessary.

The follow-up CRPropa comparison stages and pending acceptance proposals are
recorded in nuclear_decay_comparison_plan.md.
