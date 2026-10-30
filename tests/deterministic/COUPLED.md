# PD+decay deterministic energy-state reference

This extends the [PD-only control](README.md) to the complete prepared Fe56
network with nuclear decays. It retains the approved physical model and supplies
first moments, not history variances, coincidence probabilities or acceptance
verdicts. Free-neutron decay is disabled, matching the isolated retained pilot
fixtures. Cosmology, losses, GZK, PPP, thinning and EM propagation are disabled.

## Energy and rate contract

PD products retain their parent's E/A. Beta transitions retain the incoming
gamma: the nuclear E/A changes by the daughter/parent mass ratio. Compound
decays execute beta-minus, beta-plus, alpha, proton and neutron emissions in
that order. Emission after the beta steps shares the resulting E/A; identical
products retain their multiplicities. Pure spectator breakup follows the
prepared mcray contract. The currently supported pinned network is checked;
missing daughters, non-DAG transitions or incomplete decay/stability records
are errors.

For each isotope, all permitted PD and decay paths give lower and upper E/A
bounds, propagated in descending A and mass order. The bounds include channels
even when their rate at the current energy is zero. Each non-degenerate interval
receives N equally spaced energy nodes; a point interval needs only one node.
Every nuclear product is deposited at its two surrounding daughter nodes with
positive linear-energy weights. These weights preserve number, charge, baryons
and mean nuclear energy. Targets outside the proved interval fail; a 32-ulp
arithmetic guard only permits floating-point endpoint roundoff. There is no
physical overflow clipping or renormalization of escaped particles. This
discretizes nonlinear subsequent rates/rewards, so energy-grid refinement is
required even though first-moment energy is conserved.

PD coefficients are evaluated directly at each isotope/node energy, using the
same I/M/C adapters as the PD-only tool. There is no extra rate-interpolation
grid. In particular, C retains its nearest-node branch/photon discontinuities.
The M exporter accepts a complete isotope/E-per-A grid through `--grid` and
exports production total decay rates, partial channel rates and photon lines
with `--combined On --decay-photons On|Off`. Decay kernels use tabulated masses;
M and I retain mcray's 0.511 MeV electron mass and Mpc/second convention, while
C retains the pinned official constants, masses and corrected decay table.
Conventions are recorded explicitly, not adjusted to force agreement.

All decay lifetimes remain finite. The generator contains both PD and decay
rates, including their competition. Very short lifetimes make the system
stiff; sparse implicit BDF is the default path solver, with Radau and tighter
solver tolerances available as an independent numerical check. No prompt state
is replaced with an instantaneous decay. The input solver tolerances are
numerical controls and do not replace or relax physical acceptance tolerances.

## Secondary rewards

- PD lines use residual gamma and saturated Bernoulli probabilities, including
  the eight approved no-photon exceptions.
- Decay lines are attached once to the selected channel before its nuclear
  transitions and use the incoming parent gamma.
- Each line has an analytic uniform lab-energy distribution on [0,2 gamma E].
  Total, per-source and user-specified photon-bin number/energy means accumulate
  in reward coordinates.
- Beta leptons retain the existing 51-node allowed inverse-CDF approximation.
  Probability is uniform within each nonzero CDF interval, giving its exact
  interval midpoint mean energy. Isotropic boosts give exact first moments of
  this approximation. Compound beta steps use their own mass endpoints.
  Electron, positron and electron-neutrino/antineutrino count/energy rewards
  retain charge and the nuclear-plus-lepton mean energy. Lepton spectra and
  correlations are not supplied by this version.

Nuclear energy is not debited for either photon source. The reported energy
invariant deliberately excludes this additional photon energy, matching the
accepted compatibility model. Small negative populations from solver roundoff
are retained diagnostically rather than clipped into apparent rare yields.

## Run and refine

Use the existing NumPy/SciPy reference environment and normal Release build
from the public repository root; rebuild `deterministic_pd_export` after this
extension. All paths below refer to prepared data and reference files.

```bash
cmake --build build/Release --target deterministic_pd_export -j2
python3 -m unittest discover -s tests/deterministic -v
python3 tests/deterministic/run_coupled_reference.py \
  --data /path/to/rate-inputs/nuclei \
  --reference /path/to/pilot/transport-photons/cmb \
  --provider I --energy-per-nucleon 3.9810717055349857e18 \
  --nodes 5 --order 32 --distance-mpc 1 \
  --pd-photons On --decay-photons On --path-control \
  --output /path/to/new-coupled-result
```

Select M with `--provider M --native-probe build/Release/tests/deterministic_pd_export`
or C with `--provider C`; repeat on `cmb-inoue`. Each photon source can be
disabled separately. `--photon-bin-edges 0,1e14,1e15,1e16,inf` is the default
diagnostic partition, not a frozen acceptance spectrum. Runs save manifests,
source/data/binary hashes, numerical parameters, every coefficient slice,
the isotope grids, decay kernels, sparse operator and resolved energy-state
populations. Failures retain incomplete status. The native exporter has a
600-second limit per coefficient slice.

Compare new directories at different node counts while holding other controls
fixed:

```bash
python3 tests/deterministic/compare_coupled_refinement.py \
  --coarse /path/to/nodes5 --fine /path/to/nodes9 \
  --output /path/to/new-grid-comparison.json
```

The comparison verifies output hashes, input identity and physical settings;
it retains every isotope and reward difference. Results at/below solver
coordinate tolerance have no relative precision claim. An observed difference
is not a certified bound on all uncomputed refinements. Separate comparisons
must refine quadrature, the native background integration and path tolerances.
`compare_pd_pilots.py` also accepts these runs and selects the matching
`talys-decay` pilot cases with the same photon flags. No new histories are
generated; pilot batch residuals remain exploratory.

## External Linux campaign

Local runs are restricted to at most 9 nodes and order 32. Larger grids require
`--external-server` on Linux amd64. Use the prepared wrapper for campaigns
expected to exceed 30 minutes:

```bash
bash /path/to/launcher/tests/deterministic/run_external_coupled_reference.sh \
  /home2/ok/rate-inputs/nuclei \
  /home2/ok/mcray-test-data/work/PILOT_RUN/results/transport-photons \
  /home2/ok/mcray-test-data REMOTE_COMMIT full
```

Activate the NumPy/SciPy environment first. The requested commit must be
available on the public remote before running this command. The runner fetches
and resolves the exact revision, builds, runs CTest and Python controls, and
archives all results or failures in TEST_DATA_DIR. The workload uses a fresh
public checkout.

`initial` runs N=3,5,9 with I order16; `full` runs N=17,33,65 with I order32.
Both profiles cover I/M/C on both backgrounds. Additional runs double I order
and halve M background step at the finest energy grid, retaining separate
observable-level refinement reports. Every case also runs the other stiff
solver with tenfold tighter rtol/atol. There are zero primary histories and no
seeds. Neither profile freezes the final statistical protocol or issues
acceptance; the 99% simultaneous confidence, 1%/5% budgets and 200000-history
limit remain unchanged.

See [initial implementation results](../../docs/validation/deterministic_coupled_reference.md).

## Audit a completed full matrix

`analyze_coupled_matrix.py` reads an extracted server result directory with
`manifest.json` and `deterministic-coupled/`. It verifies payload, case, input,
tool and output identities, recomputes the moment/bin invariants and all
archived refinement reports, separates tighter-Radau grid differences from
path errors, and checks restricted energy-state populations. It uses Python
3.10+ standard-library modules and Git; the frozen source revision must be
available in the public checkout's local object database.

```bash
python3 tests/deterministic/analyze_coupled_matrix.py \
  --results /path/to/extracted/RUN_ID --output /path/to/new-analysis
```

This runs no transport and generates no histories. The output records the
auditor's SHA-256 and deliberately supplies no acceptance verdict. Nodal
restriction differences are not final physical spectral-bin errors. See
the [audited full-matrix results](../../docs/validation/deterministic_coupled_full_matrix.md)
for the supported conclusions and rare-output limitations.
