# Deterministic photodisintegration controls

The PD-only tool below remains available. The implemented
[PD+decay energy-state extension](COUPLED.md) supplies the joint reference;
its separate contract and refinement controls apply when decays are enabled.

These validation tools propagate independent first moments through the complete
PD-only Fe56 network. They do not run Monte Carlo, change production physics,
correct official control outputs, estimate history variances, or issue an
acceptance verdict. Nuclear decays, losses, cosmology, EM propagation and energy
debits for photons are outside this implemented scope.

For the approved PD kinematics, every nuclear product retains the injected E/A.
Thus one energy per isotope is exact in this scope. This reduction **cannot** be
used for combined PD+decay transport: decays preserve gamma instead of E/A.

## Operator and coefficient sources

For mean isotope populations n, solve dn/dx = Q n. Each parent has diagonal
-R_full; a channel contributes R_full P_channel times each product multiplicity,
where P_channel = R_channel / sum(R_channels). Full totals and exclusive sums
are intentionally distinct. Photon rewards integrate the channel-conditioned
Bernoulli line probabilities along the path. Line energies follow a uniform
distribution on [0, 2 gamma_residual E_rest]; bin counts and energies are
integrated analytically. Saturation and the eight explicit no-photon exceptions
are preserved. First moments do not require cross-line independence, but do not
determine variances or the probability of at least one photon per history.

Three sources fill the same operator:

- I: independent log-interval Gauss-Legendre convolution of the frozen raw
  cross sections and shared photon background. The photon integral is analytic
  on each background interpolation interval. No production integrator is used.
- M: full/exclusive rates and conditional line yields exported by the actual
  production mcray functions through `deterministic_pd_export`.
- C: unchanged pinned official CRPropa tables and serialized masses, linear
  total-rate interpolation in log(gamma), nearest-parent-node branches and
  nearest-residual-node photon probabilities. The remaining serialized branch
  probability goes to the last branch exactly as in the pinned runtime; rounded
  tables are not silently renormalized. Outside the open official gamma support,
  the official interaction rate is zero. Unsupported photon indices fail.

Compare MC_mcray to D_M, D_M to D_I, MC_CRPropa to D_C, and D_C to D_I
separately. An official coefficient discrepancy is not evidence of a defect in
mcray sampling; agreement with one's own coefficients does not establish their
physical accuracy.

`transport.py` solves the augmented sparse operator with `expm_multiply` and
optionally independently integrates it with DOP853. Energy reward coordinates
are scaled before solving. The coefficient schema is
`mcray-deterministic-transport-v1`: `initial`, `energy_scale_eV`, `states`, and
optional `photon_bin_edges_eV` partitioning [0,infinity), with final JSON null.
Each state has `id`, integer `a,z`, `energy_eV`, `rate_per_Mpc`, and `channels`.
Each channel has `probability`, product-id/integer-multiplicity mapping
`products`, and optional `photons` with `mean_energy_eV,probability`.
Malformed rates, probabilities, coverage, identities and nonconserving A/Z
channels fail. Small arithmetic guards in the solver are not acceptance budgets.

## Run one case

Run from the public repository root. Use the NumPy/SciPy environment used for
the frozen reference preparation; Python, NumPy and SciPy versions are recorded.
Build the new exporter through the normal Release build, or explicitly:

```bash
cmake --build build/Release --target deterministic_pd_export -j2
python3 -m unittest discover -s tests/deterministic -v
python3 tests/deterministic/run_pd_reference.py \
  --data /path/to/rate-inputs/nuclei \
  --reference /path/to/retained-run/transport-photons/cmb \
  --provider I --order 64 --energy-per-nucleon 3.9810717055349857e18 \
  --distance-mpc 1 --photons On --path-control \
  --photon-bin-edges 0,1e14,1e15,1e16,inf --output /path/to/new-result
```

For M, select `--provider M --native-probe build/Release/tests/deterministic_pd_export`
and optionally `--background-step .0025 --eps-rel 1e-6`. For C, select
`--provider C`. Repeat for `cmb-inoue`. Photon Off is also supported. References
must contain the prepared official photon tables from the transport-photon
pilots; the earlier disabled placeholder is rejected. Prepared inputs and
reference payloads are verified against their manifests. Output directories
must be new. Each run records the source revision/status, tool/input/binary
hashes, dependency versions, settings, elapsed time and output hashes, including
failure/incomplete status. The native coefficient export is bounded to 600 s.
No seed is required and no primary histories are consumed.

Compare an existing case to matching complete retained pilot batches:

```bash
python3 tests/deterministic/compare_pd_pilots.py \
  --run /path/to/new-result --pilots /path/to/retained-run/transport-photons \
  --background cmb --code mcray --level 2 --output /path/to/new-comparison.json
```

Use M for the native comparison and C for `--code crpropa`. The adapter checks
case settings and payload hashes, and reports batch-standard-error residuals.
Constant/zero samples remain unresolved. These exploratory residuals are not
simultaneous 99% intervals. The sample bin edges above are diagnostic choices,
not frozen acceptance bins. Final families, bins, numerical envelopes and seeds
remain to be fixed; existing 1%/5% budgets and 200000-history cap are unchanged.

## External Linux refinement

Checks expected to exceed 30 minutes belong on an external native Linux amd64
server. A prepared wrapper reuses the existing server clone/build/CTest/archive
runner and runs I at orders 16/32/64, M at background steps .01/.005/.0025, and
C, on both backgrounds at the initial pilot point, with independent path checks.
It generates zero Monte Carlo histories. Use a Python environment with NumPy and
SciPy on PATH, and retain the prepared reference directories from the pilot run.
The requested revision must first be available from the configured remote:

```bash
bash /path/to/launcher/tests/deterministic/run_external_pd_reference.sh \
  /home2/ok/rate-inputs/nuclei \
  /home2/ok/mcray-test-data/work/PILOT_RUN/results/transport-photons \
  /home2/ok/mcray-test-data REMOTE_COMMIT
```

The runner resolves/fetches that revision, records the exact commit, builds and
runs CTest, then executes the workload from its fresh public checkout. The
launcher and input paths must exist on the server; no private repository is
needed. Results and failures are compressed in TEST_DATA_DIR. Each coefficient
run additionally verifies frozen payloads and hashes its tools and outputs.

See [confirmed initial results](../../docs/validation/deterministic_pd_reference.md)
for the tested scope and remaining limitations.
