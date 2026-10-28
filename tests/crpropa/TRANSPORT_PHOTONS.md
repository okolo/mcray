# Server preparation for transport and secondary-photon acceptance

The implemented entry point runs independent **pre-acceptance pilots and
numerical refinements**, retaining complete primary histories. It cannot
issue an acceptance verdict or start final acceptance ensembles. Final
observable families, fixed spectral intervals, normalizations, fresh seeds,
fixed sample sizes and numerical envelopes must be frozen from these pilots
and approved before the final acceptance campaign. Near-zero absolute budgets
remain pending; unresolved rare outputs must remain inconclusive.

The existing approved budgets are simultaneous 99% family confidence,
1% for resolved integrated observables, 5% for fixed spectral intervals and
at most 200000 primary histories per program/scenario. The family correction
must include every final scenario and observable. Pilot bins and seeds are
not final acceptance bins or seeds. A zero/constant pilot sample is not a
zero-width acceptance interval.

## Implemented scope

- Fe56, shared static CMB and CMB+Inoue12, fixed physical path, no redshift
  energy losses, EM interactions, thinning, free-neutron decay, GZK or PPP.
- TALYS with PD photons Off/On; TALYS+decays with both sources Off, PD only,
  decay only and both On. Photons propagate freely to the observation surface.
  This checks intrinsic secondary production, not a propagated EM cascade.
- CRPropa3 pinned at 73fa865b3304113da2059755a83f5939c1039d68. The native
  reference adapter only toggles public module options; upstream code is
  unchanged. Corrected decay data and masses come from the frozen bundle.
- The official CRPropa convolution generates PD line/transition ratios at
  its 201 gamma nodes on the identical shared photon density. Division and
  `nan_to_num`, nearest-node selection and Bernoulli saturation follow the
  pinned generator/runtime. All eight approved missing-photon transitions
  remain omitted as in the current baseline. mcray retains direct convolution
  with `--pd-photon-yields CRPropa`. These different numerical representations
  and the known finite-upper-limit omission must be reported separately;
  this workflow does not correct upstream rates or hide them in a tolerance.
- Photons use the agreed ultrarelativistic boosts and non-debited energies.
  Exact total event-energy conservation is not an acceptance assertion for
  this baseline. Energy sums, photon energy and beta products are retained
  to measure the approximation; T5/T6 changes are outside this workflow.

The `initial` matrix has E/A=10^18.6 eV and L=1 Mpc on both backgrounds,
12 source/process cases. The `full` matrix covers E/A=10^18, 10^18.6, 10^19 eV
and L=1,10 Mpc, 72 cases. The pilot profile has 20 batches of 1000 primaries
per program/case/numerical setting. Three settings are included:

| Program | Level 0 | Level 1 | Level 2 |
| --- | --- | --- | --- |
| mcray background log10 step / tolerance | .01 / 1e-4 | .005 / 1e-5 | .0025 / 1e-6 |
| CRPropa maximum path-step fraction | .03 | .01 | .003333333333333333 |

mcray's isolated stochastic engine samples competing interactions directly;
its levels refine the background convolution. CRPropa's levels refine the
ordered propagation modules. The same seeds pair refinement batches inside
one program/case; programs/cases have independent seeds. Estimate refinement
noise from paired batch differences, not from two noisy means alone.
The initial pilot contains 1440 jobs and 1440000 primary histories including
all programs and refinements. The full pilot has 8640 jobs and 8640000
histories. They can be expensive; start with the smoke profile, then initial
pilots, and inspect results before choosing further scenarios. No runtime
estimate for the complete pilot is claimed.

## Launch on native Linux amd64

Reuse the original `rate-inputs` bundle: `nuclei/`, `cmb/`, `cmb-inoue/` and
`inputs-manifest.json`. Its hashes are verified before use. New tools come
from the requested fresh mcray checkout, not the older adapter files in that
bundle. The private workspace is never needed on the server.

Use the existing system build dependencies (Git, CMake >=3.24, make,
GCC/G++/gfortran, GSL and Boost filesystem/system), Python >=3.10 with
NumPy/SciPy/GitPython, and GitHub access for mcray, GoogleTest and pinned
CRPropa sources. Docker is not required. For the previously used layout:

```bash
source /home2/ok/mcray-validation-venv/bin/activate
python -m pip install 'cmake>=3.24,<4' numpy scipy GitPython
hash -r
cmake --version

# One-time launcher checkout; the runner creates a fresh calculation checkout.
git clone --branch talys https://github.com/okolo/mcray.git /home2/ok/mcray-launcher
```

Start inside tmux, supplying the published commit containing this workflow as
`REMOTE_REF`. An immutable commit is preferred to the moving branch `talys`.

```bash
bash /home2/ok/mcray-launcher/tests/crpropa/run_external_transport_photons.sh \
  /home2/ok/rate-inputs /home2/ok/mcray-test-data REMOTE_REF smoke initial

# After inspecting the smoke archive:
bash /home2/ok/mcray-launcher/tests/crpropa/run_external_transport_photons.sh \
  /home2/ok/rate-inputs /home2/ok/mcray-test-data REMOTE_REF pilot initial
```

The smoke profile uses two batches of four histories, one numerical level,
48 jobs/192 histories. It checks wiring and bookkeeping only. Every launch
fetches/freezes the requested source commit, builds Release, runs CTest,
builds the pinned native CRPropa adapter, prepares photon-reference copies,
runs jobs and compresses success/failure results. Long calculations are never
started locally. `--plan-only` may be used locally to inspect the matrix:

```bash
python3 tests/crpropa/run_transport_photon_pilots.py \
  --bundle /path/rate-inputs --output /path/new-plan --plan-only
```

The wrapper and `tests/server/run_tests.py` can alternatively be transferred
with their relative directory layout preserved. Only these two launch files
are needed outside the fetched public checkout.

## Results and resumption

The runner stores `.tar.gz`, `.json` and `work/<run-id>/` under the selected
test-data directory. The workload writes `results/transport-photons/`:

- `jobs.json`: frozen cases, levels, seeds, histories, source/tool/input hashes.
- Prepared reference manifests, build logs, pinned reference source and build,
  reference/native binary hashes and per-job commands/timing/logs.
- Individual particle TSVs with complete history IDs. Baryon number and charge
  (including charged leptons), completion count and Photon Off controls must
  pass; partial/discarded histories invalidate a job.
- `.done.json` files with per-primary-aggregated observables and output hashes.
- `pilot-report.json`: batch means/errors, separate photon counts/energies,
  photon-producing-history fractions, exploratory spectra, isotope yields,
  survival and paired numerical refinements, plus program differences.
  It explicitly has no acceptance verdict.

A failed outer run is archived as failed. Keep its archive. For an interrupted
ensemble, the inner driver supports `--resume` in the same retained source,
output and environment; source/tool/input identities and both binaries must
match. Completed outputs are rehashed, partial jobs restart with their frozen
seeds. Incomplete reference preparation/build requires a fresh run. The
outer wrapper starts a new run, not an automatic resume. Resumed inner outputs
need a new preserved archive before transfer; the old failed archive does not
change automatically. Avoid SIGKILL; machine failure can leave only the work
folder. Reference source/build trees are retained and can make archives large.

After the pilots, inspect rare-bin occupancy, reference systematic differences,
step/background refinement noise and photon-source controls. Freeze the final
protocol and obtain approval before preparing/running fresh final ensembles.
The statistical acceptance executor is not yet implemented: it depends on
that approved concrete protocol. A completed pilot is not physical acceptance.

## Short verification

```bash
python3 tests/crpropa/test_transport_photon_pilots.py
python3 tests/crpropa/test_transport_photon_tables.py
```
