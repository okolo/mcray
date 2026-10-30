# CRPropa comparison adapters

These tools compare the unchanged production libraries on Linux amd64.
Read [the nuclear-decay comparison results](../../docs/validation/nuclear_decay_comparison_results.md)
and [the release validation limits](../../docs/validation/nuclear_release_notes.md)
before interpreting results.
Outputs, frozen sources, generated tables and Docker images stay outside Git.
No script declares statistical acceptance from diagnostic z scores.

For the next transport and nuclear-photon gate, see
[the native server workflow](TRANSPORT_PHOTONS.md). It prepares independent
pilots and numerical refinements on both shared backgrounds. Final acceptance
families, bins and fresh seeds remain pending; it cannot issue acceptance.

## Build and prepare

From the mcray repository root, build the existing base image and then the
pinned CRPropa library and adapters:

```bash
bash bin/test-docker.sh -R 'Decay|decay|Shared|nuclear_none'
bash tests/crpropa/build_reference.sh ../CRPropa3 ../validation/reference-build
```

The reference build requires CRPropa HEAD
73fa865b3304113da2059755a83f5939c1039d68. The build archives that exact
revision, freezes adapter sources and hashes, and records both image IDs.
It does not modify the reference checkout. The resulting image is
mcray-crpropa-validation:local. Docker and a working base image are required.

Prepare the approved bundle using the commands in ../../docs/usage/nuclear_decay.md.
Examples below assume ../validation/talys-decay-fe56-34b10e2-final exists.
The Python runner verifies its manifests before running. Audit the original
and corrected generator data with audit_decay_data.py --help.

For the shared TALYS background, use the native test exporter and a Python
environment with NumPy, SciPy and GitPython (record installed versions):

```bash
build/Release/tests/photon_background_export cmb ../validation/common-cmb.tsv
python3 tests/crpropa/prepare_shared_photodisintegration.py \
  --repo ../validation/crpropa-data-34b10e2 \
  --data ../validation/talys-decay-fe56-34b10e2-final \
  --photons ../validation/common-cmb.tsv \
  --output ../validation/shared-crpropa
```

The generator extracts exact pinned source files through git show; it
computes 201-node CRPropa rate/branch curves from the approved raw 513-node
cross sections. It exports a static CRPropa photon field, corrected decay
data, hashes and dependency versions. No production background changes.
Use cmb-inoue instead of cmb to export the agreed Inoue12 combination.

nuclear_decay_probe produces the 632-row decay inventory, conditioned
single-channel events, or Fe56 full/exclusive photo-rates. The latter requires
--mode photo-rates --data DATA --photon-table TABLE --output FILE.
The external GoogleTest tools decay_reference_comparison and
beta_spectrum_comparison consume the resulting files; they are not CTest
entries because the reference data are external. Invoking them without
arguments is not supported. Use decay_reference_comparison MCRAY_TSV CRPROPA_TSV
and beta_spectrum_comparison MCRAY_TSV CRPROPA_TSV Q_MEV GAMMA EVENTS.
prepare_conditional_decay.py creates
explicitly labeled one-channel benchmark bundles; never use those as a
production data source.

To separate rest-energy quadrature from gamma interpolation:

```bash
python3 tests/crpropa/diagnose_photo_rates.py \
  --data ../validation/talys-decay-fe56-34b10e2-final \
  --reference-data ../validation/shared-crpropa \
  --mcray-rates ../validation/photo-rates.tsv \
  --output ../validation/quadrature-diagnostic.json
```

## Short pilots

```bash
python3 tests/crpropa/run_decay_comparisons.py \
  --data ../validation/talys-decay-fe56-34b10e2-final \
  --output ../validation/decay-pilot --profile pilot
```

smoke uses two batches of 1000 primaries; pilot uses twenty. A local
1800-second wall limit stops the active container and leaves completed jobs
resumable. All cases exclude cosmological losses and free-neutron decay.
Decay-only tests need no photon field. CRPropa uses native SimplePropagation
and NuclearDecay with electrons/neutrinos on and photons off.

For --family talys, also supply --reference-data. The matrix covers Fe56
at E/A=10^18, 10^18.6 and 10^19 eV, distances 1 and 10 Mpc, both TALYS
alone and TALYS+decay. CRPropa's ordered modules require a separate step
refinement study. Do not interpret the default step fraction 0.03 as
converged. The TALYS numerical-reference gate in the results document must
be resolved before large distributional comparisons.

## External Linux server: runs longer than 30 minutes

Transfer the repository (including these uncommitted adapters if applicable),
verified data directory, shared CRPropa bundle if needed, and the frozen image
to a native Linux x86_64 server with Docker. Preserve the image ID by saving
and loading it, rather than rebuilding when resuming:

```bash
docker save mcray-crpropa-validation:local -o /path/outside-git/reference.tar
# On the external server after transfer:
docker load -i /path/outside-git/reference.tar
bash tests/crpropa/run_external_linux.sh /data/nuclei /results/decay-extended
```

The external script rejects macOS and non-amd64 hosts. It defaults to
extended: twenty batches of 10000 histories per code and scenario. The
server wrapper is the entry point for long checks. Before running TALYS,
prepare a reviewable plan without calculations:

```bash
python3 tests/crpropa/run_decay_comparisons.py \
  --data /data/nuclei --reference-data /data/shared-crpropa \
  --family talys --profile extended --output /results/talys-extended --plan-only
# After resolving the numerical-reference gate and freezing the campaign:
bash tests/crpropa/run_external_linux.sh /data/nuclei /results/talys-extended \
  --family talys --reference-data /data/shared-crpropa --resume
```

On interruption, use the same command with --resume. Paths may change when
moving hosts; input hashes, runner hash, profile, family and image ID must
match. Each successful job gets a checksum completion marker; partial jobs
are rerun with the frozen seed. Changed completed output checksums fail.
Existing pilot plans from an older runner require that exact runner to resume.

jobs.json records the matrix, seeds, image ID and data manifests. Per-job
logs, TSV particle rows, completion markers and comparison-diagnostic.json
are saved in the requested directory. A sampling unit is one primary
history, with independent batch means; individual secondary rows are
correlated. The runner reports integrated species yields, energies and
primary survival, with combined standard errors. It does not yet provide
binned-spectrum confidence intervals, rare-output adequacy or simultaneous
acceptance. Long-run output can be large; retain enough server disk space.

## Network rate convergence diagnostic

The native `nuclear_decay_probe --mode photo-rates` accepts `--network 1`
(the combined TALYS+decay network from `--a`/`--charge`), a comma-separated
`--gamma-list`, `--background-step` (log10 spacing, default 0.01) and
`--eps-rel` (default 1e-4). These are adapter-only numerical controls.
Default single-isotope behavior and production settings are unchanged.
Network scans retain failed isotope/gamma cells as `# FAILED` rows, continue
with other cells, and return exit code 3 if any cell failed. Such outputs are
incomplete diagnostics, never successful acceptance runs.
After a failed cell, the probe additionally evaluates its full and exclusive
components separately and records `# COMPONENT_FAILURE` channel identifiers.

```bash
build/Release/tests/nuclear_decay_probe --mode photo-rates \
  --data /data/nuclei --a 56 --charge 26 --network 1 \
  --photon-table /data/common-cmb.tsv --gamma-list 1e9,3.9810717055349693e9,1e10 \
  --output /results/network.tsv
python3 tests/crpropa/test_network_photo_rates.py
python3 tests/crpropa/diagnose_network_photo_rates.py \
  --data /data/nuclei --reference-data /data/shared-crpropa \
  --mcray-rates /results/network.tsv --output /results/network-diagnostic.json
```

The diagnostic verifies the shared reference manifest, integrates the agreed
photon-density interpolant analytically, and uses independent Gauss-Legendre
quadrature split at every cross-section node (orders 8, 16, 32 by default).
It retains full and every successful exclusive rate, failures, direct pinned
CRPropa convolution and its 201-node gamma-table interpolation as separate
controls. The latter is a rate interpolation diagnostic; CRPropa transport
selects branches at the nearest gamma node. Refined-order agreement is an
observed convergence indicator, not a rigorous error bound or independent
validation of the underlying nuclear cross sections. The photon-integral
unit checks use SciPy adaptive quadrature including zero-ended intervals.
These tools require NumPy/SciPy and the pinned preparation dependencies.


## Full gamma-domain network scan

`run_network_rate_scan.py` scans the combined Fe56 network on 401 points from
1e6 to 1e14, including the 201 official gamma nodes and every midpoint in
log(gamma). It checks every expected full, sum and exclusive row against the
prepared network, including terminal zeros. Missing rows, duplicate keys,
negative/nonfinite rates and failed cells invalidate the matrix.

The `convergence` profile runs background spacings 0.01 and 0.005 at tolerance
1e-4, plus spacing 0.005 at tolerance 1e-5, on both shared backgrounds. The
`default` profile runs only production controls. Both compare to the pinned
official convolution and independent orders 8/16/32. `--cache-photons` in the
reference diagnostic reuses identical photon-integral samples across channels;
it preserves the formula and was checked against the uncached path.
This is finite-grid numerical evidence, not continuous-domain proof or physical
acceptance. Absolute discrepancies of weak rates must be retained separately.

```bash
python3 tests/crpropa/run_network_rate_scan.py \
  --data /data/nuclei --cmb-reference /data/shared-cmb \
  --inoue-reference /data/shared-cmb-inoue --output /results/new-rate-scan
python3 tests/crpropa/test_network_rate_scan.py
```

Local execution has a total 1800-second timeout. Estimate the workload with a
short pilot before starting it; an expected longer run belongs on an external
server. The `extended` profile uses 801 gamma points, adds spacing 0.0025 at
1e-6 and independent orders 64/128; it requires external native Linux amd64.
No profile changes production physics or replaces the input cross sections.

For long runs, transfer an input bundle containing `nuclei/`, `cmb/`,
`cmb-inoue/`, the adapter files under `tests/crpropa/`, and
`tests/server/run_tests.py`. `inputs-manifest.json` must contain SHA-256 values
for every transferred file in `payload_sha256`, relative to the bundle root.
The wrapper verifies the bundle before cloning. Install the documented native
build dependencies and NumPy/SciPy/GitPython in the server Python environment.

CMake must be at least 3.24: `bin/build.sh` uses `cmake --fresh`. On servers
with an older system CMake, install a compatible version in the same virtual
environment used for the runner; no system package replacement is needed:

```bash
source /home2/ok/mcray-validation-venv/bin/activate
python -m pip install 'cmake>=3.24,<4'
hash -r
command -v cmake
cmake --version
```

The selected CMake should come from the virtual environment. Keep that
environment active when starting the wrapper, including inside tmux.

```bash
bash /srv/rate-inputs/tests/crpropa/run_external_network_rates.sh \
  /srv/rate-inputs /srv/mcray-test-data e43a8e3b623fc2e16a865a0e82013d5f4341617b
```

The wrapper fetches the requested Git revision at startup, records its exact
commit, builds Release, runs CTest, then executes the separately hash-frozen
extended diagnostic tools. The server runner retains successful or failed
results and compresses them with provenance into the selected test-data
directory. The workload output defaults to MCRAY_RESULTS_DIR/network-rates.
There are no random seeds: these are deterministic rate calculations. The
server archive must be inspected before reporting physical agreement.


`summarize_network_rates.py` keeps full totals, exclusive channels and their
sum distinct, reports absolute discrepancies below the diagnostic 0.1/Mpc
resolution threshold, and compares each numerical refinement. It retains
zero/nonzero disagreements and never replaces a weak rate by a successful zero.
`test_summarize_network_rates.py` guards these reporting semantics.

`diagnose_rate_boundary.py` checks the finite-upper-limit term in the
integration-by-parts representation used by the pinned official convolution.
It compares the unchanged official result, independent interval quadrature,
and refined diagnostic quadrature plus the endpoint term. Its default cases
are Be11, Be9 and Fe56. This is an explanation/control of reference differences,
not a replacement production rate or a change to the zero-padded tail model.
Include both scripts and their test files in transferred tool bundles.

## Frozen initial-point acceptance

[ACCEPTANCE.md](ACCEPTANCE.md) gives the frozen PD/decay/photon contract,
384 simultaneous comparisons, history ledger and external multicore commands.
`run_external_transport_photon_acceptance.sh` fetches a published exact source,
builds/tests it, checks the retained native deterministic coefficients and runs
272 fixed independent batch jobs. A wrapper `passed` status is execution
success; physical conclusions are in `acceptance-report.json`. Insufficient
precision remains inconclusive. The local `--plan-only` path starts no histories.


## Independent weighted full-Fe56 method pilots

[RARE_TRANSPORT_PILOTS.md](RARE_TRANSPORT_PILOTS.md) documents the validation-only
native and pinned-CRPropa weighted adapters, fixed independent pilot matrix,
preflight controls, persistent history reservation and external Linux commands.
These pilots determine variance/cost and method applicability; they cannot
issue final physical acceptance or replace the unchanged production models.
