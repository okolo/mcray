# Frozen initial-point transport/photon acceptance

This executable protocol was frozen on 2026-10-20 after independent pilots
and the audited I/M/C coupled numerical matrix. It evaluates Fe56 at
E/A=10^18.6 eV, 1 physical Mpc, frozen z=0 CMB and CMB+Inoue12, with TALYS,
corrected nuclear decays and both photon sources. GZK, PPP, free-neutron decay,
redshift losses, magnetic deflection, thinning and EM transport are disabled
in these isolated acceptance adapters. Full-process CLI smoke controls are
separate. This is not acceptance across other energies, distances, fine
nuclear energy spectra or the complete Release feature set.

The canonical immutable values are
[protocols/transport_photons_v1.json](protocols/transport_photons_v1.json),
SHA-256 `e5f345a2773745790f7c59e6941b143cb70450e80a42af38a28aac11deef76a9`.
Do not change bins, seeds, samples, normalizations or stopping rules after
examining final results. A protocol revision requires a new decision and
accounting for previously executed histories, not another attempt until pass.

## References and observables

Compare mcray MC with deterministic M and unchanged official CRPropa MC with
deterministic C. The cached reference campaign is native Linux run
`20261019T205338Z-943099d8`, source `7fc04ae6ce2e76887287ba85562b721d6b31ddf8`.
Retain its full result payload. The protocol freezes input/result hashes,
N65 Radau target means and observed grid/path envelopes. I order64 means set
fixed normalizations. M/I is a separate coefficient check; C/I is an
unmodified representation diagnostic, never an empirical MC correction.
The runner checks ten current-production coefficient exports against the
cached native exports before generating any new primary histories. The
5e-12 relative comparison is a floating-point identity guard, not a new
physical acceptance tolerance; unexpected differences stop the campaign.

There are 32 observables per background/program: Fe56 survival; bound-nucleus
count, sums of A and Z, and bound nuclear energy per primary; p/n counts and
energies; electron, positron, electron-neutrino and antineutrino counts and
energies; total photon count and energy; counts of A55Z26, A54Z26, A55Z25,
A54Z25 and A53Z25; and photon count/energy in each of the four fixed intervals
[0,1e14), [1e14,1e15), [1e15,1e16), [1e16,infinity) eV. The isotope choices
come from the independent pre-acceptance reference, not final MC samples.
All quantities are means per complete primary history; descendants inside a
history are correlated. Sums of A/Z are not ratios of random variables.
Separate PD/decay photon-source spectra are outside this frozen family.

The budget is 1% of the fixed absolute I mean for integral quantities and
5% for photon spectral intervals. No near-zero absolute budget is invented.
Constant/zero and insufficiently resolved samples remain inconclusive.
Use Student batch intervals, df=B-1, with Bonferroni family confidence 99%
over all 384 comparisons (32 observables, two backgrounds, two programs and
three numerical settings). With d=MC mean minus matching deterministic mean,
h=confidence half-width, epsilon=observed deterministic numerical envelope,
and delta=fixed tolerance:

- agreement requires |d|+h+epsilon <= delta and integral SE <=1% of normalization;
- disagreement requires |d|-h-epsilon > delta;
- otherwise the result is inconclusive.

Numerical envelopes are observed refinements, not rigorous global bounds.
MC refinement outputs have their own intervals; noisy differences are not
inserted as deterministic biases or used to correct MC means. Every frozen
comparison must agree for global agreement. The report separates the main
family, refinement controls and M/I checks; C/I differences do not reject the
MC implementation. Rare outputs can keep the overall result inconclusive.
The existing pilots project that some photon integrals may also remain
inconclusive under the cap; execution success must not be called acceptance.

## Samples and parallelism

Each scenario/program has 48 main batches of 2000 histories (96000), plus
10 batches at each of two finer settings (20000+20000). Together with the
60000 prior pilot executions and a conservative 4000-history reservation for
other prior controls, this respects the unchanged 200000 cap. There are
272 new subprocess jobs and 544000 primary executions in the complete
campaign. The reservation is a conservative ceiling, not a claim that all
4000 earlier controls actually occurred.

The main mcray background spacing is .0025 decades at eps-rel=1e-6;
refinements use .00125 and .000625 at the same tolerance. The corresponding
CRPropa step fractions are .03/9, .03/27 and .03/81. Seeds are fresh,
independent between backgrounds/programs/batches, and paired across the
three numerical settings. Each subprocess has OMP/BLAS/MKL threads set to
one. `WORKERS` bounds concurrent independent subprocesses; changing CPU
scheduling does not change the frozen job seeds or sample counts.

A persistent `transport-photon-acceptance-v1-ledger.json` reserves the new
campaign before the first MC job. Another output/source cannot reuse that
reservation. Keep the ledger with the test-data directory. Never delete it
and restart a full campaign after examining results or a partial failure.

## External Linux commands in /home2/ok

The new source must first be committed and pushed to the public remote.
Replace `PUBLISHED_COMMIT` below with its exact public commit. The launcher
archive is a bootstrap, not a substitute for publication of the source.
Use native Linux x86_64 with yaml-cpp and OpenSSL development packages
(`sudo apt-get install libyaml-cpp-dev libssl-dev`) and the existing
C/C++/Fortran, Boost/GSL and CMake build environment that completed the previous
matrix. Ubuntu 22.04's stock Boost 1.74, OpenSSL 3 and CMake 3.22 are supported;
the build helper does not require `cmake --fresh`.
Use Python >=3.10. The prior reference environment recorded NumPy 2.2.6,
SciPy 1.15.3 and GitPython 3.2.0. Reuse it, or create an environment:

```sh
python3 -m venv /home2/ok/transport-acceptance-venv
source /home2/ok/transport-acceptance-venv/bin/activate
python -m pip install numpy==2.2.6 scipy==1.15.3 GitPython==3.2.0
```

The launcher can be obtained directly from the published public repository.
Use a new directory for the exact revision supplied with the launch commands;
the wrapper creates a separate clean source checkout for the calculation:

```sh
git clone --branch talys https://github.com/okolo/mcray.git \
  /home2/ok/transport-acceptance-tools
git -C /home2/ok/transport-acceptance-tools checkout --detach PUBLISHED_COMMIT
```

If using the bootstrap archive instead, transfer it and the adjacent SHA-256
file, then extract:

```sh
cd /home2/ok
sha256sum -c transport-photon-acceptance-launcher.tar.gz.sha256
mkdir -p /home2/ok/transport-photon-acceptance-launcher
tar -xzf transport-photon-acceptance-launcher.tar.gz \
  -C /home2/ok/transport-photon-acceptance-launcher
```

Preserve `/home2/ok/rate-inputs` unchanged. In particular, use its original
nuclear photon bundle, matching the frozen initial network hashes; newly
prepared GZK photon bundles belong to a separate full-process control.
The two retained server result directories must exist:

```sh
test -f /home2/ok/mcray-test-data/work/20261016T102021Z-03eb6c57/results/transport-photons/reference-build/crpropa_decay_probe
test -f /home2/ok/mcray-test-data/work/20261019T205338Z-943099d8/results/deterministic-coupled/matrix.json
```

If those working directories were removed, unpack the original successful
archives into separate retained directories and adjust the corresponding two
arguments below. Preserve nested paths and executable modes. Do not copy
partial summaries in place of full reference payloads. The supplied pinned
reference binary requires the same server ABI as the successful pilot run.

Run in tmux/screen if desired. Choose eight CPU workers, or an appropriate
positive worker count for the server's CPU and memory:

```sh
source /home2/ok/transport-acceptance-venv/bin/activate
bash /home2/ok/transport-photon-acceptance-launcher/tests/crpropa/run_external_transport_photon_acceptance.sh \
  /home2/ok/rate-inputs \
  /home2/ok/mcray-test-data/work/20261016T102021Z-03eb6c57/results/transport-photons \
  /home2/ok/mcray-test-data/work/20261019T205338Z-943099d8/results/deterministic-coupled \
  /home2/ok/mcray-test-data \
  PUBLISHED_COMMIT \
  8
```

For the Git checkout above, use
`/home2/ok/transport-acceptance-tools/tests/crpropa/run_external_transport_photon_acceptance.sh`
as the script path. Other arguments are unchanged. Each subprocess uses one
OpenMP/BLAS/MKL thread; eight workers means eight independent concurrent batch
jobs, not eight repetitions of the campaign. Keep the process in tmux/screen
until the final archive is produced.

The wrapper fetches the exact public source, records its commit, builds,
runs CTest, verifies cached inputs/references, checks coefficient identity,
then runs the fixed batch list. It retains per-job commands, seeds, complete
histories, output hashes, logs, binaries' hashes, dependency versions and
status, and archives results even on failure. No extra histories are launched
after the fixed list to obtain a pass. The final archive's `passed` label means
execution success. Read `transport-photon-acceptance/acceptance-report.json`
for agreement/disagreement/inconclusive and scope.

For an interruption, retain the source checkout, results, ledger and partial
files. Resume in that exact checkout, with its original environment and
arguments, adding `--resume` and the existing `--output` to the recorded
workload command; preserve the recorded `MCRAY_RESULTS_DIR`. This verifies
source, protocol, inputs and binary identity and skips hashed completed
batches. Complete unmarked files can be recovered after validation. A partial
batch with unknown consumed histories stops; it is never silently rerun.
After a manual resume, re-archive the retained result directory separately;
the original failed archive remains unchanged evidence.

A local `--plan-only` mode validates inputs and writes the immutable 272-job
plan without executing transport. It does not reserve the server campaign or
constitute a clean-source/server-execution check.
