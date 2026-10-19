# Bounded PP / inverse-Compton spectrum comparison

For the new numerical log-CDF candidate, use
[INVARIANT_CDF_SERVER.md](INVARIANT_CDF_SERVER.md). Its default baseline is
frozen rejection at `1cc28f3`; select `--baseline nr` to reproduce the NR
comparison described below. Statistical targets and frozen regions are retained.

This implemented experiment compares the historical NR log sampler at
`d05723f449c1eecee4f00af392014df4c898ebe9` with the sampler in the requested
current revision (default `talys`). It also compares both with an independent
integration of the secondary differential energy laws. The approved targets
are **1% relative error in core bins, 10% in tails, global significance 1%**.
This is a single-interaction test; no nuclear data or full cascade is required.

## External server invocation

For the unpublished conditional ICS/PP correction, use the companion Git
bundle and wrapper described in the
[correction report](../../docs/validation/em_sampling_fix.md). The frozen
protocol below is retained without changing bins, seeds or criteria.

Use a native Linux amd64 server with Python 3.10 or newer, Git, CMake 3.19 or newer, make,
GCC/G++, gfortran, GSL and Boost filesystem/system development packages,
CA certificates, tar and gzip. CMake downloads checksum-pinned GoogleTest
1.18.0 unless an offline FetchContent source override is provided. Python
requirements for Python 3.12+ are pinned in `secondary-requirements.txt`.
For Ubuntu 22.04 / Python 3.10, use `secondary-requirements-py310.txt` instead
(NumPy 2.2.6, SciPy 1.15.3). Actual versions and `pip freeze --all` are recorded.
When `--nr` is omitted, the runner clones
`https://github.com/blackstonep/Numerical-Recipes.git` and checks out the pinned
commit below, without following a moving branch. To use existing headers or
run offline, provide `--nr /absolute/path/to/nr`. Four SHA-256 checks enforce
the historical NR revision
`0e6af5083b426ca9423bc6da60372eba49f6ab6c`; no NR implementation is distributed
with this adapter.

Copy the **entire `tests/sampling` directory** to the server. Its portable
adapter overlays only test files on the fetched revision; it therefore works
before these tools are pushed. The physics revision still comes from Git.
No access to a private workspace repository is needed.

```sh
python3 -m venv /srv/venvs/mcray-sampling
/srv/venvs/mcray-sampling/bin/pip install -r sampling/secondary-requirements.txt
/srv/venvs/mcray-sampling/bin/python sampling/run_secondary_server.py \
  --server --repo git@github.com:okolo/mcray.git --ref talys \
  --data-dir /srv/mcray-test-data \
  --hours 24 --jobs 4 --raw-gib 20
```

The paths above are examples; choose writable server directories. On Python
3.10 replace the requirements filename in the installation command with
`sampling/secondary-requirements-py310.txt`. The optional `--nr` must point to
a real directory containing `nr3.h`, `stepper.h`, `odeint.h` and
`stepperdopr5.h`; without it NR is downloaded automatically within the same
total runtime budget. Download failures and checksum mismatches retain logs
and fail before the physics build.

Set `--ref` to an exact commit to repeat a frozen run. At startup the script
fetches the requested ref, records its resolved commit, creates source tar
snapshots and a self-contained Git bundle, builds the copied adapter and runs
CTest and Python statistical/runner tests. It rejects incompatible adapter
hooks rather than silently adapting to changed production code. Ordinary
builds do not require NR or SciPy.

The **24-hour total budget includes fetch, build, tests, reference calculation,
pilots, sampling, analysis, checksums and compression**. The default work
deadline is 23 hours; the last hour is reserved for analysis and archiving.
All subprocess groups have deadlines; a total process alarm also interrupts
Python work. Timed-out subprocesses receive TERM then KILL (up to two seconds
for shutdown). A five-second guard is withheld from the total budget for final
process termination and status writing. A deadline or disk failure retains the complete working
directory and all partial data, with explicit failure status. It cannot
promise a finished archive if storage fails or the deadline is exhausted.
The script never declares a partial run successful.

Sampling is sequential; `--jobs` controls building. The independent 1000-event
pilot per case/sampler only estimates runtime. It is excluded from inference.
The main counts are frozen from pilot timings, remaining time (65% safety
allocation) and the raw storage cap. `expected-statistics.json` records expected
category counts and approximate simultaneous precision before the main run;
these normal approximations are planning diagnostics only. Each case gets equal estimated runtime
and equal old/current event counts. Frozen chunks have at most 100000 events,
with distinct, reproducible 32-bit seeds. There is no significance-based early
stopping or post-sampling bin merging. The default raw cap is 20 GiB; require
at least twice that free space plus 2 GiB for build/source overhead. A lower
`--raw-gib` reduces the attainable statistics. Very rare tails can remain
inconclusive within 24 hours.

For infrastructure verification only, add `--smoke --hours 0.25 --raw-gib 0.01`.
This uses two cases, reduced reference grids and 1000 main events per sampler;
it always disables acceptance claims. Never run a workload expected to exceed
30 minutes on a workstation, including Docker emulation.

## Frozen observables and numerical reference

Both modes use current production CMB convolution, cross sections and
conditional secondary sampling. Only the implementation selecting invariant
`s` is replaced. This isolates the sampler update; it is not an end-to-end
build comparison of every change between the two Git revisions. The shared
cross sections and photon-background kernel are not independently validated
by this experiment.

Configuration: CMB 2.7255 K, z=0, background photon energies k/T in
[0.001,1000], grid spacing 0.01 decades, sampling tolerance 1e-3 and hard ICS
photon cutoff 1e12 eV. Electron energies are 1e13, 1e14, 1e15, 1e16, 1e18,
1e20 eV; primary photon energies are 1e14, 1e15, 1e16, 1e18, 1e20 eV.
The negligible CMB pair-production case at 1e13 eV is omitted.

There are 26 nonuniform bins in `x=E_secondary/E_primary`: logarithmic at low
x, logarithmic in `1-x` at high x and wider bins centrally. Exact edges are
saved in `edges.txt`. Additional predeclared cumulative lower/upper tails use
cuts 1e-2, 1e-3, 1e-4, 1e-6, 1e-8, plus PP endpoints near x=1/2.
PP views are the lower/higher-energy lepton, electron, positron and inclusive
lepton spectrum. ICS views are the photon and recoil electron. For each
interaction and region, counts retain multiplicity **0, 1 or 2**. The two
products are never counted as independent trials.

The oracle directly integrates the positive differential laws in
[ELMAG, equations (11), (12)](https://arxiv.org/pdf/1106.5508), with composite
Gauss-Legendre quadrature in log lepton energy. It integrates disjoint event
patterns, including both PP charge assignments, rather than subtracting CDFs
near one. Its outer integration uses the production background-convolved
`s` kernel but independent log(s) quadrature. It never calls production
secondary primitives, inverse CDFs, root solvers or either tested sampler.
Moving kinematic boundaries split quadrature cells explicitly.

Reference grids (128/16, 256/24, 512/32 outer cells/quadrature order) are
calculated and saved before each case pilot and the main sample. The numerical error
envelope is four times the largest successive change, with a 2e-14 absolute
floor per multiplicity probability. This is a convergence diagnostic, **not
a mathematically certified error bound**. Acceptance against the oracle
requires mean uncertainty below 5% of the allowed bias and 10% of the
predicted Monte Carlo standard error. Unresolved probabilities remain
inconclusive. Bins entirely outside the reference central 98% of a view are
classified as tails; all other bins use 1%. Extra cumulative tails always use
10%. This classification is frozen before sampling.

## Statistical decisions

There are 12258 predeclared category tests in the full experiment: three
multiplicities, every region/energy and three comparisons (old/oracle,
current/oracle, old/current). The whole family is retained, including
structurally empty categories. Exact equal-tail binomial tests maximize over
the oracle numerical interval; direct old/current tests use Fisher's exact
conditional test. Holm correction controls family-wise significance at 1%
under arbitrary correlations, including overlapping tails and paired products.

Simultaneous Bonferroni Clopper-Pearson category intervals are combined into
bounds on event yield `p1+2*p2`. Old-minus-current intervals combine the two
independent ensemble bounds. Relative differences use the oracle mean as
the common denominator. `within_target` requires the **whole simultaneous
interval** inside the agreed tolerance; `deviation` requires it entirely
outside. Otherwise the result is `inconclusive`. Statistical significance is
reported separately, including significant effects smaller than the tolerance.
An unobserved but physically allowed tail retains a positive upper limit.
Support-proven constant categories are assigned their exact probabilities;
forbidden regions are checked deterministically; rounded-zero
reference probabilities alone do not prove a forbidden region.

Methods: [SciPy exact binomial intervals](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats._result_classes.BinomTestResult.proportion_ci.html),
[Fisher exact test](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html),
[Holm correction](https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html).

## Outputs, failures and reproduction

The timestamped directory and `.tar.gz` archive retain source snapshots/Git
bundle, exact tools, NR headers privately, their hashes, generated C++ copies,
binary, CMake cache, compiler/library/Python versions, commands, NR download
log and `nr-provenance.json`, reference
grids, frozen reference/plan, seeds, raw pilot/main events, per-chunk counts,
logs, progress and failure status. `sha256.json` covers archived inputs and
results; `final-status.json` is retained in the directory; an identical `.status.json`
sidecar and the archive log remain beside the archive. The post-archive status
is not part of the archive itself.
Do not publish the archive containing licensed NR headers without checking
those license terms. The uncompressed directory is deliberately retained.

`report.md` summarizes results; `report.tsv` contains every bin/tail yield and
simultaneous interval; `report.json` also includes every unadjusted/adjusted
p-value. `diagnostics.json` preserves proposal/envelope/narrowing/argmax counts
and maximum kinematic/conservation residuals. Any nonfinite or out-of-support
energy, invalid product identity, lost event, conservation violation or output
write error fails sampling. A failing event's values and index are logged,
and any available product record is preserved. Envelope exceedances and
fallbacks are recorded for diagnosis without changing sampling decisions.

Raw `.events` format: magic eight bytes `MCRAYEM1`, then 26-byte little-endian
records `<dddBB`: s/m_e^2, first energy fraction, second energy fraction,
first particle type (0 electron, 1 positron, 2 photon), audit bits (1 envelope
exceeded, 2 range narrowed, 4 argmax returned). PP first product is the lower
lepton; ICS first product is the photon. Event identity is chunk plus record
index, with the chunk's frozen seed. Per-event s and both fractions allow
independent rebinning/diagnosis; the original family remains the acceptance
family. Raw files from failed chunks may be partial and are not analyzed.

Reanalyze without repeating MC:

```sh
python adapter/analyze_secondary.py analysis-input.json report-repeated.json
```

For a source rebuild, import the two explicit bundle refs:

```sh
git init replay
git -C replay fetch ../frozen-history.bundle \
  refs/reproducibility/current:refs/heads/replay \
  refs/reproducibility/historical:refs/heads/historical
git -C replay checkout replay
```

Verify the recorded current commit, copy `adapter` into `tests/sampling`, append the include of
`secondary_targets.cmake` as in the launcher, and configure with the archived
`nr-private` directory. Use the archived generated sources to inspect exactly
which historical function and audit instrumentation were built. Bitwise RNG
reproduction across compilers/library versions is not promised; preserve the
recorded environment or the saved events for exact reanalysis.

Exit codes: 0 all frozen precision checks established; 2 material deviation;
3 inconclusive (including smoke); 124 deadline; 1 execution/archive failure.
Missing tail precision is a result, not a justification to loosen the agreed
criterion. A timeout before a complete frozen ensemble preserves raw data but
does not issue acceptance or significance claims from a selected partial run.
