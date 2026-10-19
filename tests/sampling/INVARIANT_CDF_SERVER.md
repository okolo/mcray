# Numerical log-CDF invariant sampler: bounded server comparison

This experiment compares the numerical inverse CDF with the frozen rejection
sampler from `1cc28f381e961f3d4947d08f0149ac72d63c3635`. Both use the corrected
ICS/GammaPP conditional secondary spectra. The independent secondary-spectrum
oracle, primary energies, log-energy bins, endpoint regions, seed namespace,
paired-product inference and acceptance targets are unchanged.

The default labels are `old` = frozen rejection and `current` = numerical
log-CDF inversion. The archive records this mapping in `sampler-labels.json`,
the frozen plan, analysis input and report. Raw control TSV headers identify
`sampler=rejection`. The old label no longer means NR by default. To repeat
the historical NR comparison instead, explicitly use `--baseline nr`.
The NR source remains pinned to `d05723f449c1eecee4f00af392014df4c898ebe9`;
its headers are downloaded at the existing verified revision when needed.

## Run on Ubuntu 22.04 / Python 3.10

The companion package contains a Git bundle with the committed candidate,
both frozen control revisions and every test tool. Publication is unnecessary.
Extract the package, enter its directory, then run:

```sh
PYTHON=python3.10 bash sampling/run_fixed_server.sh \
  --baseline rejection \
  --data-dir /home3/ok/mcray-test-data \
  --hours 24 --jobs 4 --raw-gib 20
```

The wrapper installs the pinned Python 3.10 requirements in `.venv`. System
dependencies are Git, CA certificates, CMake, make, GCC/G++, gfortran, GSL,
Boost filesystem/system development packages, Python venv, tar and gzip.
Use at least 42 GiB free space for the 20 GiB raw-event setting. See
[SECONDARY_SERVER.md](SECONDARY_SERVER.md) for budget and archive details.

After publication, the same runner accepts a remote repository and exact
candidate commit through `--repo` / `--ref`. The exact resolved commit is
always frozen at startup. Run from the candidate's `tests/sampling` tools;
an older portable adapter does not know the CDF and rejection-control hooks.

For a short infrastructure smoke check, add
`--smoke --hours 0.25 --raw-gib 0.01`. Smoke runs cannot establish acceptance.

## Unchanged statistical acceptance

- Core bin relative accuracy: 1%; tail relative accuracy: 10%.
- Global significance level: 1%, with the same Holm correction and
  simultaneous exact confidence intervals as the previous experiment.
- Full run: 11 cases, 1362 regions, 12258 category hypotheses. Overlapping
  views and paired products retain their original multiplicity treatment.
- No data-dependent bin merging, tail deletion, sample-size adaptation after
  the pilots or significance-based early stopping.
- The 24-hour total budget covers fetch, compilation, unit tests, deterministic
  controls, pilots, frozen ensembles, analysis, checksums and compression.
  Sampling is sequential; `--jobs` controls building only.

## Additional deterministic accuracy and timing records

Before the pilots, each physical invariant kernel is checked at 19 frozen
quantiles spanning both tails down to probability 1e-12. Prefix and suffix
probabilities are computed independently with composite fixed Gauss-Legendre
quadrature on 8192/16384 log intervals, orders 16/24. The production method
uses adaptive Gauss-Kronrod quadrature and safeguarded Newton inversion.

The extra numerical control is a 1e-3 relative small-tail probability error,
with a 1e-4 coarse/fine reference-change guard and 1e-4 normalization control.
These are deterministic implementation diagnostics, not an alternative
Monte Carlo acceptance criterion. Unresolved controls stop the run and retain
the data. Their convergence guard is empirical, not a rigorous error bound
for an arbitrary callable density.

Each kernel also records three interleaved timing blocks of 128 draws for
CDF inversion, frozen rejection and historical NR. Every CDF draw includes
rebuilding its normalization and inversion; no prebuilt-CDF shortcut is used.
The timing loop excludes event-file writing. The ordinary pilots additionally
measure the complete secondary-event path and determine attainable statistics.
The same maximum runtime can yield fewer events with the more expensive CDF.

Inspect `invariant-diagnostics/*.tsv`, `invariant-diagnostics.json`,
`pilot-timing.json`, `expected-statistics.json`, `report.json` and `report.tsv`.
Raw events, seeds, toolchains, source bundles, checksums and failure logs are
retained as before. `deviation` and `inconclusive` must be interpreted separately
for candidate, control and direct comparisons.

## Production scope

Background invariant sampling uses `MathUtils::SampleLogDistributionCDF`.
It integrates both tails directly and does not narrow the support, substitute
argmax or retry a numerical failure as another event. Known interpolation
knots are propagated through the background kernels in physical coordinates.

The old generic `SampleLogDistribution` remains for the pre-existing noisy
nested proton-pair secondary PDF. Linear-coordinate `SampleDistribution` is
also unchanged. This correction does not claim to validate or replace those
separate conditional samplers.
