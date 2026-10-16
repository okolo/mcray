# Historical electromagnetic sampler comparison

For the approved 1% / 10%, global 1% single-interaction experiment with a
24-hour total server limit, use [SECONDARY_SERVER.md](SECONDARY_SERVER.md).
The exploratory kernel/cascade tools below have a different purpose.

This opt-in diagnostic substitutes the exact log-distribution NR sampler from
`d05723f449c1eecee4f00af392014df4c898ebe9` into the current transport code.
Production source files are never edited. Both modes use current background
integrals, cross sections, secondary-energy laws, rates and propagation controls.
This isolates the sampler change; it is not a full executable comparison of
all changes between the two revisions.

Supply Numerical Recipes headers separately. They are neither bundled nor
needed for ordinary builds/tests. The historical revision must be available
in local Git history (fetch it if using a shallow checkout).

```sh
cmake -S src/app/crbeam -B build/Release \
  -DMCRAY_NR_REFERENCE_DIR=/path/to/licensed/nr
cmake --build build/Release --target em_sampling_probe -j 4
build/Release/tests/em_sampling_probe kernels 200000 16384 1e-3 > kernels.tsv
build/Release/tests/em_sampling_probe cascade old 100 20261012 old.tsv 0.1
build/Release/tests/em_sampling_probe cascade current 100 20261013 current.tsv 0.1
python3 tests/sampling/summarize.py old.tsv current.tsv
```

Use the documented platform build configuration first. On macOS, initialize
with `bash bin/build.sh Release` before enabling the optional target.
`kernels N grid tolerance` uses a CMB at 2.7255 K, z=0, photon-background
spacing 0.01 decades and a hard ICS photon threshold of 1e12 eV.
The diagnostic checks six electron energies (1e13--1e20 eV) and five photon
energies (1e14--1e20 eV); the negligible CMB pair rate at 1e13 eV is skipped.

The reference CDF uses composite 16-point Gauss-Legendre quadrature in log(s),
including direct integration inside the final CDF cell. Check convergence by
increasing `grid`. The old sampler is tested at 2048 stratified quantiles;
its reported maximum CDF residual is a lower bound on the continuum maximum.
The current sampler uses N random samples and reports a one-sample KS distance.
These two errors have different statistical meanings. For K kernels, the
simultaneous 95% DKW sampling bound is `sqrt(log(2*K/0.05)/(2*N))`.
No physical acceptance threshold is imposed by this exploratory tool.

`envelope_shortfall` is the ratio of the largest log-PDF value found on a
4*grid scan to the maximum found during the sampler's GSL61 quadrature,
minus one. This identifies an imperfect rejection envelope; it is not the
CDF bias itself. The sampler can update its maximum during rejection attempts.
Logs must be inspected for range-narrowing and argmax fallback warnings.
`secondary_fraction_mean` uses the current sampled s and the production
conditional secondary law; it is descriptive, not an old/new comparison.

The cascade uses production PropagationEngine, GammaPP, hard ICS and IcsCEL,
with E0=1e15 eV photons, Emin=1e12 eV, source z=0.01, CMB and configurable
path <=1 Mpc. No thinning, magnetic field, EBL, nuclear interactions or
redshift energy-loss process is enabled. The CMB evolves along the short path.
The sampler/integration tolerance is the CRbeam default 1e-3 and the standard constant-rate propagation mode is used (`Run(0)`, default
CEL fractional step 0.01). This is a short-path exploratory diagnostic, not a
precision-converged propagation reference. An optional final argument sets
maximal rate variation for a slower propagation convergence control.
The two modes consume different random-number sequences even at equal seeds.
Use independent seeds for reported difference standard errors.

Each cascade output records sums and sums of squares of **per-primary**
counts and energies in 0.25-decade bins, separately for photons and combined
leptons. Empty primary bins contribute zeros. Descendants are correlated and
must not be treated as independent samples. `summarize.py` reports count and energy
means, independent-ensemble errors and descriptive z-scores; no multiple-test
acceptance decision is made. Nonzero propagation failure counts fail the run.

Generated source provenance and NR header hashes are in
`build/Release/tests/sampling_reference/manifest.json`. The extraction disables
only the historical optional function-call logger. Review the generated source
when changing the historical revision or sampler signature.

## External Linux amd64 ensembles

Never run a calculation expected to exceed 30 minutes locally. Use the existing
server runner, which resolves/fetches the requested revision, freezes its commit,
builds, runs CTest and archives logs, failures and results in the configured
server test-data directory. Publish the diagnostic code before requesting a
remote revision that includes it; an unpushed working tree is not available.

```sh
export MCRAY_NR_REFERENCE_DIR=/srv/licensed/nr
export MCRAY_CASCADE_N=10000
python3 tests/server/run_tests.py --server --ref talys \
  --data-dir /srv/mcray-test-data -- \
  bash tests/sampling/run_comparison.sh
```

The workload uses old seeds 20261012--20261014 and current seeds
20261112--20261114 and retains per-seed summaries. Equal seeds across
implementations would not give paired histories.
Set `MCRAY_RATE_VARIATION` to a positive value (for example 0.01) for an
additional propagation-step control; compare it separately before interpreting
small spectral changes.

Large runs refine spectral sensitivity; the default sample size and seed set
are exploratory and do not establish a preapproved acceptance criterion.
