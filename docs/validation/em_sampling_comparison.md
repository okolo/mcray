# Historical EM sampler accuracy comparison

The completed large-ensemble server run is reported in
[Server PP/ICS spectrum results](em_sampling_server_results.md). The sections
below retain the earlier exploratory controls and adapter preparation.

## Scope and frozen sources

Diagnostic date: 2026-10-12. Current public code: `cb5b2d9210846c6ba6b46ff2c4bb80ffb47b9751`
(`talys`). Historical sampler: `d05723f449c1eecee4f00af392014df4c898ebe9`.
The user's `/Users/ok/git/nr` checkout is clean at
`0e6af5083b426ca9423bc6da60372eba49f6ab6c`, exactly the historical NR gitlink.
The transition commit is `f62e40fcb1a62fd26377a2c416b6b9332499e8a9`.

The active production path after that transition is **rejection sampling with
GSL quadrature**, not the Boost ODE sampler retained in MathUtils.cpp.
In an electromagnetic cascade it changes the conditional distribution of
interaction invariant s. ICS and GammaPP conditional secondary-energy laws
at fixed s are unchanged between the historical and current revisions, apart
from particle identity assignment. Nuclear processes are irrelevant here.

The adapter extracts the historical NR function from Git, removes only its
optional logger and substitutes it into current Background.cpp in a generated
build file. It keeps current rates, background integration, secondary laws and
transport code identical in both modes. This measures the sampler's causal
effect. It is **not** a whole-program comparison of the two revision trees.
No production physics or sampler was modified.

## Independent CDF measurement

CMB: 2.7255 K, z=0, background energy range 1e-3--1e3 kT, spacing 0.01 decades.
Hard ICS photons: Emin=1e12 eV. Sampler tolerance: **1e-3**, the CRbeam default.
The oracle integrates the *same background-convolved production kernel*,
using composite 16-point Gauss-Legendre quadrature in log(s), independent of
both the adaptive NR ODE and GSL QAG. It does not independently validate the
underlying cross sections or background interpolation.

NR is evaluated at 2048 stratified uniform quantiles. Its CDF residual is
`max |F_reference(s_NR(u))-u|` at those points; it is a lower bound on the
continuum maximum. Current rejection sampling uses 200000 draws per kernel
with seed 20261012. Its number is an empirical one-sample KS distance,
which includes Monte Carlo noise and must not be read as deterministic bias.

| Process | Primary E [eV] | NR CDF residual [%] | Current empirical KS [%] | NR integral error [%] | GSL integral error [%] |
| --- | ---: | ---: | ---: | ---: | ---: |
| ICS | 1e+13 | 2.1531 | 0.1779 | -0.12960 | +0.003218 |
| ICS | 1e+14 | 2.1674 | 0.1492 | -1.02930 | -0.000596 |
| ICS | 1e+15 | 1.6041 | 0.2524 | -0.26058 | +0.000152 |
| ICS | 1e+16 | 1.4843 | 0.1825 | -0.49488 | +0.000095 |
| ICS | 1e+18 | 1.1856 | 0.2202 | -0.20387 | +0.000083 |
| ICS | 1e+20 | 1.8515 | 0.1441 | -1.14142 | -0.000190 |
| PP | 1e+14 | 0.9651 | 0.2388 | -1.34854 | -0.001258 |
| PP | 1e+15 | 2.7475 | 0.1500 | -2.07704 | +0.000012 |
| PP | 1e+16 | 1.7286 | 0.1883 | -1.41623 | +0.000941 |
| PP | 1e+18 | 1.6617 | 0.1588 | +1.04627 | +0.000335 |
| PP | 1e+20 | 0.9633 | 0.1674 | -0.29937 | -0.000347 |

Doubling the oracle grid from 16384 to 32768 cells changes the NR CDF
residuals by at most 1.41e-9 and relative normalization errors by 2.49e-9.
The results above are stable far below their reported precision.

The simultaneous 95% DKW sampling bound for 11 kernels and N=200000 is
0.00390089 (0.3901 percentage points). Every current KS distance is below it;
the largest is 0.00252411. Thus these data show no statistically resolved
current-sampler bias. Conservatively, they bound the underlying current CDF
error by 0.006425 simultaneously, conditional on the oracle and RNG
assumptions. They do not establish 1e-4 sampling accuracy.
At the default tolerance, all tested NR CDF errors exceed that conservative
current bound. The replacement improves accuracy in this tested family.

The largest historical discrepancy occurs for GammaPP at 1e15 eV:
NR CDF residual 2.7475%, integral error -2.0770%; current empirical KS 0.1500%,
GSL integral error +1.16e-5%. At the tighter tolerance 1e-4, the initial
GammaPP diagnostic still gives NR CDF residual about 2.4643% and integral
error -1.6718%. Merely tightening the ODE tolerance does not automatically
control inverse-CDF interpolation error or the error over table breaks.

**Integral errors in this table are internal sampler normalization errors.**
The propagation engine uses separately evaluated rates for the free path;
it does not use the normalization returned while generating s. Do not
interpret the table as an equivalent change in the mean free path.

## Remaining rejection-sampler limitations

The quadrature's largest evaluated log-PDF value is not a certified envelope.
Dense scans found underestimates of up to 0.7378% in this family. This is an
envelope-height discrepancy, not a measured CDF bias; the maximum can also
increase while rejection attempts are made. The code additionally narrows
the range after 1000 unsuccessful attempts and eventually returns argmax.
Neither mechanism has a general unbiased-sampling guarantee. No argmax
fallback warnings occurred in the completed kernel ensembles. These results
support improvement for the tested kernels, not universal correctness for
arbitrary sharply peaked distributions.

## Cascade pilots and limitations

Both modes use production GammaPP, hard ICS and IcsCEL, without thinning,
EBL, magnetic fields or nuclear processes. Primaries are 1e15 eV photons,
Emin=1e12 eV, source z=0.01, L=0.1 Mpc. The CMB follows native redshift
scaling along this short path. No separate redshift energy-loss process is
added. Per-primary bin moments preserve descendant correlations.

The default sampler tolerance is 1e-3. The short-path production
constant-rate mode (`Run(0)`, CEL fractional step 0.01) is shared by both
modes. Step convergence has not been established for percent-level spectral
claims. Independent seeds are used for the difference standard errors.

| Ensemble | Old NR | Current rejection |
| --- | --- | --- |
| 100-primary pilot | seed 20261013; 8.42 s; completed | seed 20261012; 7.94 s; completed |
| 2000-primary comparison | seed 20261020; 164.05 s; completed | seed 20261120; 163.63 s; completed |

No histories were discarded. In the 2000-primary comparison the largest
absolute count z-score is 2.372, in the lowest photon bin (1e12--10^12.25 eV):
old 15.479 and current 15.009 photons per primary; difference
-0.470 +/- 0.198 photons per primary (one standard error). This is a -3.04%
observed difference relative to the old mean, with about 1.28 percentage points
of statistical error. It is an exploratory, uncorrected bin comparison.
The largest energy-sum z-score is 1.632. There is no convincing family-wide
spectral discrepancy in this short ensemble, but it does not establish
percent-level spectral equivalence or absence of a small systematic effect.
Rare lepton bins are especially statistics-limited. Larger independent seeds
and a propagation convergence control belong in the supplied server workflow.

Two preliminary adaptive-step pilots were interrupted for runtime and are
not spectral evidence. At tolerance 1e-4, a separate current 100-history
constant-rate pilot finished in 9.62 s, while the old seed-20261013 pilot hit
a 120-second cap without a complete spectrum. This is not a controlled speed
benchmark and is not used to infer the default-tolerance spectral difference.
A 1e-5 kernel scan was also interrupted after its ICS controls; no GammaPP
result is claimed for that scan. No workstation calculation exceeded 30 minutes.

## Reproduction and platform verification

- macOS arm64: Homebrew Clang 22.1.5, GSL 2.8; 129/129 short CTest checks passed.
- Docker Linux amd64 emulation: GCC 13.3, GSL 2.7.1; 129/129 short CTest checks passed.
- The optional NR probe was independently compiled with GCC inside the Linux
  container and all 11 default-tolerance kernels were checked (N=2000 for its
  Monte Carlo smoke control). Maximum macOS/Linux difference in the deterministic
  NR CDF residual: 1e-13. This small Linux ensemble is not a second
  high-statistics rejection-sampler validation.

Runnable adapter and external Linux workload:
[tests/sampling/README.md](../../tests/sampling/README.md).
Local raw summaries, logs, source/header hashes and run status remain outside
Git in the adopted workspace's `validation/em-sampling-20261012/` directory.
The server workflow freezes the requested revision, runs tests and archives
results to the configured test-data directory. No remote ensemble was launched.
No commits or pushes were made for this diagnostic.

## Bounded secondary-spectrum server test (2026-10-12)

The implemented [PP/ICS server adapter](../../tests/sampling/SECONDARY_SERVER.md)
uses the user-approved 1% core / 10% tail criteria and global 1% significance.
It freezes 26 nonuniform energy-fraction bins, cumulative tails, exact numerical
references and sample plans before main sampling. Event multiplicity retains
paired-product correlations. A total 24-hour deadline includes build, analysis
and archiving; insufficient rare-tail precision remains inconclusive.

Preparation checks: all 132 CTest cases pass on macOS arm64 and Linux amd64
Docker emulation, including three oracle checks. Twelve Python checks cover
exact inference, tail bias, zero-count limits, incomplete-run rejection,
process termination, deterministic seeds and retention after preflight failure.
All 11 physics cases were exercised with three reference grids and 1000 events
per sampler on macOS (22000 interactions total). Independent rebinning of
saved binary events reproduced every multiplicity count table. These short
runs validate the adapter and data format, not the approved precision targets.
Two complete Linux infrastructure smoke runs fetched/froze the requested revision,
built, tested, sampled PP/ICS, analyzed and compressed outputs in 251 and 235 s;
both deliberately small ensembles reported inconclusive as required.

Short preflight data are retained outside Git in
`validation/em-secondary-preflight-20261012/` of the adopted workspace.
No 24-hour calculation or remote ensemble was started. Production sampling
and physics files were not modified by this test preparation.

The server runner now downloads NR automatically from
`https://github.com/blackstonep/Numerical-Recipes.git` when `--nr` is omitted.
It checks out the historical NR gitlink exactly, validates all four header
hashes before building and preserves download logs/provenance. An actual
remote download reproduced the expected commit and header hashes. Local
`--nr` overrides remain available. Fifteen Python checks pass, including
pinned-revision selection despite a moving branch and checksum failure
retention. Python tests also pass in Ubuntu 22.04 amd64 Docker with Python
3.10.12, NumPy 2.2.6 and SciPy 1.15.3; a separate pinned requirements file
is provided for that environment. This does not constitute a full C++ build
verification on Ubuntu 22.04 or a long acceptance calculation.
