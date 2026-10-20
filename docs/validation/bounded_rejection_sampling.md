# Analytically bounded invariant rejection sampling

Branch `codex/bounded-rejection-sampling` starts at `fa562d2` on the CDF
candidate. Conditional ICS/GammaPP energy laws are unchanged. Continuous and
monochromatic invariant-S sampling now use analytic rejection bounds when
the cross section supplies them. Other invariant kernels retain CDF inversion;
the separate proton-pair conditional sampler is unchanged.

## Envelope contract

With `t=log(S/S0)`, the sampled density is `g(t)=S*f(S)`. Existing GAUSS61
quadrature computes normalization and an observed lower bound on the maximum.
An adaptive priority queue covers the entire physical support with intervals.
Every interval has an analytic upper bound U; sampled values never replace U.
The largest U is refined until `max(U)<=1.01*max(observed g)`. The proposal
remains uniform in t, with constant height M=max(U). Sampling uses g/M.

Two bounds are combined. Positive-factor bounds remain valid at thresholds,
zero-density endpoints, and interpolation knots. Away from zeros, if
`(log g)'' >= -K` on smooth pieces and the sum of downward log-slope jumps
is at most J, then on an interval of width h:

`g <= max(g_left,g_right)*exp(K*h*h/8 + J*h/4)`.

The smooth term is the linear-interpolation remainder; the jump term follows
from the Dirichlet Green function, whose maximum at a jump is h/4.
See [NIST DLMF 3.3.5](https://dlmf.nist.gov/3.3.E5).
Curvature is bounded analytically, without finite differences or fitting.
Known table slope jumps are accumulated explicitly, rather than treating
the interpolant as twice differentiable across its knots.

Arithmetic uses long-double intermediates, a 1e-10 relative evaluation guard
and outward conversion with nextafter. Background subtraction has an additional
additive cancellation guard. The mathematical inequalities enclose the physical
formulas; these guards are not a formal interval-arithmetic proof of every
platform's libm or of upstream numerical cross-section evaluation. Independent
envelope audits and runtime violations remain mandatory.

The interval construction is limited to 4096 refinements. Unsupported or
unresolved bounds select CDF before any proposal variate is consumed. EM
validation treats such a fallback as a failure and records counters.
A negative/nonfinite density, an invalid bound or `g>M` is fatal; M is never
silently increased after a proposal. Normalization errors are non-retryable,
so a failed numerical event cannot discard just one cascade history. Exact
physical endpoints are preserved; collapsed absolute-log coordinates select
CDF before quadrature instead of returning an empty interaction.
One million unsuccessful attempts also
fail the calculation. No support restriction or argmax return is permitted.

## Cross-section inequalities

For GammaPP, let beta=sqrt(1-4m_e^2/S), and write sigma=C*A(beta)/S.
The positive bracket is increasing: its series divided by beta is

`A/beta = 2 + 4 beta^2 - (4/5) beta^4 + sum_{n>=3} a_n beta^(2n)`,
`a_n=(8n-20)/((2n+1)(2n-3))>0`.

Its derivative is positive on (0,1). Thus
`sigma(S)<=S_hi*sigma(S_hi)/S_lo` on the whole interval. Differentiating
the same bracket provides a lower bound on log-sigma curvature. No finite
curvature is required at the pair-production threshold.

For hard ICS, use w=S/m_e^2-1, c=r_min/(1-r_min), y=1/(1+w*u),
and C=SigmaCoef/m_e^2. The implemented cross section has the positive form

`sigma/C = integral_{c/w}^1 F(w,u) du`,
`F = y^3+y-4u(1-u)y^2`.

For fixed u, F decreases with w, is bounded by 2, and its log-w derivative
is in [-3,0]. Consequently, for a<=w<=b:

`sigma(w) <= min(sigma(a)+2*C*c*(b-a)/(a*b), sigma(b)*(b/a)^3)`.

The corresponding lower bound is
`max(sigma(a)*(a/b)^3, sigma(b)-2*C*c*(b-a)/(a*b))`.
Also `F=y*(1-y)^2+(2-4u(1-u))*y^2` is a positive sum whose log-w
curvature is at least -3/4. Integration adds a nonnegative variance;
the moving lower-limit terms are bounded explicitly in ICS::Sigma.
The change from log w to log S and the kinematic factor S*(S-m^2) are
included in the final bound.

Positive cumulative background tables are nonincreasing. Their log-linear
segments are powers; positive redshift mixtures have nonnegative smooth
log curvature after composition with log(S-m^2). Downward slope jumps are
bounded using immutable prefix sums, rebuilt by Clone. If the finite-beta
subtracted integral is nonzero, only the factor bound is used. This avoids
assuming log-convexity of a difference. Negative/nonfinite background increments
are rejected at construction.

## Validation and reproduction

The physical diagnostic retains all 19 CDF-control quantiles, including both
tails down to 1e-12. The bounded envelope is checked on independent fixed
Gauss-Legendre nodes, at and next to table knots, and locally in every final
cell. Forced tail proposals check that the production acceptance path retains
nonzero density throughout the support. These deterministic checks do not
assert Monte Carlo precision for 1e-12-probability events.

Three interleaved 128-draw timing blocks include normalization and full envelope
construction on every draw. No cross-event cache is used. Native transport
profiling separately includes rates and continuous losses. Emulated amd64
times are never extrapolated to a physical server.

The fixed-size local spectrum runner samples rejection, CDF and bounded
rejection and saves every raw interaction, source, executable and seed.
Server criteria, commands and failure retention are described in
[the portable server guide](../../tests/sampling/BOUNDED_REJECTION_SERVER.md).
Confirmed measurements follow below; all raw results remain separate from Git.


## Confirmed short checks on 2026-10-20

The computational revision frozen in the server package is
`837d0108dd5995a2f7b3b080ba5ec28e7228f07d`. Its native macOS arm64 build
passes all 155 CTest entries. The exact package independently fetched pinned
NR, built from scratch on Ubuntu 22.04/Python 3.10 under amd64 emulation,
and passed all 155 C++ and 25 Python tests. All 25 Python tests also pass in
the dependency-equipped native environment.

Both platforms pass 209 deterministic CDF quantile controls. The maximum
observed relative small-tail probability residual is 1.75e-5, with a maximum
coarse/fine reference change of 2.68e-5. Each platform audits approximately
4.36 million independent global-envelope points, plus every final cell and
all interpolation-knot neighborhoods. There are no observed global or local
enclosure violations, and all 19 forced tail proposals are retained in each
case. These are numerical audits of the analytic bounds, not statistical
validation of a 1e-12-probability bin.

Native fixed-kernel timings include normalization, complete envelope construction
and sampling on every draw, in three interleaved blocks per method. Bounded
sampling costs approximately 16--50 microseconds per draw. Its remaining
construction cost makes it 1.20--3.15 times slower than the former unverified
rejection, while it is 3.79--31.58 times faster than CDF. No cross-event CDF
or envelope cache is used.

| Process / energy, eV | Bounded, microseconds | Bounded / old rejection | CDF / bounded |
| --- | ---: | ---: | ---: |
| ICS-1e+13 | 32.40 | 2.01 | 8.01 |
| ICS-1e+14 | 34.85 | 1.79 | 12.40 |
| ICS-1e+15 | 18.06 | 2.66 | 31.58 |
| ICS-1e+16 | 19.52 | 2.83 | 28.47 |
| ICS-1e+18 | 18.63 | 3.15 | 29.46 |
| ICS-1e+20 | 27.87 | 1.86 | 20.04 |
| PP-1e+14 | 50.22 | 1.20 | 3.79 |
| PP-1e+15 | 27.25 | 1.56 | 9.97 |
| PP-1e+16 | 16.03 | 2.37 | 22.20 |
| PP-1e+18 | 19.41 | 2.96 | 27.47 |
| PP-1e+20 | 25.25 | 1.70 | 21.04 |

Full-transport timings use three blocks of 100 photon primaries per method
and distance, on CMB at 2.7255 K, E0=1e15 eV, threshold=1e12 eV, source z=0.01,
no thinning, `max_rate_variation=0` and integration tolerance 1e-3. Entries
are medians; initialization and final disk output are excluded.

| Distance | Old rejection, s | Bounded, s | CDF, s | Bounded / old | Old S fraction | Bounded S fraction | CDF S fraction |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.1 Mpc | 7.708 | 8.274 | 10.585 | 1.073 | 1.60% | 3.27% | 24.48% |
| 1.0 Mpc | 25.297 | 26.065 | 29.494 | 1.030 | 0.74% | 1.51% | 12.29% |

These ratios include differences in realized trajectories/rate and continuous-loss
work, because changing the sampler changes uniform consumption and samples.
They are descriptive timing ensembles, not confidence bounds on overhead.
Transport profiling was frozen at de3b3c5. The final endpoint guards at 837d010
reproduced all 100-primary output-bin moments and S-call counts at both
distances; all 33 fixed-case raw streams also replayed byte-for-byte.
Final fixed-kernel timings and numerical/spectral controls were rerun at 837d010.
Linux emulation is used for compatibility/numerical checks; its timings are not
a forecast for a physical Linux amd64 server. Server pilots measure that CPU.

The final spectrum check contains 660000 interactions: 20000 per method and
case, with independent seeds for the three methods. All 220000 bounded events
use the analytic enclosure, with zero CDF fallbacks, violations, range narrowing
or argmax returns. An independent raw-file audit reproduces every bin and
checks kinematic support, particle types, energy conservation and all 161
saved hashes. Reanalysis and the independent raw audit under Python 3.10 /
NumPy 2.2.6 / SciPy 1.15.3 reproduce all precision states and significance
judgments from the native Python 3.12 analysis. Both frozen controls match
1000-event prefixes from the prior fa562d2 executable in all eleven cases.

At global alpha=1%, allocated as 0.5% to each comparison family, no significant
category comparisons or material spectral deviations are found. Among the
1362 bounded observables (including structural zeros), 690 meet the precision
target and 672 remain inconclusive. Direct bounded-versus-rejection differences
have 662 within-target and 700 inconclusive observables; bounded-versus-CDF
has 668 and 694. No claim of full 1%/10% equivalence follows from this short
sample, especially for rare tails. Development replays reuse seeds and are
not pooled to increase these event counts.

For example, at ICS E=1e20 eV the recoil-fraction tail below 1e-4 contains
4177/20000 bounded events. Its simultaneous relative-effect interval is
[-6.66%, +8.08%], inside the unchanged 10% target. Below 1e-6 there are only
21 events and the result remains inconclusive. Zero events below 1e-8 also
remain inconclusive rather than being interpreted as a physically absent tail.

The exact portable package completes its Ubuntu 22.04/Python 3.10 smoke in
285.45 seconds, with all three main streams (6000 events), both comparison
reports, auto-downloaded NR provenance and a retained result archive. Exit
code 3/inconclusive is expected for smoke. All 595 manifest entries and all
archive file contents were verified. The full 24-hour external-server experiment
has been prepared but has not been run locally.

All repeated timings, numerical controls, counters, source/binary hashes,
replay evidence and package identity are in
[bounded_rejection_20261020.json](data/bounded_rejection_20261020.json).
The package freezes 837d010; subsequent result-documentation commits do not
change its computational source or the archived experiment contract.

## Header readability follow-up, 2026-10-20

The constructor and nontrivial BoundedLogEnvelope methods now live in
src/lib/BoundedLogEnvelope.cpp; the header retains declarations and short inline
helpers. The analytic formulas, refinement limits, acceptance decisions and
random-number consumption are unchanged. Comments explain the Jacobian,
efficiency target, roundoff allowance and optional cross-section contracts.

The refactor passes all 155 CTest cases on native macOS arm64 and emulated
Ubuntu 22.04 Linux amd64, plus all 25 Python tests on native Python 3.12 and
Ubuntu Python 3.10. Replaying 20000 bounded interactions in each of the eleven
ICS/GammaPP cases reproduces all 220000 saved native raw events byte-for-byte.
This replay reuses the release seeds and supplies no additional independent
statistical evidence. The full server comparison remains pending; the server
pilot measures performance of the newly compiled source on the actual CPU.
