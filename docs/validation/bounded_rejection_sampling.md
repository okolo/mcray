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
silently increased after a proposal. One million unsuccessful attempts also
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
Confirmed measurements will be recorded below after the checks complete.
