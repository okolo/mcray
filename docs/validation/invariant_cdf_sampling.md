# Numerical CDF inversion for background invariant sampling

The candidate branch `codex/inverse-cdf-sampling` starts at the corrected
ICS/GammaPP revision `1cc28f381e961f3d4947d08f0149ac72d63c3635`.
The cross sections, conditional secondary spectra and photon-background
interpolation are unchanged. Background invariant sampling now calls
`MathUtils::SampleLogDistributionCDF` instead of rejection sampling.

## Implemented numerical contract

`NumericalLogCDF` integrates `x*f(x)` in a logarithmic coordinate over the
entire requested physical support, clipped only to `Function::Xmin/Xmax`.
Known interpolation knots and physical features can be supplied through
`Function::IntegrationBreakpoints()`. Table functions expose their physical
knots, and background kernels map the background knots into invariant/energy
coordinates, including both kinematic arguments and both redshift tables.
This metadata does not modify the underlying function or interpolation law.

Normalization uses adaptive positive Gauss-Kronrod 7/15 integration. Its
partition and cell integrals are retained for the quantile's initial guess
and partial integration. The guess is never returned without a probability
check. Partial integrals preserve known zero-interval boundaries and refine
the selected tail as needed; they are not obtained by linearly interpolating
a sampled CDF table.

Each sample consumes one open-uniform variate. For `u<=0.5`, inversion uses
the directly integrated lower tail; otherwise it integrates from the upper
endpoint with target probability `1-u`. The coordinate is anchored at that
endpoint, avoiding subtraction of almost equal CDF values and loss of log
coordinate resolution near a distant endpoint. Safeguarded Newton iteration
retains the full bracket and uses bisection whenever its Newton step is unsafe.

The small-tail relative probability tolerance is `min(requested,1e-5)`.
Quadrature receives 1/32 of this tolerance; the stopping condition includes
the estimated partial-integral error and normalization uncertainty. Probes
far below the target use a target-scaled absolute quadrature budget, avoiding
unnecessary refinement of negligible intermediate probabilities. The returned
normalization is also evaluated more tightly than the previous loose rate
quadrature; exact trajectory/rate bit identity is not promised.

Negative/nonfinite densities, quadrature failure, nonconvergence or exhausted
coordinate resolution throw a non-retryable `Utils::Exception`. No support
narrowing, maximum substitution or rejection fallback is used. A zero
normalization returns the established empty-interaction result, not a sample.
The work is bounded by 4096 quadrature cells per integral and 128 inverse
iterations. Unsupported sub-resolution tolerances fail explicitly.

Gauss-Kronrod errors remain numerical estimates, not rigorous enclosures for
arbitrary black-box functions. Unannounced narrow spikes can defeat numerical
quadrature; known features must be supplied as integration breakpoints.
Returning a binary64 energy imposes an unavoidable probability resolution
of order `f(x)*ulp(x)/integral` at a nonzero-density endpoint. Analytic tests
include that representation limit explicitly.

## Scope and regression checks

Both continuous and monochromatic background invariant paths use the new
method, for `s` and for background-photon energy in the particle rest frame.
The previous `SampleLogDistribution` remains for the separate nested proton
pair secondary PDF. That PDF can itself be evaluated at percent-level
quadrature accuracy near threshold; forcing the tighter CDF contract on it
caused nonconvergence in the existing threshold regression. Its conditional
model and sampler are outside this EM invariant correction. The separate
linear-coordinate rejection sampler is also unchanged.

Five new GoogleTest cases cover analytic power laws over 8.7--34.7
decades and amplitudes from 1e-200 to 1e200, both probability tails to 1e-12,
a separated 1e-7-probability mode, zero-density intervals, an extremely narrow
supplied feature, empty densities and fatal invalid inputs. All 143 short
CTest entries pass on macOS arm64 and Ubuntu 22.04 Linux amd64 under emulation.

The portable Python suite includes tests for missing/nonfinite deterministic
records, a failing quantile control, unresolved reference precision and
complete construction-cost timings; all 19 tests pass locally.

## Deterministic physical-kernel controls

For each of the original six ICS and five GammaPP cases, the diagnostic
checks 19 frozen quantiles spanning both tails down to 1e-12. The reference
uses independent composite fixed Gauss-Legendre integration, with separately
accumulated prefixes/suffixes at 8192/16384 intervals and orders 16/24.
It does not use the production adaptive quadrature or inverse solver.

Native macOS controls find a maximum relative small-tail probability residual
of 1.75e-5 and a maximum coarse/fine reference change of 2.68e-5. These are
observed numerical diagnostics, not certified error bounds. The diagnostic
acceptance control is 1e-3, with a 1e-4 reference-change guard; this does not
replace the spectral Monte Carlo 1%/10%, global alpha=1% criteria.

Timings include normalization construction and inversion on every draw,
using three interleaved 128-draw blocks for CDF, frozen rejection and NR.
Native CDF costs range from approximately 0.23 to 0.68 ms per draw in this
family, approximately 5--101 times the former rejection cost. This first
implementation prioritizes explicit probability-error control; no CDF cache
across different particle energies, redshifts or process parameters is used.
CPU emulation and native server timings must not be treated as interchangeable.

## Server experiment

Use [INVARIANT_CDF_SERVER.md](../../tests/sampling/INVARIANT_CDF_SERVER.md).
The default comparison is candidate CDF versus frozen rejection at `1cc28f3`,
with both compared against the unchanged independent secondary-spectrum
oracle. `--baseline nr` selects the historical control instead. Control labels
and commits are frozen in each archive; plots read the recorded labels.

The 11 cases, 1362 bin/tail regions, 12258 category hypotheses, exact paired
product inference, seeds and 1%/10% relative accuracy targets at global
significance 1% are retained. Setup, deterministic controls, pilot forecasting,
sampling and final archive all share the maximum 24-hour budget. More costly
sampling can leave more rare regions statistically inconclusive.

Full-spectrum server acceptance is not inferred from the deterministic
controls or short local smoke runs. All reference-convergence records, raw
events, source versions, seeds, timing blocks and failures are retained.
