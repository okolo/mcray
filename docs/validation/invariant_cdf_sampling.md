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

Six new GoogleTest cases cover analytic power laws over 8.7--500
decades and amplitudes from 1e-200 to 1e200, both probability tails to 1e-12,
a separated 1e-7-probability mode, zero-density intervals, an extremely narrow
supplied feature, empty densities and fatal invalid inputs. All 144 short
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
The rejection control retains its GSL workspace between calls, as production
did. Its first 1000 raw events match the earlier control byte-for-byte in all
11 cases after this workspace-lifetime adjustment.
Native CDF costs range from approximately 0.23 to 0.70 ms per draw in this
family, approximately 5--107 times the former rejection cost. This first
implementation prioritizes explicit probability-error control; no CDF cache
across different particle energies, redshifts or process parameters is used.
CPU emulation and native server timings must not be treated as interchangeable.
The final emulated Linux measurements span approximately 0.93--30.0 ms per
CDF draw; these are not a forecast for a physical amd64 server. All 209 Linux
and 209 native quantile controls pass, with maximum observed relative
small-tail residual 1.75e-5 and coarse/fine reference change 2.68e-5.

Machine-readable measurements and environment records are in
[invariant_cdf_20261020.json](data/invariant_cdf_20261020.json); all 418
quantiles are in [the TSV table](data/invariant_cdf_quantiles_20261020.tsv).
The `linux` rows refer to amd64 emulation on the arm64 workstation.

## Fraction of transport and generation time spent selecting the invariant

The opt-in `em_sampling_probe` adapter now accepts `rejection` as a cascade
control and `MCRAY_PROFILE_S=1` to time invariant selection, including its
normalization. The control calls the unchanged production rejection function;
its body is byte-identical to revision `1cc28f3`. A thread-local MathUtils
retains the integration workspace. Production sources are unchanged by this
profiling addition. Kernel construction outside the sampler is not included
in the invariant timer. The loop timer excludes background initialization
and final disk output, and includes the profiling clock calls.

Native macOS arm64 measurements used three repetitions of 30 independent
primary histories per mode and distance, seeds 20261020--20261022. Cascades
start with 1e15 eV photons at z=0.01 on CMB at 2.7255 K, with GammaPP, hard
ICS and ICS continuous losses, a 1e12 eV transport/secondary threshold,
no thinning, and `max_rate_variation=0`. Background integration tolerance
is 1e-3. Entries below are medians over the three runs.

| Distance | Rejection loop, s | CDF loop, s | Rejection S fraction | CDF S fraction | Loop time ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0.1 Mpc | 2.218 | 3.349 | 1.71% | 27.10% | 1.51 |
| 1.0 Mpc | 7.574 | 8.706 | 0.76% | 13.76% | 1.15 |

These are small timing ensembles, not spectral acceptance tests. The sampler
changes random-number consumption and trajectories, so the loop ratios also
include differences in realized work. The earlier 107x result refers to a
particular fixed-energy invariant kernel, not a full transport calculation.
Most transport time in this configuration lies outside invariant selection.
The fraction depends on the transport setup and must be profiled on the
actual server workload before making server performance claims.

The separate `interactions` adapter measures SampleS plus SampleSecondaries,
vector construction and a checksum, without transport, histogram accumulation
or raw-event I/O. Across the original eleven physical cases, three blocks of
1000 events per mode (66000 events total), at integration tolerance 1e-4,
the median S fractions range from 71.65--94.68% for rejection and
98.95--99.67% for CDF. These fractions are not the fraction of an entire
server acceptance run, which also builds references, writes raw events and
accumulates/analyzes observables. Turning profiling on/off preserves the
1000-event secondary checksums exactly in four ICS/GammaPP controls.

Machine-readable timings, all repeated measurements, source/binary hashes
and configuration are in
[invariant_sampling_profile_20261020.json](data/invariant_sampling_profile_20261020.json).
The bounded runner preserves exact source, local source changes, adapter,
executable, seeds, commands, output/logs and SHA-256 hashes:

```bash
python3 tests/sampling/profile_cascade_sampling.py \
  --binary build/native/tests/em_sampling_probe --output /tmp/cascade-profile \
  --distances 0.1 1 --primaries 30 --repeats 3 --seconds 480
python3 tests/sampling/profile_cascade_sampling.py \
  --binary build/native/tests/em_sampling_probe --output /tmp/generation-profile \
  --scope interactions --primaries 1000 --repeats 3 --seconds 180
```

Each subprocess has an additional 60-second limit. A variable-rate pilot was
stopped because it exceeded the intended short profiling cost; it supplies
no performance conclusion. The full native runs above completed. Ubuntu
22.04/Python 3.10 under amd64 emulation also rebuilt the adapter and completed
all eleven generation cases (132 events), and all 19 Python tests passed on
both platforms. The emulated run is a compatibility smoke check; its timing
is not a physical server estimate.

## Bounded local spectrum check

The completed fixed-size native check uses production commit
`a960550f53ba7ce9c04f9b930183d294716418d0`, with 20000 interactions per sampler
per case, 440000 total, and a 1200-second limit. It completes in 129.46 seconds.
The subsequent commits adjust only control workspace lifetime and validation
documentation; the production source hashes in the measurement record match.

The independent 128/256/512-cell secondary-spectrum references, frozen bins,
paired-product inference and global 1% significance level are retained.
There are no significant category tests and no material deviations in either
sampler or the direct comparison. Candidate states are 688 within target and
674 inconclusive; the control has 690 within target and 672 inconclusive.
These region counts include structurally zero and overlapping views.

The ICS 1e20 recoil region `E_e/E_0<1e-4` contains 4211/20000 candidate events,
with simultaneous relative bias interval [-5.684%, +8.687%], within its 10%
target. Rarer regions remain inconclusive in this short ensemble. It does
not establish 1% precision in every core bin or 10% in every tail.

Every raw record passes the independent support/type/energy check, and every
frozen count table is reconstructed from the raw events. All local manifest
hashes match. The compressed local result archive has SHA-256
`2f583449a640910fa7d3c7f9ae05e9f162c1ef4bc1b394c65b0d34cb75f4dad1`.
No full electromagnetic cascade or nuclear acceptance ensemble was run.

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
