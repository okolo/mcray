# ICS endpoint error: fixed-invariant diagnosis

This document records the pre-fix diagnosis. The subsequent
[correction and repeat-server workflow](em_sampling_fix.md) are documented
separately. The current diagnostic executable reports corrected production
quantiles alongside the retained legacy primitive variants; the numerical
tables below describe the archived pre-fix executable.

Diagnosed on 2026-10-17. Two independent defects in conditional secondary
sampling explain the missing low-energy recoil electrons and the high-energy
photon cutoff. Both defects are present in `d05723f449c1eecee4f00af392014df4c898ebe9`
and in the server physics revision `cb5b2d9210846c6ba6b46ff2c4bb80ffb47b9751`.
The relevant production implementations are unchanged at the local diagnosis
revision `82df04adc32420fd68c12656e23873d7977d17fa`.

This diagnosis supplements the [server results](em_sampling_server_results.md)
and [ICS spectra and primary formulas](ics_secondary_spectra.md). No production
physics or sampling code was modified in this investigation.

## 1. Invalid domain for the integral's Taylor expansion

In [ICS::PartialSigma](../../src/lib/ICS.cpp), define
`S=s/m_e^2`, `yMin=1/S`, `yMax=1-r`, and `delta=yMax-yMin`.
The exact analytic expression is selected when `delta > 1e-3`; otherwise the
code evaluates a fifth-order Taylor polynomial about `yMin`.

The polynomial contains powers of `delta/yMin`. In particular, the logarithm
in the exact integral expands as `log(1+delta/yMin)`. Its Taylor series requires
`abs(delta/yMin)<1`, and a fifth-order approximation needs a considerably
smaller ratio for accurate probabilities. An absolute bound on `delta` does
not enforce this condition. At large S, the current branch evaluates the
polynomial far outside its convergence domain.

For S=1e6 and photon threshold rMin=1e-8, normalize the integral by its value
at rMin. It should equal the recoil-electron CDF, P(Y<=y).

| y | delta/yMin | Production integral ratio | Direct positive-density integral |
| --- | ---: | ---: | ---: |
| 2e-6 | 1 | 0.05471944 | 0.04841942 |
| 1e-5 | 9 | 725.17535 | 0.16084561 |
| 1e-4 | 99 | 131208745.47 | 0.32169100 |
| 5e-4 | 499 | 431169366911.25 | 0.43411710 |
| 1.1e-3 | 1099 | 0.48919423 | 0.48919423 |

The ratio exceeds one by many orders of magnitude, then drops when the exact
branch resumes. Consequently the function passed to the inverse sampler is
not a valid monotone CDF. Increasing root-solver accuracy cannot repair it.
The independent reference integrates equation (12) of
[ELMAG, arXiv:1106.5508v2](https://arxiv.org/pdf/1106.5508v2), without this
primitive or its Taylor series.

## 2. Inverse sampling stops at insufficient endpoint resolution

[ICS::SampleSecondaryGammaFracE](../../src/lib/ICS.cpp) solves for the photon
fraction x using bisection and relative coordinate tolerance 1e-3.
[MathUtils::SolveEquation](../../src/lib/MathUtils.cpp) stops with
`gsl_root_test_interval(x_lo,x_hi,0,relError)`.
For x near one this permits an absolute interval of approximately 1e-3.
The relative error of the recoil fraction y=1-x is amplified as
`abs(delta_x)/y`; it is not controlled by this stopping rule.

The nearest returned midpoint to the upper endpoint is one part in 2048 of
the initial interval. Thus the artificial recoil floor is approximately
`1/S + (1 - 1/S - rMin)/2048`, approaching 4.8828e-4 for large S.
This is exactly the dyadic endpoint pattern found in the saved server events.
The invalid Taylor branch changes the probability assigned to those positions;
the loose solver independently prevents resolving the spectrum below the floor.

[GammaPP::SampleSecondaryFracE](../../src/lib/GammaPP.cpp) uses the same 1e-3
coordinate tolerance with upper bound 0.5. It similarly misses a narrow
equal-energy splitting region. The ICS Taylor defect is separate and does
not occur in that PP primitive.

## Deterministic reproduction and independent controls

The [diagnostic executable](../../tests/sampling/ics_endpoint_diagnostic.cpp)
links unmodified production code. It calls public `ICS::SampleSecondaries`
at fixed S, with seeds 1 through 12 and S=10,100,1000,1e4,1e6,1e8: 72
quantiles, at parent energy 1e20 eV and photon threshold 1e12 eV. Neither the
NR sampler, rejection sampler nor background convolution is used. Its replica
of the original inversion agrees with all 72 production energy outputs to
floating-point rounding (maximum native difference 1.12e-16 in energy fraction).
The replica uses the existing production
`PartialSigmaSeries` and `PartialSigmaExact` helper functions and the same
branch condition, normalization and root solver.

The reference CDF integrates the positive differential density in log(y).
Gauss orders 16 and 32 agree within 1.34e-15 in absolute CDF probability
over the checked primitive grid. Four variants isolate the two defects:
original production, tolerance alone reduced to 1e-12, exact primitive alone,
and both diagnostic changes together. These substitutions exist only in the
diagnostic, not the production implementation.

Example: S=1e6, seed=2, uniform deviate u=0.480968228982789:

| Variant | Returned recoil fraction y | Reference CDF at returned y |
| --- | ---: | ---: |
| Production | 0.000489280756836 | 0.432603241121 |
| Tighter tolerance only | 0.0000031423930964 | 0.079982237903 |
| Exact primitive only | 0.001465842270508 | 0.509250805310 |
| Both diagnostic changes | 0.000977800911080 | 0.480968228959 |
| Independent inverse CDF | 0.000977800911411 | u |

Tightening the tolerance alone can make the physical answer substantially
worse because it more accurately inverts the invalid polynomial. Replacing
the polynomial alone leaves the endpoint quantization. Both must be addressed.
Even the combined diagnostic's fixed photon-coordinate tolerance is not a
general endpoint solution: at S=1e8 its maximum tested CDF residual is about
6.46e-7. A production fix should use a stable monotone integral and explicitly
control the endpoint variable and probability residual across the full domain.
Do not globally change the semantics of the shared root solver for this fix.

Build and run without NR:

```bash
cmake -S src/app/crbeam -B build/Release -DBUILD_TESTING=ON
cmake --build build/Release --target ics_endpoint_diagnostic -j4
build/Release/tests/ics_endpoint_diagnostic /path/to/existing/results/ics-fixed-s
```

The executable writes `PREFIX-primitives.tsv` and `PREFIX-quantiles.tsv` and
prints replica and quadrature control errors. Exit zero means that these
diagnostic controls passed, not that production sampling is accurate.
Native macOS arm64 and Linux amd64 reproduction tables and source/binary
hashes are retained with the external server analysis under
`validation/em-secondary-analysis-20261016/` in the adopted workspace.
These files are not dependencies of public builds or workflows.

After a production fix, rerun the frozen full spectral test, including tails
and the unchanged 1%/10% criteria at family significance 1%. The short
fixed-S diagnosis is not a replacement for that acceptance campaign or for
the separate rejection-envelope investigation.
