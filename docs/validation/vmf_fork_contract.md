# Restored historical vMF contract

On 2026-10-31 the user required maximum fidelity to the fork's algorithms,
weights and recognizable code structure, except for demonstrable defects.
This supersedes the normalized conditional-vMF decision implemented in PR #12.
The [current usage contract](../usage/vmf_source.md) defines original weights,
parameter ranges, activation/zero branches and numerical exceptions.
The [old normalized-model report](vmf_source.md) remains a frozen historical
record and is not evidence for this restored contract.

## Implementation scope

Restore Jet.h/Jet.cpp and the original `vMF_weight` function name/angle-based
entry point, with a small normalized-axis overload for the current CLI.
SourceGeometry follows the original sequence: uniform-cone proposal, then
multiplication of the incoming primary weight by the full-sphere density.
Remove the direct vMF inverse CDF, conditional normalization mass/cache and
small-log series. No alternative normalized mode is added. CLI/YAML names,
current particle model/RNG and ordinary pencil/cone behavior are retained.

Finite signed kappa is accepted, including the old abs(kappa)<1e-14 unit
multiplier branch and its discontinuity. The vMF angular activation threshold
is the original opening>1e-14 radians. Inactive opening leaves weight/RNG
unchanged. The current explicit-axis orientation remains an approved interface
extension over the fork's ignored zero-opening rotation.

Documented arithmetic fixes use expm1 for cancellation, cosine reflection for
negative kappa, a log fallback for representable density tails, a bounded
roundoff correction at unit-vector poles and fatal checks for nonfinite or
unrepresentable results. They preserve the original formula rather than
renormalizing its source strength. No field, observer, thinning, energy loss,
background, Gaussian source or nuclear transport algorithm is changed.

## Software and direct-reference checks

The self-contained C++ tests check:

- Same modern seeded directions and next RNG value as the uniform cone;
  weight is multiplied only after drawing the direction. Zero/small-kappa
  and inactive-angle paths retain their documented behavior.
- Independent 100-digit uniform-proposal point checks and full-sphere density
  checks, including signed kappa, the zero-branch discontinuity, extreme
  concentrations and a representable tail lost by the original exponential.
- Preservation of identity, nuclear metadata, energy, ancestry and times;
  weight changes by exactly the prescribed multiplier.
- The existing 40000-draw/20-bin/six-sigma software controls now test the
  uniform proposal CDF u/U, not the superseded direct-vMF CDF. Sample size,
  seed and grid remain explicit; this is not physical acceptance.
- CLI/YAML/snapshot/replay equivalence, signed parameters, weighted neutral
  ballistic output, invalid input and weight-overflow failure.

The optional hash-pinned `tests/fork_compat/vmf.py` compiles the actual original
Jet.cpp and compares weights at identical directions with modern Jet.cpp.
It freezes 87 points before execution: 70 ordinary points (two axes, five
unit directions and seven signed/zero parameter values), twelve cancellation
controls near the old threshold and five range controls. Nonfinite original
values remain visible as such. Numerical controls are not pooled into an
ordinary maximum or a claim that defects were reproduced.

```sh
mkdir -p build/checks/vmf-tmp
TMPDIR="$PWD/build/checks/vmf-tmp" PYTHONDONTWRITEBYTECODE=1 \
  python3 tests/fork_compat/vmf.py --fork-root /path/to/pinned/reference \
  --build-dir build/Release --output build/checks/vmf-points
```

Use the immutable public reference pinned in AGENTS.md. Python standard
library plus the configured C++/GSL/Boost/Fortran toolchain is sufficient.
Output must be new/empty and outside the reference. `--command-prefix` supports
the existing frozen offline container. Inputs, reference/production/adapter/
core hashes, compiler command/version, raw points and execution logs are
retained. Normal CI runs the self-contained Python/design guards and C++ tests
without the historical checkout; it does not require external reference data.
No new source/observer/cascade Monte Carlo acceptance ensemble is run here.
The earlier 48-million-history pencil/cone result remains scoped to its
original frozen source revision and has not been rerun.

## Results on 2026-10-31

Production restoration is committed as `8f9f681`; final point probes ran from
clean `7d2c924`. [Machine-readable evidence](vmf_fork_contract.json) records
source/reference/input/adapter/core hashes, raw/log artifact hashes and the
complete numerical controls. All three environments completed the same 87
point design. On the 70 ordinary points the maximum relative difference from
the actual historical function is **3.59268e-15** in each environment.

Numerical controls retain intentional departures explicitly:

- Near abs(kappa)=1e-14, subtraction in the fork denominator loses up to
  0.07993% relative precision; expm1 agrees with the independent high-precision
  formula. The original unit branch below the threshold is still exact.
- In the representable density-tail control, the fork returns zero while the
  corrected result is approximately 3.03e-321. For negative concentration the
  original can instead become nonfinite.
- At the opposite pole for kappa=-1000 or -1e308 the original formula is
  nonfinite, while the equivalent stable expression gives a finite peak.
  The positive kappa=1e308 pole agrees with the fork exactly in this design.

These arithmetic controls are not passed off as matching original values or
as reasons to change source normalization. All ordinary point results remain
separate from them.

| Environment | Complete CTest suite | Final targeted vMF/configuration suite |
| --- | --- | --- |
| macOS arm64 Release, Clang 22.1.5 | 250/250, 36.19 s | 39/39, 0.72 s |
| Emulated Ubuntu 22.04 Linux amd64 Release, GCC 11.4 | 250/250, 221.38 s | 39/39, 15.31 s |
| Same Linux build with ASan/UBSan and leak detection | 250/250, 343.89 s | 39/39, 20.60 s |

Seventeen self-contained Python guards passed and use the existing Release
CI entry point. The frozen offline image and exact compiler/link recipes are
recorded in the evidence. No sanitizer diagnostics or suppressions occurred;
system libraries are not rebuilt with instrumentation. Linux execution was
under emulation, not a native Linux performance benchmark.

An additional tail assertion initially converted a subnormal Boost reference
to double zero. The final assertion compares in least-subnormal units, keeping
the same four-ulp bound. The initial failed assertion log is retained. The
final targeted suites above cover this strengthened independent check; no
production change, source resampling or acceptance-band widening was needed.

Raw points, compiled original fragments, executables and logs remain under
ignored `build/checks/vmf-restore-points-final-*`. Complete build/test logs and
the final independent-reference checks are retained separately under
`build/checks/vmf-restore-linux-*` and `build/checks/pilot-tmp`. The input hashes
and final source revisions make the point comparison reproducible without the
private workspace. These results do not establish full weighted cascade or
arbitrary-parameter statistical equivalence with the fork.
