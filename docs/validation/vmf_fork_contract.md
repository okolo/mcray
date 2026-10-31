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
