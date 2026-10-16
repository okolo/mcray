# Conditional ICS and pair-production sampling correction

Implemented on branch `codex/fix-em-secondary-sampling`, forked from
`82df04adc32420fd68c12656e23873d7977d17fa` on `talys`, 2026-10-17.
The package manifest records the exact committed correction revision.

## Change and numerical checks

The [fixed-S diagnosis](ics_endpoint_diagnosis.md) identified an invalid Taylor
domain in the ICS integral and an artificial endpoint gap from the relative
photon-coordinate root tolerance. Both are corrected. The differential
physics laws, photon background and invariant-s sampling implementations
are unchanged.

`EMSecondarySampling.h` evaluates the analytic ICS integral with `log1p` and
a stable logarithmic remainder. For `q=(y-yMin)/yMin < 0.01`, the logarithm and its
remainder use 18-term convergent series. The small-argument logarithm also
avoids measured loss of accuracy in x87 log1pl under Linux amd64 emulation. The
invalid absolute-width Taylor switch is removed. Both branches use the
same analytic integral; this is a numerical correction, not a new spectrum.

ICS now inverts the recoil CDF in log(y), evaluates that integral directly
from y, and creates the electron with E0*exp(log(y)) and the photon with
E0*(-expm1(log(y))). Thus neither secondary's small energy is reconstructed
by subtracting nearly equal large energies.

GammaPP uses a stable lower bound `2/[S*(1+sqrt(1-4/S))]`, a cancellation-safe
analytic integral of the same pair differential law, and inversion in log(r).
Its primitive uses `log1p` differences and factored rational differences,
avoiding subtraction of primitives at the two endpoints.

Both inversions target absolute CDF residual
`1e-12*min(u,1-u)`, with termination at adjacent representable logarithmic
coordinates when floating-point resolution is exhausted. The shared
MathUtils root solver is unchanged. Floating-point saturation prevents an
unconditional relative-residual guarantee for arbitrarily small probabilities.

Five new GoogleTest cases check the analytic integrals against direct positive
differential integration, production quantiles and conservation/support,
threshold and extreme-energy cases, and inverse probabilities down to 1e-12
in both tails. Narrow reference integration intervals are mapped to [0,1]
before calling GSL to avoid losing their width in double endpoint conversion.

Short results:

- macOS arm64: all 137 CTest entries pass, including the five new tests.
- The Python statistics/runner suite passes all 15 tests.
- Fixed-S ICS diagnostic: maximum tested production CDF residual 4.71e-13;
  independent quadrature orders 16 and 32 differ by at most 1.34e-15.
- A fixed-seed 10000-event CMB ICS control at E0=1e20 eV contains **2109**
  recoil electrons with y<1e-4 (independent server reference approximately
  2084 expected). The previous full server run had none.
- A fixed-seed 100000-event CMB PP control at E0=1e14 eV contains **37**
  lower-energy leptons with 0.49995<r<0.5 (reference approximately 37.4
  expected). The previous current-sampler run had none in this region.

These controls demonstrate restored endpoint support, not acceptance of all
spectra at 1%/10%. Full ensembles remain necessary. Local logs, events,
count tables and numerical diagnostics are retained outside Git in the
adopted workspace at `validation/em-secondary-fix-20261017/`.

## Frozen server rerun without publishing the branch

Create the self-contained package after committing the correction:

```bash
python3 tests/sampling/package_secondary_rerun.py --output /path/to/new/package
```

Copy the resulting tar.gz to the native Linux amd64 server, extract it,
enter its directory, and run:

```bash
bash sampling/run_fixed_server.sh \
  --data-dir /home3/ok/mcray-test-data --hours 24 --jobs 4 --raw-gib 20
```

The wrapper verifies SHA-256 hashes, uses the companion Git bundle and selects
the pinned Python 3.10 requirements on Ubuntu 22.04. Python 3.12 is not needed.
It installs dependencies in a package-local virtual environment. Environment
setup precedes the bounded experiment; the Python runner enforces the same
24-hour total including source fetch/build, references, sampling, inference
and archiving. NR is downloaded at its pinned revision unless `--nr` overrides
it. No push or private-workspace access is needed.

All 11 cases, 26 nonuniform bins, cumulative and endpoint tails, seeds,
reference convergence checks, the 1% core / 10% tail targets and family-wise
significance 1% are retained. Statistics, raw-data limits, event records,
freezing and archive retention follow the unchanged
[server protocol](../../tests/sampling/SECONDARY_SERVER.md).

Both `old` and `current` modes use the corrected conditional secondary sampler;
those names still select only the historical NR versus current invariant-s
sampler. They do not mean pre-fix versus post-fix conditional sampling.
The original server archive remains the pre-fix evidence. The aggregate
report can still flag deviations from historical NR or the separate rejection
envelope defect; assess the corrected `current`/oracle regions individually.
Rare regions without sufficient precision remain inconclusive.
