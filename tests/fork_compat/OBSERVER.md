# Historical observer and ballistic stopping diagnostic

This is a deterministic prerequisite to the proposed source/transport/observer
ensemble comparison. It does not sample sources, add magnetic fields or
interactions, compare source luminosity, or assign physical/statistical
acceptance. The [transport diagnostic](TRANSPORT.md) remains independent.

## Reference and reproduction

Use the immutable historical snapshot pinned in AGENTS.md, or a read-only
checkout matching `observer_reference.json`. That manifest pins seventeen
files: the previous transport inputs plus Observer.cpp/h, ParticleStack.cpp/h
and the audited PropagationEngine.cpp. The published snapshot contains all
seventeen with matching original bytes.

From a compatible configured modern Unix Makefiles build with grid tests:

```sh
mkdir -p build/fork-observer/tmp
TMPDIR="$PWD/build/fork-observer/tmp" PYTHONDONTWRITEBYTECODE=1 \
  python3 tests/fork_compat/observer.py \
  --fork-root /read-only/path/to/CRbeam_v0 \
  --build-dir build/Release --output build/fork-observer/run
```

Python >=3.8 standard library and the configured build's C++/GSL/Boost/Fortran
toolchain are needed. `--command-prefix` accepts the same offline container
argv prefix as the existing adapters. Compilation and probe execution use it;
the probe has a 120 s timeout. Use the instrumented build to check extracted
historical functions and modern code with ASan/UBSan too.

The adapter copies the original Observer constructor, membership and step
functions verbatim. It also extracts Particle::PropagateFreely and the actual
Result continuation/time-left functions, supplying only their storage shell
with filter disabled. Existing hash-checked particle and shared time/unit
definitions come from the transport adapter. The old ballistic driver is the
interaction-free, no-abort, no-deflector reduction of the audited constant-rate
engine loop: original continuation, patched remaining time and original free
motion. A 100000-step budget rejects a stalled/expensive historical case.
This is not a full historical engine build or a variable-rate reference run.

The modern side runs the actual Result/PropagationEngine with a capture output
in both constant- and variable-rate modes, and rejects reported failed histories.
Both receive the same initialized photon positions, directions, energy and
time, with H0=70 km/s/Mpc and OmegaLambda=0.7. Normal ballistic time starts at
z=0.01. No source sampler or source-axis conversion is exercised here.

## Frozen cases and numerical interpretation

`inputs.json` and the case hash are written before compilation/evaluation.
There are 210 cases:

- 168 normal-scale point/step queries: radii 0.01/0.1 Mpc, origin/translated
  centers, strict inside/outside modes, seven positions (center, interior,
  surface, exterior, adjacent-representable points and a diagonal near the
  surface), and three positive remaining-time fractions (10,0.01,1e-9 of R).
- Six range controls use radius 1e165 Mpc, whose double-precision squared
  internal value overflows. These are arithmetic tests, not physical runs.
- 36 ballistic cases: two radii, two centers, and nine scenarios covering
  birth hits, entry, exit, misses, surface directions, a short end-time remainder
  and no observer. Each is run through both modern rate paths.

Point queries compare raw kernels, even where continuation would already stop
at birth. Therefore a point-kernel step difference is not automatically a
realized trajectory difference. The diagonal and adjacent-representable points
probe floating-point classification; their represented coordinates need not
lie exactly on the mathematical sphere. Range controls preserve the fork's
overflow behavior rather than concealing it.

`observer.jsonl` retains every result. `report.json` lists all membership,
step-cap, ballistic-hit and rate-path differences without a verdict. Final-time
differences are reported in Mpc and in ulps of the represented end time; the
`past_end` flag is an exact comparison, not an arbitrary tolerance. A one-ulp
excess must be distinguished from the old minimum step substantially exceeding
remaining time. The original limit and both modern engine paths remain visible.

`manifest.json` records all reference/adapter/modern source and linked-core
hashes, revision/dirty state and execution prefix. Inputs, compiler command,
compiler version and compilation/execution logs are retained under `--output`.
Failed reruns remove stale success reports. All outputs must remain outside
the reference directory and inside any workspace-required writable root.

Self-contained guards run through `python3 tests/test_fork_transport.py -v`
and the existing Release CI step, without an external checkout. The guards
cover case coverage, reference hash rejection and preservation of differences.
