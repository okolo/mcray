# Frozen source/observer ensemble protocol v1

Approved by the user on 2026-10-31, including the observable family,
confidence method, counts and seed reservations below. The user additionally
authorized local execution with a hard four-hour limit, overriding the usual
30-minute workstation limit for this specific comparison. All other scope,
failure and preservation rules remain in force. Reproduction on native Linux
amd64 remains the primary-platform follow-up; native macOS results must be
identified as such.

### Model and sample reservation

Retain the three pilot geometries without extending the physical scope:
photon, 1e12 eV, initial weight 1, z=0.01, H0=70/OmegaLambda=0.7, origin
source, inside sphere R=0.01 Mpc at distance d=0.1 Mpc along the declared
axis, a 0.2 Mpc horizon, no fields/interactions, constant-rate propagation.
Use the same represented observer center for both implementations. Preserve
strict endpoint membership and the original source/observer algorithms,
including recorded historical numerical differences.

| Case | Half-angle, degrees | Axis (theta,phi), radians | Fork seed | Modern seed | Histories per implementation |
| --- | --- | --- | --- | --- | --- |
| pencil_z | 0 | (0,0) | 2004101 | 2004102 | 8000000 |
| cone_z | 20 | (0,0) | 2004103 | 2004104 | 8000000 |
| cone_rotated | 20 | (0.4,0.7) | 2004105 | 2004106 | 8000000 |

There are 48000000 complete primary histories in total. These source seeds
are reserved for this ensemble; their consumption is recorded in each run manifest. Exclude all pilot seeds 2003101 through 2003106.
Maintain one uninterrupted source RNG stream per implementation/case; batching
must not reseed it. Engine seed 17 is retained: there are no interaction draws
in this model. The historical reference remains snapshot e1fc71f and its
seventeen pinned input files. Record full source revisions and hashes before
execution; resolve the requested modern Git revision at server startup.

### Frozen observables: 33 differences plus 6 hit fractions

Use each complete primary once, including misses. In each of the three cases
compare the modern mean minus the historical mean for these eleven quantities:

1. Hit indicator H, bounded in [0,1].
2. Four unconditioned arrival indicators H * I(tau <= t), with
   t=0.092, 0.094, 0.096 and 0.098 Mpc of c*time, each bounded in [0,1].
   Tau is the actual final time minus initial time, including represented
   absolute-time subtraction. Misses contribute zero. These are not CDFs
   conditioned on arrival. Thresholds are chosen from the fixed d/R geometry
   before final data; they avoid the pencil atom at d-R=0.09 Mpc.
3. Three final coordinate components X_j / L for all histories, including
   misses, with L=0.25 Mpc and bounds [-1,1]. Report intervals also in Mpc.
   This is a statistical scaling only, not an output-coordinate change or
   the legacy per-file T scale. The geometry/horizon and historical minimum
   step leave margin below L; a range violation invalidates the run.
4. Three emitted direction components Pdir_j, bounded in [-1,1], without
   rotating or renormalizing either implementation's sampled direction.

This gives 3*11=33 difference intervals. Add six individual hit-fraction
intervals, one per implementation/case, for a total simultaneous family
M=39. This makes the 0.001 precision target explicit for both the differences
and individual hit fractions. Dependence among observables from the same
history is allowed; do not treat them as additional independent histories.
No conditional moments, adaptive bins, angular distribution scan or other
post-hoc statistic enters this family. Additional descriptive diagnostics
must be labelled exploratory and supply no acceptance claim.

### Confidence construction and fixed counts

Use a two-sided Hoeffding bound for each mean/difference with equal error
allocation alpha/M=0.05/39, followed by the union bound. This conservative
finite-sample construction avoids normal approximations or pilot variance
estimates, including at zero/all-hit boundaries. Its assumptions are the
usual Monte Carlo independent-primary/independent-stream model and the fixed
observable bounds. Correlation across different observables does not affect
the union bound. RNG quality is an assumption, not established by this test.

For a shared observable range [a,b], write w=b-a and let N be the number of
histories in each implementation. Independence of the two samples gives:

    Pr(|(mean_modern - mean_fork) - true_difference| >= h)
        <= 2 exp(-N h^2 / w^2).

Therefore choose the difference radius:

    h_difference = w * sqrt(log(2*M/alpha) / N).

For an individual hit fraction choose:

    h_individual = sqrt(log(2*M/alpha) / (2*N)).

Use closed, untruncated intervals centered on the observed mean or difference.
In particular, an all-hit pencil does not receive a zero-width interval.
Do not clip intervals to the physical parameter range before checking their
precision. For N=8000000 and M=39, the deterministic radii are:

| Quantity | Half-width |
| --- | --- |
| Hit-fraction or arrival-indicator difference | 0.000958673635 |
| Individual hit fraction | 0.000677884628 |
| Direction-component mean difference | 0.001917347269 |
| Coordinate mean difference | 0.000479336817 Mpc (0.479337 kpc) |

The worst-case minimum N for a difference half-width at most 0.001 is
ceil(log(1560)/0.001^2)=7352442; round up to 8000000 before data collection.
The moment widths are stated separately; the user-approved 0.001 target
applies to hit fractions, not to every dimensional moment. If narrower moment
precision or a practical equivalence margin is desired, that needs a separate
decision before the run.

This formula is derived for this two-sample design from Hoeffding's bounded
independent-sum inequality; the factor of two differs from a one-sample mean.
Primary source: W. Hoeffding (1963), Probability Inequalities for Sums of
Bounded Random Variables, JASA 58, 13-30,
https://doi.org/10.1080/01621459.1963.10500830 (Theorem 2).
No distribution-free method can establish complete physical equivalence from
this finite family of observables alone.

### Decision wording and numerical diagnostics

If all data/precision checks pass and all 33 difference intervals contain
zero, report compatibility with zero difference for the frozen observables
at the reported resolution. This is non-detection of a difference, not proof
of exact equality or a practical equivalence margin. If any difference
interval excludes zero, report a detected discrepancy and retain its case,
observable, interval and raw provenance for investigation. Do not widen the
interval, change the family, increase N after inspecting data, or rerun the
same reservation to seek a different outcome.

Weights, identity, energy, direction/coordinate bounds, finite times and
positions, contiguous primary IDs/counts, failure count and source/input hashes
are validity checks. Initial/final weights must remain 1 for this agreed
ordinary-cone model; preserve an unexpected original value in diagnostics
and invalidate the statistical report instead of renormalizing it. A failed
history, missing/duplicate record, nonfinite value or range violation makes
the comparison invalid/inconclusive. Never count it as a miss or replace it.
A partial run produces no successful comparison report.

Keep historical step counts, signed final-time errors and modern failure
counts as separate diagnostics. The known historical step-floor overshoot
must remain visible, not be reclassified automatically as a sampled miss or
silently repaired. Counts/intervals refer to the actual endpoint membership
of each algorithm. No timing-overshoot acceptance tolerance is introduced.

### Execution and resource limits

Use an external native Linux amd64 server for primary-platform reproduction. The existing pilot under amd64
emulation extrapolates to about 3.1 hours of sampling/propagation alone for
this count; it is not a native-server benchmark. The user explicitly authorized a local run of this protocol with a hard
14400-second collection limit; that exception does not authorize other
long workstation calculations. Normal C++/Python regression tests and a fixture-based
adapter check precede it; do not draw the final seed reservations as a smoke
test. Reference/hash failure or a failed required check prevents collection.

The public runner must record the exact resolved code/reference revision,
compiler, dependencies, numerical controls, seeds, counts and hashes; preserve
exit codes/logs and compress the final or partial artifacts into the configured
test-data directory. It must not require this private repository. An output
directory/reservation must not be silently overwritten or resumed/reseeded.
A failed attempt remains spent and requires a separately reviewed recovery.

Use streaming or bounded blocks for reduction and compressed per-history
retention. The pilot-size extrapolation is about 17 GiB of uncompressed JSONL
for 48 million histories; compression and throughput are unmeasured. Reserve
sufficient disk and retain all raw histories, manifest and summary hashes.
Keep counters/aggregates per case/implementation; do not load the full ensemble
or millions of primary IDs into Python sets. Validate IDs in monotone order.
Block boundaries do not change statistical sample count or reset RNG state.
Freeze this public protocol and implementation before the
first final history. This protocol authorizes no field coupling, translations,
outside-observer ensemble or interacting cascades.

## Runner and reproduction

Use Python >=3.8 and the native compiler/link recipe from a configured build
with `grid_magnetic_field_tests`, GSL, Boost and Fortran runtime dependencies.
Use the immutable public snapshot retrieval in AGENTS.md; the reference is
read-only and is not needed by ordinary CI. From the public repository root:

```sh
mkdir -p build/checks/ensemble-tmp
TMPDIR="$PWD/build/checks/ensemble-tmp" PYTHONDONTWRITEBYTECODE=1 \
  python3 tests/fork_compat/source_observer_ensemble.py \
  --fork-root /path/to/pinned/reference --build-dir build/Release \
  --output build/checks/ensemble-fixture --mode fixture --max-seconds 180
```

The default fixture mode replays only the six pilot source seeds, at 1000
histories per stream. It reports no acceptance intervals. Run the normal C++
and self-contained Python regressions before final collection. Commit the
protocol and runner, then run from that clean frozen revision:

```sh
TMPDIR="$PWD/build/checks/ensemble-tmp" PYTHONDONTWRITEBYTECODE=1 \
  python3 tests/fork_compat/source_observer_ensemble.py \
  --fork-root /path/to/pinned/reference --build-dir build/Release \
  --output build/checks/ensemble-final --mode final --workers 6 \
  --max-seconds 14400
```

Output must be a new/empty directory, outside the reference and within the
permitted writable root. The resolved revision and dirty state, complete
adapter/reference/production hashes, linked-core hash, compiler command and
version, executable hash, design and start time are retained. The manifest
marks final reservations spent before the first process starts. Do not use
`--mode final` for a benchmark or fixture. The command executes the current
resolved HEAD, without switching or fetching a branch.

Six isolated child processes maintain the six independent source streams.
Historical logging/RNG state is never shared between threads. Parallel worker
count affects scheduling only. Python streams 1 MiB blocks through gzip at
level 1 and hashes/counts raw records, without storing millions of dictionaries
or IDs. The C++ loop enforces contiguous IDs, per-primary validity and bounded
observables, and accumulates eleven sums with compensated long-double addition.
Every history retains its full source/flight record. Block boundaries do not
reset RNGs. `*.log` reports progress every 100000 histories.

A run-wide deadline kills child process groups; stream readers and archive
verification also check cancellation. Compiler execution uses the remaining
time budget. Failed/partial data are retained with failure status and no
successful confidence report. No history is replaced and no reservation is
resumed. Before reporting success, read back every gzip artifact to verify its
CRC and uncompressed SHA-256, require the expected raw-record count and six
valid summaries, and confirm source/core/reference hashes have not changed.
`report.json` contains the 39 intervals and separate numerical diagnostics.
Coordinate intervals are in units of L=0.25 Mpc; multiply their centers and
half-widths by L for Mpc. `metadata.tar.gz` compresses logs, source shims,
manifest, executable and summaries; per-history gzip files remain separate
in the same output directory, with raw/compressed hashes in verified summaries.
The archive and files stay usable after a failed collection as well.

The approved local execution yields native macOS evidence if run on macOS.
It does not substitute for a later native Linux amd64 execution. Automatic
Release CI runs the self-contained design, confidence and deadline guards
through `tests/test_fork_transport.py`; it does not spend these reservations.
