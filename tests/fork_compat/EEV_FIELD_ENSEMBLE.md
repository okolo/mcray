# Frozen EeV field/source ensemble protocol

On 2026-11-06 the user requested implementation of the proposed test and a
single Ubuntu 22.04 remote-server launcher, with an **eight-hour total limit**.
This implements the preceding 162-million-history proposal; the new limit
supersedes its four-hour resource proposal. Final streams remain unconsumed
until the user starts the launcher. No workstation exception is granted.

## Shared model

Retain the approved electrons, photons and protons at 1 EeV, initial weight 1,
source at the origin, three source profiles (+z pencil, axial 20-degree cone,
20-degree cone rotated by theta=0.4/phi=0.7 radians), and the identical
represented sphere center 0.1 Mpc along each declared axis with radius
0.01 Mpc. Use the existing complete divergence-free 3x4x5 field with zero
extension, z=0.01, H0=70, OmegaLambda=0.7 and flight horizon 0.2 Mpc.
Constant-rate propagation has no random/CEL interactions, adiabatic losses,
filters or thinning. Keep actual original source/field/deflector/move/observer
functions, both RNGs and the modern APIs; no physical model is changed.

Use **accuracy=40 in both implementations** for final collection, rather
than the cost pilot's 20. Existing 10/20/40 controls remain numerical diagnostics;
they do not prove a universal continuum error or supply an acceptance tolerance.
The comparison concerns the two algorithms at this fixed numerical setting.
Preserve the strict endpoint observer and known old/new precision/step-cap
exceptions; do not repair the reference to force matching results.

Public canonical prerequisites:

- [EeV deterministic refinement](../../docs/validation/fork_eev_species_comparison.md).
- [EeV cost-pilot results](../../docs/validation/fork_eev_field_source_pilot.md).
- [Original photon ensemble construction](SOURCE_OBSERVER_ENSEMBLE.md).

## observable family

Each complete primary is one sample. Include misses in all means, without
conditioning on arrival. For each of the nine profile/species cases, compare
modern mean minus historical mean for these fourteen bounded quantities:

1. Hit indicator H in [0,1].
2. Four arrival indicators H*I(tau<=t), in [0,1], with the original geometrical
   thresholds t=0.092,0.094,0.096,0.098 Mpc of c*time. Tau is final minus initial
   absolute time. These thresholds are reused from the frozen photon protocol,
   not fitted to pilot counts. Misses contribute zero.
3. Three final coordinates X_j/L in [-1,1], L=0.25 Mpc, for all histories.
4. Three final direction components in [-1,1], for magnetic deflection.
5. Three emitted direction components in [-1,1], retaining the source check
   from the photon protocol alongside the new final-direction moments.

This gives **126 difference intervals** (9*14). Add **18 individual hit-fraction
intervals** (nine cases, two implementations), for total simultaneous family
M=144. Observables from the same primary may be dependent; no extra independent
histories are inferred from their count. No pilot-derived bins, conditional
moments or adaptive statistics enter this family. Extra descriptive results
remain exploratory and cannot change the decision.

Record source weight/direction, final position/direction, actual elapsed time,
primary ID, hit flag, signed end-time error and numerical counters for each
history. Source/final weights remain exactly 1, energy/identity are checked,
source and inputs are hash-verified. Bounds are strict validity checks; never
clip/renormalize an unexpected value or count a failed primary as a miss.
A missing, duplicate, nonfinite, failed or out-of-bounds history invalidates
the comparison. Preserve partial data; do not retry or replace histories.

## confidence construction and counts

Reuse two-sided Hoeffding bounds with equal error allocation alpha/M and the
union bound, alpha=0.05, assuming the standard independent-primary/independent-
source-stream Monte Carlo model. Independence across observables is unnecessary.
RNG quality remains an assumption. For range width w=b-a and N histories in
each implementation:

    h_difference = w * sqrt(log(2*M/alpha)/N)
    h_individual_hit = sqrt(log(2*M/alpha)/(2*N))

Use closed untruncated intervals; even deterministic/all-hit cases keep the
prescribed positive width. Target the earlier absolute hit-fraction difference
half-width <=0.001 (0.1 percentage point). The worst-case minimum N is
ceil(log(5760)/0.001^2)=8658693. Use **N=9000000 per source stream**, eighteen
streams and **162000000 total histories**. Do not tune N from pilot variance.

| Quantity | Proposed half-width |
| --- | ---: |
| Hit or arrival-indicator difference | 0.0009808552251 |
| Individual hit fraction | 0.0006935693810 |
| Emitted/final direction-component mean difference | 0.0019617104502 |
| Final coordinate mean difference | 0.0004904276126 Mpc (0.490428 kpc) |

The 0.001 target applies to hit/arrival fractions; moment precision is stated
separately. These are statistical uncertainty bounds for the frozen numerical
algorithms, not integration-error bounds or practical equivalence margins.
The existing fixed-state refinement does not justify claiming that every
numerical bias is below these widths over the full source distribution.

If every validity/precision check passes and **all 126 difference intervals
contain zero**, report compatibility with zero difference for this frozen
family at the stated resolution. Do not claim exact or whole-program physical
equivalence. If any interval excludes zero, retain/report that discrepancy;
do not widen intervals, change accuracy/bins/family, add histories or rerun a
reservation to seek a preferred outcome. An incomplete or invalid run is
inconclusive, never successful. Known old observer time overruns remain
separate numerical diagnostics rather than rejection thresholds.

## reservations and execution

All following seeds are frozen, unconsumed and disjoint from the previous
pilot/final and software-fixture reservations in this migration:

| Profile | Species | Fork seed | Modern seed |
| --- | --- | ---: | ---: |
| pencil_z | electron | 2010101 | 2010102 |
| pencil_z | photon | 2010103 | 2010104 |
| pencil_z | proton | 2010105 | 2010106 |
| cone_z | electron | 2010107 | 2010108 |
| cone_z | photon | 2010109 | 2010110 |
| cone_z | proton | 2010111 | 2010112 |
| cone_rotated | electron | 2010113 | 2010114 |
| cone_rotated | photon | 2010115 | 2010116 |
| cone_rotated | proton | 2010117 | 2010118 |

Use separate fixture seeds 2011101-2011118 with 1000 histories per stream for
software checks. Fixtures are explicitly outside the final family and may be
replayed for software validation. Never smoke-test final reservations.
Maintain one uninterrupted source RNG per case/implementation; block boundaries
or scheduling do not reseed it. Modern propagation seed 17 consumes no draws
from interactions in this model. A failed final attempt remains spent and
needs separately reviewed recovery; there is no silent resume/replacement.

Use an external **native Linux amd64** server, at most **16 isolated worker
processes** for the eighteen source streams, and a **eight-hour whole-run limit**.
This is not a new workstation-run exception. Reuse the existing process/block
streaming approach: original global logging/RNG state is never shared across
threads. Capture each complete history directly into streaming gzip level 1
and compensated long-double summaries. Do not retain millions of Python
objects or primary IDs; verify contiguous IDs online and output counts exactly.

At startup resolve and record the exact requested modern Git commit and the
immutable historical snapshot, compiler/dependency/image versions, adapter/
production/core/reference/input hashes and seeds. Run required C++/Python
regressions and a disjoint fixture before final collection. Require a committed
clean candidate, new empty output under configured test-data storage and
sufficient disk (propose at least 64 GiB free). The runner must be public and
self-contained, without dependency on this private workspace. Preparing this workflow does not push code or start a remote calculation.

Hash/count raw bytes while streaming; keep compressed final or partial records,
logs, exit status, exact manifests, extracted source and summaries. A global
deadline must stop process groups and readers. Before success, decompress each
archive to verify gzip CRC and raw SHA-256, recheck source/input hashes and
require all eighteen complete summaries and expected counts. Compress metadata
and retain everything in the configured test-data directory. A resource/error
failure produces no successful statistical report.

## Cost and storage scope

Scaling the short native macOS pilot at accuracy 20 to the frozen counts
gives roughly 41.6 serial compute minutes or 68.3 process minutes including
its text output. Accuracy 40, compression, scheduling and server hardware were
not benchmarked; neither a factor-of-two accuracy scaling nor ideal 16-worker
speedup is established. The eight-hour server cap is a frozen resource limit,
not a measured runtime prediction. Do not launch this campaign locally.

Pilot-size scaling gives about 103.3 GiB of raw JSONL and 12.9 GiB compressed
with the pilot archive's gzip level 9. Proposed streaming gzip level 1 and
longer primary IDs change storage/CPU costs; these figures are estimates,
not capacity guarantees. No full raw file needs to be materialized. Fail
cleanly with preserved partial data if disk/time limits prevent completion.


## Entry points and reporting

Use [the Ubuntu server instructions](EEV_SERVER.md) and the single standalone
run_eev_field_server.py entry point. It resolves an exact supplied commit,
configures/builds Release, runs the complete short CTest suite and Python
guards, then a disjoint 18000-history fixture before final collection.
Dependencies are prepared once as documented, outside the timed invocation.
The launcher includes source/reference downloads, build/checks, collection,
archive verification and metadata packaging within its eight-hour budget.

The lower-level eev_field_ensemble.py defaults to fixtures and supports a
prepared build; final mode requires a persistent exclusive ledger. The
launcher uses one stable ledger in its configured test-data directory. Do
not remove/change that directory to bypass a spent reservation. A report
is written only after complete collection/verification and metadata packaging.
No partial report is accepted; preserve failures and obtain review before
recovery. Final mode draws all eighteen streams; it is not a benchmark option.

A validated collection can detect a discrepancy. Server status completed
means the collection/verification finished; inspect report.json's
compatible_at_reported_resolution and detected_differences for the scientific
result. Fixture reports contain no acceptance intervals. The raw gzip files
remain separate from compressed metadata to avoid duplicating large payloads.
