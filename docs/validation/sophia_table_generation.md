# Independent SOPHIA table generation for GZK validation

## Discovery and provenance

On 2026-10-07 the original predecessor table generator was found at
/Users/ok/propagation/sources/sophia. Earlier missing-generator statements
referred to the TransportCR checkout and adjacent project copies; this
additional source location removes that dependency.

The parent propagation repository HEAD is
c27c442ce86eff5585c609557fbdda73c641a10a. The frozen source-file SHA-256
manifest is authoritative for this build; the generator banner contains
VERSION_STR rather than a useful Git identifier. Do not infer a clean or
matching historical binary from the repository HEAD or banner.

The implementation generates secondary energy-fraction histograms by
calling its own SOPHIA Fortran sources. It is not derived from the mcray
sampler. Both codes still share the underlying SOPHIA physical model, so
this is not an independent validation of that model itself.

## Reproducible build

Use tests/transportcr/build_sophia_generator_docker.sh with a source-root
argument and a fresh output directory outside Git. It copies TParticle.h,
Sophia.cpp/h, GitVersion.cpp, sophia2cpp.f, and the five original Fortran
sources; records their hashes and the Docker image ID; and builds with
GCC/GFortran in Linux amd64. Assertions are enabled with -D_DEBUG.
No source, table, or binary in the original propagation/TransportCR trees
is modified. The tested working copy is
../validation/sophia-generator-20261007.

The generator's seed initializes C rand and selects a deterministic sequence
of 1000 warm-up SOPHIA calls; it does not directly set the Fortran RNG seed
(the direct seed call is commented out in the original source). Preserve
this behavior. Reproducibility is demonstrated in this Linux image; it is
not claimed across different C libraries or architectures.

## Diagnostic runs

Separate proton and neutron pilots use 200 bins per decade, kmax=1e7 GeV,
100 events per photon-energy cell, and seeds 2026100701/2026100702.
Both completed successfully in Docker Linux amd64, taking 19.29/19.43
seconds wall time while run concurrently. The input photon grid contains
1564 proton cells and 1565 neutron cells. All 18 secondary files passed
format/range/count validation; no inconsistent empty rows were detected.
A repeated proton pilot produced byte-identical nine secondary files and
the cross-section file with the same seed.

Manifests, logs, sources, executable, and diagnostic data are retained in
../validation/sophia-generator-20261007, outside Git. This avoids relying
on temporary files that may disappear between sessions.

A native TransportCR loader/propagation smoke run with the diagnostic
200-bin tables is prepared in ../validation/gzk-smoke200-20261007. It is
strictly a technical format test. It completed successfully in 149 seconds
(130 seconds of propagation). The table statistics are far too small
for spectral acceptance or a 0.2% convergence claim.

## Approved full-statistics generation

Existing proton secondary tables have the following header statistics:

| Bins per decade | Events per photon-energy cell | Photon cells | Approximate maximum k, GeV |
| --- | ---: | ---: | ---: |
| 20 | 1380000 | 157 | 9.6681e6 |
| 50 | 2000000 | 356 | 1.9158e6 |
| 100 | 1000000 | 712 | 1.9559e6 |

The old 100-bin generation log reports 232064 seconds for the proton run
on its historical machine. The current pilot roughly extrapolates to
45--55 hours for concurrent p/n runs with 200 bins per decade, one million
trials per cell, and kmax=2e6 GeV. This covers the existing fine table's
upper range. The estimate is not a benchmark or completion guarantee;
event cost varies with photon energy and I/O overhead does not scale
linearly with trials.

The user approved the full p/n generation on 2026-10-07. Both detached
Linux amd64 jobs were launched at 20:33:59 UTC with 200 bins per decade,
one million trials per cell, kmax=2e6 GeV, and seeds 2026100711/2026100712.
Container names are mcray-sophia200-p-20261007 and
mcray-sophia200-n-20261007. They write full200-p/full200-n below the frozen
working directory. full-generation.json records the commands and image/
binary identities; each run publishes exit-code.txt on exit.

They are long-running jobs, not completed references. The generator has
no checkpoint/resume protocol. Interrupted runs must not be resumed by
mixing partial data or overwritten without preserving their provenance.
Do not turn the diagnostic tables into a reference, increase tolerances,
or replace absent counts with interpolated pseudo-events.

The generation is itself a Monte Carlo calculation. Increasing energy-grid
density also changes table noise; a difference between two table densities
is not pure discretization error. Before claiming reference convergence,
separate spatial/grid changes from the stochastic uncertainty of generated
histograms, using independently generated replicas or a separately approved
method. The necessary replication/statistics will be assessed from pilot
measurements; one million trials does not by itself prove 0.2% accuracy.

After generating suitable independent tables, run the existing unchanged
TransportCR with 100/200-bin references and the fixed 28 GZK observables.
Keep the approved log-cell reconstruction, 0.2% reference criterion,
99% Bonferroni intervals, and 5% spectral budget. Agree final MC acceptance
statistics only after reference convergence is resolved.


## Completion and validation

Run tests/transportcr/validate_sophia_tables.py only after both exit-code
files record success. It validates primary/secondary identities, density,
trials, nonnegative counts, energy-fraction index range, identical photon
grids within each primary, coverage through the requested kmax, and matching
cross-section grids. It writes a new checksum manifest only after both
runs pass. The tool was checked on the complete diagnostic runs and rejects
the running full-generation directories. A valid format manifest does not
establish table sampling accuracy.

Different supplied seeds select warm-up sequences; they do not directly
reseed Fortran's RNG. Byte reproducibility is demonstrated, but independence
of future same-primary replica tables must be established separately. Do
not estimate table variance assuming such replicas are independent solely
because their integer seeds differ.
