# Electron/photon/proton field/source cost pilot at 1 EeV

The user-authorized pilot completed on 2026-11-04 with **18000 complete primary
histories and zero failures**, from clean revision
`f60b415316e1bb998fdbd73aaeda3ef0bfb53538`. The
[frozen protocol](../../tests/fork_compat/EEV_FIELD_SOURCE_PILOT.md) owns the
model, profiles, counts, new source reservations and failure/deadline rules.
[Machine-readable evidence](fork_eev_field_source_pilot.json) retains the
raw results, descriptive groups, timings, source/input/core/reference hashes
and software-check metadata. No statistical-equivalence verdict is assigned.

## Scope and runtime

Electrons, photons and protons have energy 1 EeV and initial weight 1. Reuse
the +z pencil, axial 20-degree cone and rotated 20-degree cone, the existing
tabulated magnetic field, accuracy 20, common cosmology, source at the origin
and sphere at distance 100 kpc/radius 10 kpc. The constant-rate path has no
random/CEL interactions or energy losses. Historical functions execute in the
previously documented reduced driver; modern uses the actual current engine.
Only validation tooling changed for species/energy selection.

The native macOS arm64 process took **0.455 seconds** including startup/output.
Sampling/propagation timers sum to **0.277 seconds**. Compilation, execution
and reduction took **7.084 seconds** in total, below the 900-second limit.
The linked evidence retains exact timings. These short measurements are
rough cost estimates on this workstation, not native Linux benchmarks or
stable speed comparisons between different driver architectures.

## Descriptive counts and costs

Each row has exactly 1000 histories per implementation. Compute timings
include source sampling, per-primary engine/deflector setup, propagation and
counters. They exclude grid loading, source-stream setup, JSON output,
compilation and startup. No field values, source directions or weights were
renormalized for this comparison.

| Profile | Species | Fork hits | Modern hits | Fork compute, ms | Modern compute, ms |
| --- | --- | ---: | ---: | ---: | ---: |
| +z pencil | electron | 1000 | 1000 | 16.972 | 22.505 |
| +z pencil | photon | 1000 | 1000 | 1.170 | 6.015 |
| +z pencil | proton | 1000 | 1000 | 16.840 | 21.088 |
| Axial 20-degree cone | electron | 90 | 93 | 21.125 | 20.574 |
| Axial 20-degree cone | photon | 75 | 87 | 2.841 | 9.858 |
| Axial 20-degree cone | proton | 61 | 77 | 16.952 | 21.395 |
| Rotated 20-degree cone | electron | 72 | 75 | 19.849 | 23.505 |
| Rotated 20-degree cone | photon | 80 | 85 | 3.221 | 10.359 |
| Rotated 20-degree cone | proton | 72 | 75 | 20.031 | 23.039 |

Cone streams have independent seeds and different RNG implementations.
Unequal counts are descriptive pilot results; no confidence family or final
acceptance criterion was applied. Pencil directions are deterministic, so
their repeated histories do not measure source-distribution uncertainty.
All ordinary-cone weights remain one; no zero or negative emitted weights
occur. Maximum initial direction-norm residual is 2.2205e-16.

There is one historical end-time overrun, in the rotated-cone proton stream,
by **1.1509361544 pc** of flight parameter. Modern has none. This retains the
known [historical observer step-floor behavior](fork_observer_comparison.md)
as a diagnostic rather than dropping or replacing the primary. Maximum outer
engine moves are 2182 for electrons, 2042 for protons and 2032 for photons;
these counts do not include internal magnetic integration stages. No old
zero-scale calls occur. Per-history energy/identity and modern weight checks
pass. The minimal historical shell retains source weight separately and has
no weight-changing propagation operation.

## Software checks and retained data

Twenty-three self-contained Python guards pass through the existing Release
CI entry point. They cover distinct EeV/PeV and fixture/pilot reservations,
unchanged grid reuse, revised species/counts and complete-primary accounting.
The new --design eev selects the explicit EeV model; the default PeV model
and its frozen earlier results remain separate.

Before pilot collection, separate 36-history fixtures passed native macOS
Release, emulated Linux amd64 Release and Linux ASan+UBSan with leak detection.
No sanitizer diagnostics occurred. Fixture seeds 2008101–2008118 are disjoint
from the pilot and are not statistical pilot evidence. All checks used clean
f60b415. Production algorithms were unchanged, so no full CTest suite or
previous statistical ensemble was rerun for these validation-tool changes.
Linux checks use the pinned offline image recorded in the evidence and run
under amd64 emulation on the arm64 workstation.

Pilot seeds **2007101–2007118 are now spent**. Do not rerun them or include
them in a future final ensemble. The new/empty output policy preserves every
attempt. Inputs, extracted historical fragments, compiler commands/version/
logs, manifest, executable, raw JSONL and archive remain under ignored
build/checks/eev-field-source-pilot-native. The 12325677-byte raw JSONL
compresses to 1534820 bytes; full gzip decompression reproduces SHA-256
`f0f44b415cd603af0fdd68873b9b80bbb844a183577a6e8cca04c5013825ba4a`.
Startup interaction-table warnings are retained; the probe does not use those
tables. The successful report requires all 18000 primaries and eighteen
source-stream timing records.

## Cost implications and remaining decisions

A linear compute-only extrapolation gives roughly **17–24 seconds per million
charged histories per stream**, and **1–11 seconds per million photons per
stream**, on this workstation. Scaling the whole eighteen-stream matrix to
one million histories per stream gives about 4.6 minutes of serial compute;
linearly scaled whole-process cost is about 7.6 minutes. These short-pilot
extrapolations have no scaling guarantee and do not select a final count.
Compression, storage, different hardware and native Linux execution change
costs. Runs expected to exceed 30 minutes belong on an external Linux server.

The [fixed-state EeV refinements](fork_eev_species_comparison.md) retain their
numerical limits, including the shared integration/observer stopping sensitivity
and the distinct Linux 10 pc boundary-step difference. This pilot does not
choose a new accuracy tolerance or turn those diagnostics into an error bound.
Final statistical observables, bounds, confidence/precision choices, counts,
new seeds and stopping/failure rules still require agreement before collection.
