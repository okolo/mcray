# Electron/photon/proton field/source cost pilot at 1 EeV

On 2026-11-04 the user authorized this pilot: 1000 primary histories per
profile/species/implementation for electrons, photons and protons at 1 EeV.
There are three profiles and two implementations, hence **18000 histories**.
This is a separately reserved cost/failure diagnostic without intervals or
an equivalence verdict. The completed 1 PeV pilot remains frozen and separate.

## Frozen model and reservations

Reuse [FIELD_SOURCE_PILOT.md](FIELD_SOURCE_PILOT.md) for geometry, source
sampling, immutable tabulated field, cosmology, accuracy=20, constant-rate
engine, per-history guards, retained artifacts and the 900-second whole-run
budget. There are no random/CEL interactions or energy losses. Source profiles
are +z pencil, axial 20-degree cone and rotated 20-degree cone with axis
(theta,phi)=(0.4,0.7) radians. Energy is 1e18 eV for all three species.

| Profile | Species | Fork seed | Modern seed |
| --- | --- | ---: | ---: |
| pencil_z | electron | 2007101 | 2007102 |
| pencil_z | photon | 2007103 | 2007104 |
| pencil_z | proton | 2007105 | 2007106 |
| cone_z | electron | 2007107 | 2007108 |
| cone_z | photon | 2007109 | 2007110 |
| cone_z | proton | 2007111 | 2007112 |
| cone_rotated | electron | 2007113 | 2007114 |
| cone_rotated | photon | 2007115 | 2007116 |
| cone_rotated | proton | 2007117 | 2007118 |

Each row has exactly 1000 primaries per implementation. Each case uses one
uninterrupted source RNG stream. Modern propagation seed 17 consumes no
interaction draws. No resampling, retries, replacement histories, adaptive
counts or reuse of previous pilot/final reservations is allowed. Any failed,
missing, duplicate or nonfinite primary invalidates the whole pilot.

The separate software fixture has two primaries per stream (36 total), with
seeds 2008101–2008118. It may be replayed for software checks; it is never
pilot or statistical acceptance evidence. Both reservations are disjoint
from the earlier 1 PeV and photon-baseline streams.

## Reproduction

Use a compatible configured build and the pinned read-only historical snapshot
as described in AGENTS.md. From the public repository root:

```sh
mkdir -p build/checks/eev-pilot-tmp
TMPDIR="$PWD/build/checks/eev-pilot-tmp" PYTHONDONTWRITEBYTECODE=1 \
  python3 tests/fork_compat/field_source_pilot.py --design eev \
  --fork-root /path/to/pinned/reference --build-dir build/Release \
  --output build/checks/eev-field-source-fixture
```

Explicit **--mode pilot** selects the authorized 18000-history reservation;
use a different new/empty output directory and a committed clean source tree.
The default --design pev retains the old model and reservations; do not use
its spent pilot streams again. --command-prefix supports prepared offline
Linux compiler/runtime checks. Preserve compiler, executable, source/core,
reference/input hashes, raw histories, timings and logs. A successful report
requires all 18000 histories and eighteen timing records.

The [1 EeV fixed-state refinement](../../docs/validation/fork_eev_species_comparison.md)
retains its numerical limits and known boundary-step differences. This pilot
adds real source sampling but does not turn that refinement into a numerical
acceptance bound. Final ensemble observables, confidence/precision choices,
counts, new seeds and stopping rules still require agreement.
