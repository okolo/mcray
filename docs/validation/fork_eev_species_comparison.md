# Electron, photon and proton comparison at 1 EeV

The user requested this species/energy check on 2026-11-01. All **729 fixed-state
cases** completed in each of native macOS arm64 Release, emulated Linux amd64
Release and Linux ASan+UBSan with leak detection. Hit flags agree in every
old/new pair. No sanitizer diagnostics occurred. The check reports numerical
refinement and historical differences without a statistical-equivalence verdict.

The [coupled protocol](../../tests/fork_compat/COUPLED.md) owns fields, geometry,
cosmology, extraction and stopping/failure rules. Select
**--species electron photon proton**; energy remains 1 EeV, with accuracy
10/20/40. This separate 729-case matrix replaces positrons with photons.
Earlier frozen results and the completed 1 PeV source cost pilot remain intact.
There is no source RNG, new statistical reservation, interaction or energy loss.

[Machine-readable evidence](fork_eev_species_comparison.json) records the full
clean source revision `03bb431ce588bffba2aa01c6bd52fdbf4f3203ff`, pinned reference,
input/adapter/core hashes, compiler versions, offline execution prefixes,
per-species/field differences and both implementations' refinements. Production
code was not changed. Default coupled inputs preserve the original matrix;
22 self-contained Python guards passed, including the requested species design.

## Old/new agreement and retained stopping difference

Each species has 243 cases. Both implementations record 120 electron hits,
99 photon hits and 120 proton hits in each environment. These are counts over
fixed geometrical probes, not estimates from a Monte Carlo source distribution.
No hit flag changes between accuracy 10/20 or 20/40. No end-time overruns or
historical zero-scale calls occur in this matrix.

| Environment | Hit differences | Maximum old/new position difference, Mpc | Maximum flight difference, Mpc |
| --- | ---: | ---: | ---: |
| macOS arm64 Release | 0/729 | 5.3919e-16 | 0 |
| Linux amd64 Release, emulated | 0/729 | 1.0000e-5 | 1.0000e-5 |
| Linux ASan+UBSan, emulated | 0/729 | 1.0000e-5 | 1.0000e-5 |

On macOS, photons agree exactly in all fields. With an active field, the
maximum charged-particle position difference is 5.3919e-16 Mpc on macOS and
1.9527e-17 Mpc on Linux. Linux photon differences in installed-field cases
remain free-motion stopping differences: the neutral deflector applies no force.

The Linux 10 pc maximum is the previously explained
[first-step rounding amplification](fork_field_observer_comparison.md#linux-one-step-boundary-amplification).
Double subtraction in the old observer lands slightly inside the axial sphere;
long-double subtraction in the modern observer lands exactly on its surface.
Under strict inside membership, modern then takes the 0.001*R=10 pc minimum
step. This changes the endpoint/time, while both record a hit.

There are 30 different outer engine-move counts: six no-field axial cases
per charged species and eighteen photon cases (six in each field selection).
Photon motion is identical across field selections. All other move counts
agree. This preserves the known precision/step-floor distinction rather than
changing production arithmetic to manufacture matching boundary endpoints.

## Refinement inside each implementation

The table gives maximum modern **position changes in pc**, over the 27 initial
states for each field/species. Each environment retains all 486 accuracy-pair
records, including directions, time changes and the corresponding old results.

| Species | Field | Accuracy pair | macOS | Linux Release / ASan+UBSan |
| --- | --- | --- | ---: | ---: |
| Electron | Constant | 10→20 | 10.122 | 10.122 |
| Electron | Constant | 20→40 | 3.629 | 3.629 |
| Electron | Tabulated | 10→20 | 32.269 | 32.269 |
| Electron | Tabulated | 20→40 | 16.594 | 16.140 |
| Proton | Constant | 10→20 | 10.222 | 10.222 |
| Proton | Constant | 20→40 | 3.602 | 3.602 |
| Proton | Tabulated | 10→20 | 33.711 | 33.711 |
| Proton | Tabulated | 20→40 | 17.986 | 17.514 |
| Photon | Any | Both pairs | 0 | 0 |

No-field charged motion also has zero refinement changes. Old refinements
follow the same scales; the exact values and per-implementation hit-change
flags are in the evidence. Photon invariance is expected because there are
no interactions and no magnetic force.

Charged endpoint changes decrease over the two refinement pairs in these
aggregate maxima. This is not a continuum-error bound, universal convergence
rate or proof that accuracy=40 is sufficient for an unspecified observable.
Observer stopping and magnetic integration are coupled: even with unchanged
hit flags, different stopping times affect endpoint changes. Platform rounding
and the shared minimum observer step remain visible. No acceptance threshold
or additional accuracy level was chosen after seeing these results.

## Reproduction and retained evidence

```sh
mkdir -p build/checks/eev-species-tmp
TMPDIR="$PWD/build/checks/eev-species-tmp" PYTHONDONTWRITEBYTECODE=1 \
  python3 tests/fork_compat/coupled.py --species electron photon proton \
  --fork-root /path/to/pinned/reference --build-dir build/Release \
  --output build/checks/eev-species-check
```

Use the immutable public snapshot and compatible configured build described
in AGENTS.md. The compiler is limited to 240 seconds and execution to 180
seconds. Ordinary builds/CI do not need this optional external snapshot.
Raw 21-digit JSONL, verified gzip archives, exact fragments, inputs, compiler
commands/versions/logs and manifests remain under ignored
build/checks/eev-species-{native,linux-release,linux-sanitizers}. Archive
round trips reproduce their raw SHA-256 values. Startup interaction-table
warnings are retained; this interaction-free probe does not read those tables.
Linux timings are emulated, not native server measurements.

This diagnostic does not rerun spent source seeds or establish sampled-source,
interacting-cascade or whole-program equivalence. In particular, it checks
1 EeV as requested; it does not relabel the previous 1 PeV pilot or supply
its still-unperformed refinement study.
