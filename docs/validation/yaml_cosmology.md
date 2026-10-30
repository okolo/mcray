# YAML cosmology verification

This change exposes the existing H0 and OmegaLambda values through the shared
YAML configuration without changing the flat matter-plus-Lambda model or its
defaults. The [usage contract](../usage/cosmology.md) defines units, domains,
global-state handling and replay.

## Checks

Seven new CTest cases cover:

- CLI/default values and partial overrides in both YAML schema versions;
- invalid, nonfinite, quoted, singular and unsupported fields and duplicate
  keys, rejected before creating output;
- real executable runs in fresh processes with (H0,OmegaLambda) equal to
  (67.8,0.692), (70,0.692), (70,0.7) and (70,0), plus snapshot replay;
- rejection of conflicting already-initialized global physical settings;
- parameter identity and redshift-range guards on a local cosmology object;
- fatal internal H0 conversion underflow without a usable initialized state;
- an actual executable rejecting a finite positive H0 whose internal rate
  underflows the double-precision loss interface, before output creation.

The executable control propagates two massless neutrino primaries without
interactions or magnetic fields from redshift 0.01 to zero. A zero-opening cone
selects the existing spatial output, whose final column is actual FlightMpc.
The measured flight is compared with the library time mapping initialized
explicitly with the requested parameters. This checks application wiring,
not an independent validation of the cosmological integration formula.
Tolerance scales with time precision and the universe age to include rounding
from subtracting large absolute times. The H0-only flight ratio is checked
against 67.8/70; the Lambda-only change must affect flight. Particle output
must replay byte-for-byte from each saved snapshot in a fresh process.
The snapshot values and run.yaml configuration hash are checked too.

An existing CLI/YAML/default-snapshot seeded-history check now also verifies
that the default cosmology is exported. The complete commented version-1
example contains the new mapping and remains part of the existing input tests.

## Results

On 2026-10-30, the complete short suite passed on the production change:

| Environment | Complete CTest | Elapsed |
| --- | --- | --- |
| macOS arm64, Clang 22.1.5, Release | 241/241 | 17.22 s |
| Ubuntu 22.04 Linux amd64, GCC 11.4.0, Release | 241/241 | 115.50 s |
| Same Linux build with ASan/UBSan and leak detection | 241/241 | 166.82 s |

The final test harness additionally rejects sanitizer diagnostics in child
logs even when a nonzero application exit is expected. After this harness-only
change, all seven affected cosmology cases passed again on macOS Release
(0.34 s), Linux Release (7.81 s), and Linux ASan/UBSan (12.12 s). No production
code changed after the full-suite results above.

The Linux runs used the existing offline Ubuntu image
`sha256:37682238eeaf9f01f7eddea17d4e55551d2e317307dfa5c8ea8dd0dc8084baf6`
under amd64 emulation on the arm64 workstation. These are not native Linux
performance measurements. The sanitizer configuration instruments the modern
C/C++/Fortran project with address/undefined sanitizers, fail-fast recovery,
frame pointers and non-PIE executables. It uses leak detection and no
suppressions. All checks were bounded below 30 minutes; build/temp/cache
artifacts stayed inside the project.

A preliminary sanitizer run exposed a stale Cosmology object built before
the final rate-range guard. The new unchanged underflow regression failed
because the child incorrectly completed and created output. Rebuilding the
final source fixed that same test and the full suite; no conditions or
tolerances were relaxed. The original failed log is retained locally.

A separate before/after default regression used the previously verified
Release executable and the new executable with the same omitted-cosmology
input: 12 seeded, noninteracting, undeflected neutrinos, a rotated 20-degree
cone, redshift 0.01 and a power-law source spectrum. Particle files were
byte-identical before/after on each platform (not between platforms):

- Native output SHA-256: `9a9a1504c5bfba7ab52069d8160e421396a12c236087bfdb32a0d9b56fe07575`.
- Linux output SHA-256: `cf2e763c91194b9a517d89927e76f66aff6936ab1b24589d7503d4fab721c47a`.

No astrophysical acceptance or historical-fork whole-cascade equivalence is
asserted by these deterministic controls. Local ignored logs, reference
executable hashes and inputs are retained under `build/checks`; ordinary CI
uses the self-contained C++ tests and does not require the previous binary.

## Reproduction

Configure/build with the supported toolchain, then run:

```sh
ctest --test-dir build/Release -R Cosmology --output-on-failure
ctest --test-dir build/Release --parallel 2 --output-on-failure
```

The checked production-source hashes are:

| File | SHA-256 |
| --- | --- |

| `src/lib/Cosmology.h` | `23aaa2a6e66907393b80f54f369faad87f2a49102314aec57ea6d1ee34663c18` |
| `src/lib/Cosmology.cpp` | `72f0246254e150a3384b3c14ab8549debc6c689e162f031cac19f0213ecbef1d` |
| `src/app/crbeam/Configuration.cpp` | `bcd0cfbf67e4941439bd287129b42e2c7f709a340ed3ec4e923507b7ab301f22` |
| `src/app/crbeam/CRbeam.cpp` | `d768d23f772f837de0be160fad61109fcc7e1f06cfb8446a416c26e500613388` |
