# Grid magnetic field implementation checks

On 2026-10-30, the grid-field and composition C++ API passed the short
implementation checks below. The [usage contract](../usage/grid_magnetic_fields.md)
describes units, interpolation, extensions, ownership and current application
scope. These results do not establish physical acceptance of an imported
field realization or turbulence model.

The initial boundary and step-scale rules described below were subsequently
aligned with the fork in the [direct kernel comparison](fork_kernel_comparison.md).
The usage contract describes the current behavior; the results below retain
their original tested-source scope.

## Source and environments

The implementation series starts from main
`ff55a22db8cd65c2ea89435276c28401f0217dce`. The checked final files have SHA-256:

| File | SHA-256 |
| --- | --- |
| `src/lib/GridMagneticField.cpp` | `56bf2c4d2980f850498a8b0c4abd9679a4de97c4fa09f84f457040dce517ed02` |
| `src/lib/GridMagneticField.h` | `023f265668b667d9e540b9a201f3707046848b7c1fa35999d1400c981832f2e3` |
| `tests/grid_magnetic_field.cpp` | `917bc83e6ede933cd103656180944be135af18de486caea1ffa179de67a6ccb1` |

Native environment: macOS arm64, Homebrew Clang 22.1.5, GNU Fortran 15.2.0,
C++17/libc++, Release build and the existing Homebrew dependencies. The
configuration uses the toolchain and search paths from `bin/build.sh`, with
an offline `FETCHCONTENT_SOURCE_DIR_GOOGLETEST` override pointing to the
already extracted, pinned GoogleTest 1.18.0 source. The build directory for
these retained local checks is `build/verified`.

Linux environment: Ubuntu 22.04 linux/amd64, GCC 11.4.0, Docker image
`sha256:37682238eeaf9f01f7eddea17d4e55551d2e317307dfa5c8ea8dd0dc8084baf6`.
The container ran on the macOS arm64 workstation under emulation, without
network access and with the source checkout mounted read-only. This is not
native Linux performance evidence.

## Results

| Check | Result |
| --- | --- |
| Native Release build of all targets | Passed. |
| Native final CTest suite, parallelism 2 | 164/164 passed; 15.03 seconds elapsed. |
| New magnetic-field tests | 9/9 passed natively. |
| Linux CMake build of production library and new test target | Passed. |
| Linux targeted CTest checks | 9/9 passed; 4.38 seconds elapsed under emulation. |
| Final-source macOS ASan/UBSan checks | 9/9 passed; no reported sanitizer error. |
| Final-source Linux ASan/UBSan checks | 9/9 passed; no reported sanitizer error. |

The Linux sanitizer check covers the final fixed-underlying-type boundary
mode enum and the complete final test file. The sanitizer builds instrument
`GridMagneticField.cpp` and `grid_magnetic_field.cpp`, link the existing
normally compiled mcray/GoogleTest dependencies, and use
`-O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined -fopenmp`.
They run with `ASAN_OPTIONS=detect_leaks=1` and
`UBSAN_OPTIONS=halt_on_error=1`. They do not instrument the entire historical
library or external runtimes. Concurrent-read tests use four OpenMP workers;
no ThreadSanitizer result is asserted.

The nine cases cover:

- Node, face, edge and affine interior interpolation on unequal grid spacings.
- Positive, negative and fractional box origins; unchanged values at different
  redshifts; finite positive minimum spacing even outside a zero-extended box.
- Zero, clamped, periodic and mirrored extension, including negative wraps,
  multiple periods and distant finite coordinates.
- Missing files, invalid modes, malformed/nonfinite/extra columns, incomplete
  and duplicate nodes, compensating missing/duplicate rows, irregular axes
  and degenerate dimensions; all reported errors are nonretryable.
- Invalid query coordinates and vector sizes.
- Component ownership, empty/null compositions, invalid scales and sum
  overflow, including a clone outliving the original deflector and caller's
  field references.
- Concurrent queries of the shared immutable composition.
- Uniform-field propagation compared with a separate constant-field control
  for an electron and an explicit Fe56 identity, including clone behavior,
  energy, weight, coordinates, direction and time.

All inputs are small synthetic files created in temporary directories and
removed by the tests. No external magnetic maps, long cascades, FFT fields,
new backgrounds or changed physical processes are involved. Existing output
schemas, CLI/YAML defaults, RNG and nuclear identity are unchanged.

## Reproduction and retained logs

Configure/build with the supported toolchain as described in
[tests/README.md](../../tests/README.md), using its offline GoogleTest option
when needed. After configuration:

```sh
cmake --build build/verified --parallel 2
ctest --test-dir build/verified --output-on-failure --parallel 2
ctest --test-dir build/verified -R '^unit.magnetic\.' --output-on-failure
```

The local checkout retains build, CTest and sanitizer logs under the
untracked `build/grid-field-verification` directory. Those logs and binaries
are local artifacts and are not dependencies of the public tests. There is
no claim of a native Linux server ensemble or statistical physics acceptance.
