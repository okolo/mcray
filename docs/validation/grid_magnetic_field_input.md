# Grid magnetic field CLI and YAML checks

The application input extension adds the grid-field CLI and magnetic-field
schema version 2 described in the
[usage contract](../usage/grid_magnetic_fields.md). Existing CLI and version 1
remain available. The implementation starts from library commit
`02a8bed9115ff7f99c3ad594568c041ab5fc9849`. The independent buffer-release
repair is commit `1150a1a`.

## Short checks on 2026 10 30

| Check | Result |
| --- | --- |
| Native macOS arm64 Release build | Passed. |
| Complete native CTest suite | 170/170 passed, parallelism 2, 15.23 seconds. |
| Ubuntu 22.04 linux/amd64 CMake build | Production library, application, input tests and grid tests built. |
| Linux targeted CTest checks | 51/51 passed, parallelism 2, 7.34 seconds under emulation. |
| Native application-input ASan/UBSan checks | 16/16 configuration cases passed with no reported sanitizer error. |
| Final-source Linux application-input ASan/UBSan checks | 16/16 passed after the independent buffer-release repair. |

Native tooling is Homebrew Clang 22.1.5, libc++, GNU Fortran 15.2.0, C++17
and the supported Homebrew dependency paths. Linux uses GCC 11.4.0 in image
`sha256:37682238eeaf9f01f7eddea17d4e55551d2e317307dfa5c8ea8dd0dc8084baf6`.
Docker runs linux/amd64 under emulation on macOS arm64, with network disabled
and the source checkout mounted read-only. Linux native performance is not
measured. GoogleTest uses the existing pinned 1.18.0 source offline.

The native sanitizer binary instruments the configuration test file,
CRbeam.cpp, CmdLine.cpp, Configuration.cpp, MagneticConfiguration.cpp,
Sha256.cpp GridMagneticField.cpp and Background.cpp. It links the normally built core
library, GoogleTest and external dependencies. Flags are C++17, `-O1 -g`,
`-fno-omit-frame-pointer -fsanitize=address,undefined -fopenmp`, with
`ASAN_OPTIONS=detect_leaks=1` and `UBSAN_OPTIONS=halt_on_error=1`.
Its build metadata identifies a diagnostic binary, not a published release.
No whole-library sanitizer or ThreadSanitizer coverage is asserted.

## Coverage and limitations

Six new configuration cases cover CLI/YAML/snapshot equivalence for a grid,
absolute effective paths and load-time hashes, incompatible/duplicate/invalid
settings, grid failure before output creation, per-primary stochastic sums
and their seeded replay, and equivalence of explicit wave-model selectors
with the legacy CLI. They also parse both new examples. Existing tests cover
version 1 snapshots, legacy input, nuclear data, failure recording and hashes.
The nine grid library controls remain part of both platform checks.

Grid input is frozen in memory and those exact bytes are hashed and parsed
through the stream overload. Per-primary construction reuses the loaded grid.
Only configuration is written to config.yaml; file identity belongs in
run.yaml. No actual field files, field realizations or external source paths
are committed. Published examples require a user-supplied grid.

The checked C++ input sources have these SHA-256 values:

| Source | SHA-256 |
| --- | --- |
| MagneticConfiguration.cpp | `b61cdd6fdc38bbb26b2a554abea9146376e394286edecd3b7cdc14a58aa0bd77` |
| Configuration.cpp | `c3e4639f6ba2f8f0d0bc47b2edcd5d241246421292e1703cd849edcdbf9f7f44` |
| CRbeam.cpp | `a16bd378d8e3e4df098f5b41689e452f48cac879004a227e9f62fe8d738a9e0e` |
| CmdLine.cpp | `3370ede4d0bcd657acaf22ebbf230fc2743bcd6f7bfd1e1bc92e4b0031cdeec1` |
| Sha256.cpp | `7064bff016a428a05d792c9eea99a62f3293fcb1f9da3b9b44695d419b373829` |
| GridMagneticField.cpp | `ff1e6b85805387e2d756bca0e7178afbf99c2de42a428df4f97c3568b4ae14a6` |
| configuration_input.cpp | `221e742dd281386a3fc63d5aadd2e207698b154765e25a65389cd427bf2782e7` |

The first Linux sanitizer execution stopped in the pre-existing
MatrixBackground loader: an array allocated with new[] was freed with scalar
delete. Commit `1150a1a` corrects both releases to delete[], preserving the
background data and numerical model. The existing valid-additional-background
configuration case exercises this path. No sanitizer check is disabled to
bypass that failure. Background.cpp after the repair has SHA-256
`d7bded21ce882eeb7ba88fa0dc050bd1aa4eac6af92d71404301169f3cb362c3`.

Retained local artifacts are under untracked `build/grid-interface-verification`.
These short implementation checks do not validate physical maps, boundary
convergence or propagated distributions. Source geometry, observers, new
thinning algorithms, FFT turbulence and new backgrounds are separate stages.
