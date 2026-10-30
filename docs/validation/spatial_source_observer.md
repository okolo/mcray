# Spatial source, sphere observer and output checks

Checks on 2026-10-30 apply to the local source/observer/output extension on
`codex/grid-magnetic-fields`, based on `1c8229a9b442164aa4e9ea7acfb791280eeb7ac7`.
These results cover implementation and regressions, not physical agreement
with the historical fork or an independent propagation code.

| Check | Result |
| --- | --- |
| macOS arm64, Clang 22.1.5 Release build | Passed |
| Complete native short CTest suite | 184/184 passed |
| Final-source native configuration/geometry/output cases | 30/30 passed |
| Published example configuration and delay command | Passed locally |
| Ubuntu 22.04 Linux amd64 C++ ASan/UBSan build and targeted cases | 30/30 passed |

The Linux run uses GCC 11.4.0 under amd64 emulation on macOS, with network
disabled and a 1200-second limit. Its pinned local image is
`sha256:37682238eeaf9f01f7eddea17d4e55551d2e317307dfa5c8ea8dd0dc8084baf6`.
The current C++ library, configuration library and tests were rebuilt with
AddressSanitizer and UndefinedBehaviorSanitizer; Fortran was not instrumented.
Leak detection was enabled and UBSan configured to stop on errors. All 30
targeted tests passed without sanitizer reports. This is an emulated Linux
check, not a native-server performance measurement.

The targeted cases cover finite/nonzero axes, coordinate units, normalized
CLI/YAML settings and malformed vectors/regions/angles/end times; cone moments
at a fixed seed; preservation of the pencil RNG stream and particle identity;
strict membership at birth and on the surface; sphere/end-time step caps;
photon and Fe56 propagation in both interaction-rate paths; and birth/exit
registration through the real application energy filter. End-redshift-only
runs stop at the requested time. Cone CLI/YAML/saved-config runs have identical
seeded output within a platform. The existing grid-only replay and legacy
input checks remain in the full native suite.

Output tests open several files, including two isotope groups, and flush
again after another file selects a different scale. They verify that each
file retains its own scale. Zero-flight registration uses the agreed fallback,
retains finite coordinates and records zero actual flight; a later positive
flight does not change the scale. Tiny positive flight retains the original
choice and a non-exponential header. The optional fork columns retain their
positions and FlightMpc is last.

Local diagnostic artifacts are retained outside Git under
`build/spatial-verification`: build/test logs, source-manifest.json, the Linux
check script, the example's config/run/table files, and calculated delays.
The manifest SHA-256 is
`ce3f702140e99aee8e94e69b8fba02795bde16bfff5dacc6cb3e82cee9531386`.

The sphere algorithm preserves sampled endpoint membership, including its
small step floor; it is not an exact surface crossing algorithm. Normalization
and text precision do not remove subtraction cancellation in 1-z/T. These
checks do not change or validate the engine's existing cosmological
coordinate/time convention, nor establish a cosmological arrival-time
conversion for the geometric excess-time tools. Ordinary output remains the
existing modern schema; fork-compatible output is selected by nondefault
spatial settings.
