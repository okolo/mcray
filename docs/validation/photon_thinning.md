# Photon-threshold thinning verification

The implementation adds an opt-in PhotonThinning policy to modern main at
4f5a351f1cb0fd1927a6d7db05d0ff4752f9bc09. Its
[usage and numerical contract](../usage/photon_thinning.md) defines the
probability, inverse-probability weight, CLI/YAML fields and failure behavior.
The existing energy-based policy remains the default. Nonuniform jets, FFT
fields and new backgrounds are outside this change.

## Scope of the short checks

Eight new unit cases cover finite constructor validation in Release,
disabled/inapplicable RNG bypass, singleton selection, stable mixed-list
compaction, photon metadata and unaffected nuclear identity, invalid data,
probability/weight range failures before sampling, extreme unit scales,
frozen analytic moment controls and one-secondary engine failures in both
constant/variable-rate paths. Engine probes count actual interaction calls
to exclude retrying an unrepresentable history.

Four configuration cases cover explicit alpha>1 photon configuration,
CLI/YAML/config.yaml replay, invalid/inapplicable fields and schema-version
rules, legacy energy-policy equivalence, actual EM cascade reweighting,
and failed run metadata for an underflowing probability. The real-cascade
case uses four independent primaries; statistical checks are in the separate
frozen Bernoulli test, not inferred from this small interface fixture.
Two EM interface/regression cases have a 180-second wall-clock limit for
sanitized/emulated execution. Other input cases retain 30 seconds. An initial
emulated run reached the former 30-second limit in the cascade replay case;
its assertions and the statistical thresholds were not relaxed.

## Local results on 2026-10-30

- macOS arm64, Clang 22.1.5 Release: complete final suite 197/197 passed.
  The final input changes also passed all 25 configuration cases separately.
- Emulated Ubuntu 22.04 Linux amd64, GCC/GFortran 11.4.0 Release:
  complete suite 197/197 passed.
- The same Linux amd64 source built with ASan/UBSan and leak detection:
  complete suite 197/197 passed, without sanitizer reports or suppressions.

Linux uses the existing local image
`sha256:37682238eeaf9f01f7eddea17d4e55551d2e317307dfa5c8ea8dd0dc8084baf6`.
The local container has a read-only root with /tmp bound to a project-local
folder for Rosetta temporary files. An initial attempt without that writable
binding failed in the emulator during test discovery, not in application code.
The sanitizer build uses the same C/C++/Fortran instrumentation and runtime
settings as [CI](../usage/ci.md), including leak detection, fatal sanitizer
errors, frame pointers, non-PIE executables and no sanitizer suppressions.
System libraries are not rebuilt with instrumentation. This is not TSan
coverage or a new physical acceptance result.

Build logs, CTest JUnit files, source SHA-256 values, the timeout/emulator
diagnostics and the bounded local Linux command are retained outside Git in build/checks of the feature
worktree. These results are local; no hosted CI outcome for this new change
is claimed. New tests are discovered automatically by both required CI jobs.
