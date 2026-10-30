# Fork integration branch review

On 2026-10-30, origin/main was fetched at
`2dd26db5444458c06872baf50ce630bb2eb7a1ac`. Its two commits since the feature
branch's base updated contributor/documentation guidance and CI review guards;
they did not modify propagation code. Merge commit `e36dae9` incorporates them
into codex/grid-magnetic-fields without conflicts or history rewriting.

## Scope and corrected findings

Local review covered field loading/composition, source CLI/YAML normalization,
observer stop/step/output consistency, initial membership, file normalization,
ordinary-run compatibility, and Python postprocessing. It found these issues:

1. **Legacy modern continuous-source coordinates.** RawOutput3D normalizes
   each row by its own flight interval but puts only the first interval in
   the header. Applying one header scale silently reconstructs wrong positions
   for different primary birth times. The loader now rejects this layout when
   config.yaml/params declares continuous emission or z_src varies. The new
   spatial format uses an immutable file scale and is unaffected; historical
   fork tables have their distinct z_production layout.
2. **Delay command output precision.** printf-style %.17g formatting narrowed
   extended scalars on Linux amd64. The next value above longdouble(1), printed
   as 1.0000000000000000001 with scalar formatting, became 1. Direct NumPy scalar
   formatting now preserves it; a subprocess round-trip regression covers it.
3. **Known incomplete results.** The loader ignored adjacent run.yaml status.
   It now rejects failed, interrupted, running, unknown/missing statuses and
   malformed YAML roots, while allowing completed runs. Historical standalone
   files without run metadata remain supported; absence of metadata is not
   proof of a completed run.

The executable round-trip fixture now controls threads with OMP_NUM_THREADS
and asserts the recorded count. It previously passed an unsupported --threads
option, which the compatible CLI ignored. Tests keep temporary data inside
the project. No propagation or accepted physical model changed in these fixes.

## Final checks

| Check | macOS arm64 | Linux amd64, emulated |
| --- | --- | --- |
| Release CMake build and complete short CTest suite | 184/184 passed | 184/184 passed |
| mc_rotate numerical/block checks | 18/18 passed | 18/18 passed |
| Spatial loader/delay checks | 9/9 passed, including executable round trip | 8 passed, 1 executable case skipped |
| CI review publication-guard unit tests | 15/15 passed | 15/15 passed |

The Linux C++ suite includes real CRbeam CLI/YAML runs for cone replay,
birth membership, sphere exit and end-time capping. Its Python container
lacks the compatible Linux CRbeam binary, so the Python executable fixture
was explicitly skipped there; it was exercised on macOS, not reported as a
Linux Python pass. The new metadata, normalization and precision regressions
ran on both platforms. CI unit tests mock publication; they do not send GitHub
reviews or establish that the hosted workflow ran.

The C++ merge build used Clang 22.1.5 on macOS and GCC 11.4.0 / Ubuntu 22.04
under amd64 emulation for Linux. C++ sources are unchanged by the Python review
fixes; earlier sanitizer scope remains in the
[spatial implementation results](spatial_source_observer.md). The new check
is a normal Release build, not a new sanitizer run.

Linux C++ image:
`sha256:37682238eeaf9f01f7eddea17d4e55551d2e317307dfa5c8ea8dd0dc8084baf6`.
Linux Python image:
`sha256:8d38c4e97825896aa9f2ae1d4bc66f2c5e5131925bb39fd3c544136fb8ba2f74`,
Python 3.12.3, NumPy 1.26.4 and PyYAML 6.0.3. The PyYAML cp312 manylinux x86_64
wheel SHA-256 is
`ba1cc08a7ccde2d2ec775841541641e4548226580ab850948cbfda66a1befcdc`;
the NumPy wheel is recorded in the [rotation checks](mc_rotate.md).
Containers had read-only roots, no network, project-local build/temp/cache
paths, and a 1200-second bound on the full Linux build/test script.

Commands, logs, downloaded dependency artifacts, source hashes and the
precision reproduction are retained outside Git under build/premerge-review.
Build artifacts are not part of the branch. No push or PR publication was
performed as part of this local integration check.

## Remaining limits

The observer retains the agreed sampled membership semantics, not an exact
crossing solver. The geometry/time convention and physical acceptance limits
of existing propagation remain unchanged. Header-only legacy data cannot
recover missing individual scales or prove the run configuration; the loader
rejects identifiable ambiguous modern continuous-source cases rather than
inventing coordinates. Historical targeting/Healpy scripts and later FFT,
nonuniform-jet, thinning and background extensions are outside this stage.
A formal PR review and hosted CI remain subsequent publication steps.
