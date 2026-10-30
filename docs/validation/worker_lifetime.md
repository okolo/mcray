# OpenMP worker lifetime regression

The pre-existing worker cleanup released shared intrusive references from
several OpenMP threads at once. PropagationEngine serialized clone creation,
but each thread's SafePtr destroyed its clone outside that critical section.
Fields, filters and background objects can retain non-atomic intrusive counters;
concurrent const field queries do not make those lifetime mutations safe.

RunMultithread now explicitly clears the worker pointer under the same named
OpenMP critical section used for construction. Its automatic destructor then
sees an empty pointer. Propagation and the RNG operations remain outside this
cleanup section. This fixes the engine clone lifecycle without claiming that
arbitrary concurrent SmartPtr copies/releases elsewhere are safe.

## Regression evidence on 2026-10-30

The new GoogleTest case
`unit.concurrent.PropagationLifetime.WorkerCleanupSerializesSharedReferenceRelease`
creates two and four workers over three empty-stack batches each. Interaction
clones retain a shared intrusive token. An atomic observer records overlapping
clone destructors, with a short sleep widening the overlap window. The test
also checks clone/destructor counts, the restored token refcount after each
batch and exactly one final token destruction.

Before the fix the case failed on macOS because destructors overlapped for
both team sizes. After the fix it passed on macOS arm64 and emulated Ubuntu
22.04 Linux amd64, including 20 consecutive CTest repetitions on each platform.
The complete Release CTest suite passed 185/185 on both. Native Python checks
also passed 18 rotation, 9 loader/executable and 15 CI-publication guard cases.

These are focused structural/refcount and regression checks, not a claim of
general ThreadSanitizer coverage or new physical acceptance. The empty-stack
case deliberately isolates cleanup from physics and random history sampling.
Local logs and source snapshots are retained outside Git in
build/lifetime-ci-verification. The Linux C++ image is
`sha256:37682238eeaf9f01f7eddea17d4e55551d2e317307dfa5c8ea8dd0dc8084baf6`.

The new [automated build workflow](../usage/ci.md) runs the lifetime test with
the full short suite and the Python executable round trip on hosted Linux.
Workflow YAML passed actionlint and all shell steps passed bash syntax checks
locally. Hosted results are reported separately by the PR's build/test check;
this document records the local checks, not an assumed hosted success.
