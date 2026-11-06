# Eight-hour EeV server workflow preflight

On 2026-11-06 the single-entry Ubuntu 22.04 launcher and the frozen EeV
streaming ensemble were prepared. This preflight consumes **no final seed**
and makes no statistical acceptance claim. Use the
[server instructions](../../tests/fork_compat/EEV_SERVER.md) and
[frozen protocol](../../tests/fork_compat/EEV_FIELD_ENSEMBLE.md).
[Machine-readable evidence](eev_field_server_preflight.json) records checked
source hashes, compiler environments, compressed-stream hashes and log hashes.

## Confirmed checks

- The new 18000-history fixture completed on native macOS arm64 Release,
  emulated Ubuntu 22.04 Linux amd64 Release and the corresponding Linux
  ASan+UBSan build with leak detection. Each used clean adf5900 and produced
  eighteen complete, fully decompressed/hash-verified gzip streams. No
  sanitizer diagnostics occurred. All used only seeds 2011101-2011118.
- The old photon runner's 6000-history fixture passed after extending its
  shared collection class with optional validator/runtime-prefix arguments.
  Default behavior and the old frozen confidence family are preserved.
- Twenty-eight Python guards pass on Linux. Native macOS passes twenty-seven
  and skips the Linux-only process-tree test. The Linux Python test runs in
  an existing arm64 container; that is separate from the amd64 numerical
  environment above. It is not an Ubuntu amd64 server benchmark.
- The real Linux watchdog test kills a worker that starts a separate process
  session when its parent workload times out. Other guards check exclusive
  final ledgers, traversal/link rejection, the fourteen-observable family,
  disjoint reservations, confidence widths, final-direction-only differences
  and invalid/incomplete/nonfinite summaries.
- The server workflow's phase sequencing and fixture-only isolation are tested
  with mocked external clone/configure/build/check commands. Real metadata
  packaging and status handling execute in that test. The launcher also
  verifies that its own bytes match the candidate commit before building.

The entire fresh Ubuntu 22.04 clone/build/CTest workflow was **not replayed
end-to-end in this session**. The available frozen amd64 compiler image lacks
Python/Git/curl. Instead, its actual compiler/core/sanitizer execute the new
numerical runner, and separate Linux Python/phase tests exercise orchestration.
The full short CTest suite and Python guards are mandatory automated phases
of the delivered server launcher before final collection. Production code
was unchanged; no prior large ensemble or new final reservation was run.

## Retained outputs and transfer

Raw fixture archives, summaries and logs remain under ignored
build/checks/eev-ensemble-verified-{native,linux-release,linux-sanitizers}; the
legacy regression is build/checks/eev-server-legacy-fixture. These fixture
reports intentionally contain no acceptance intervals. Candidate source
hashes and tested revisions in the evidence distinguish subsequent launcher/
documentation guards from the unchanged compiled numerical code.

The prepared transfer package contains the standalone launcher, a Git bundle,
REVISION, the pinned reference/GoogleTest archives, SHA256SUMS and the server
instructions. It supports unpushed committed code without this private
workspace or a remote branch update. The full eight-hour test must be started
on the external server; its final reservations remain unconsumed locally.
