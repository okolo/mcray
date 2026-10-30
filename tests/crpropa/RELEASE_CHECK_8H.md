# Eight-hour limited-release check

On 2026-10-26 the user selected one eight-hour ordinary check followed by
release with documented unresolved validation limits. This replaces the
planned 20-hour resource allocation for this release; it does not declare
the rare-method freeze complete or change physical tolerances.

The check reuses the frozen initial Fe56, 1 Mpc, E/A=10^18.6 eV case,
CMB and CMB+Inoue12, both programs, both photon sources and all three existing
numerical settings. It uses original ordinary production transport only.
The ICS/pair fix remains an ancestor; EM propagation, GZK and PPP are off in
this isolated PD/decay comparison. No production model or table is changed.

## Fixed work and time limit

The previous complete v1 payload is hash-verified and used for timing only.
Before new histories, the planner freezes 1024 jobs: 128 main batches and
64 at each refinement per background/program. Allocation is 50%/25%/25%
across settings; program sample sizes match within a cell. A 7.5-hour planning
budget with a 1.25 safety factor estimates approximately six hours of worker
time, leaving preparation, reporting and runtime slack. On the original
16-worker timing records this reserves 7935872 primary executions.

Seeds start at 3300000000, separate from all prior campaigns. Only resolutions
of the same primary batch are paired. All old reservations remain spent;
no previous sample is pooled. Histories are fixed before outcomes, not
extended until a comparison agrees. A persistent ledger and launcher lock
prevent duplicate campaigns; there is no automatic replay or resume.

The eight-hour wall limit includes restoration, clone/build/test and the
workload. The workload leaves three minutes for reporting, and the outer
runner terminates its whole active process group at the deadline. Logs and
completed/partial raw files are archived even on timeout. Compression runs
after the calculation limit and can add time. An incomplete campaign remains
incomplete; completed subsets are not a replacement acceptance ensemble.

## Single-command server launch

Use the exact published commit containing these files. Activate the existing
`/home2/ok/transport-acceptance-venv` environment and run:

```bash
bash tests/crpropa/run_external_release_check.sh /home2/ok PUBLISHED_COMMIT 16
```

The wrapper locates `/home2/ok/rate-inputs`, the retained pilot reference,
the 20261019 coupled deterministic cache and the 20261020 v1 results. If
retained work directories were removed, it restores the original archives
from `/home2/ok/mcray-test-data` into `/home2/ok/release-check-inputs`.
Keep `20261020T153324Z-d619dbcb.json` beside its original archive.
No new nuclear data is required. System `flock` and
the existing Linux amd64 compiler/scientific environment are used.

The external harness fetches the exact published revision into its own fresh
checkout. The persistent ledger is
`/home2/ok/mcray-test-data/release-check-8h-ledger.json`. Do not remove it or
rerun after a failure. Main outputs are in the usual timestamped archive
and JSON sidecar. The report is `release-check-8h/release-check-report.json`.

## Interpretation and release policy

All original 32 outputs remain recorded. The 24 ordinary outputs retain
1% integral and 5% spectral tolerances, fixed independent normalizations,
the original native/official/independent reference roles and a nominal
0.008 primary failure-probability allocation over the original 384 rows.
The eight rare outputs are explicitly diagnostic; their applicability and
weight-tail problems are not certified by this ordinary run. A rare row
cannot establish agreement, even when an exploratory interval is narrow.
The unresolved 0.002 method-check allocation is not claimed to have passed.
Paired step diagnostics are retained without correcting empirical means.

`passed` means the workflow completed; inspect the report for actual
agreement, disagreement, incomplete work and uncertainty. Release under the
user's resource-limit decision is separate from a claim of complete 99%
physical acceptance. Keep unresolved rare statistics, reference-step effects,
the beta-CDF approximation and photon energy accounting explicit in release
notes. This tool neither tags nor publishes the release automatically.

## Preparation verification

On 2026-10-26, all 35 CRPropa-tool Python tests and seven server-runner
tests pass locally, including actual descendant cleanup on timeout and
archive creation after deadline expiry. Shell syntax checks pass. A plan-only
run over the fully verified original server payload confirms 1024 jobs and
7935872 reserved primary executions, with a 5.9997-hour runtime estimate.
No new nuclear histories or long calculations were run on the workstation.
