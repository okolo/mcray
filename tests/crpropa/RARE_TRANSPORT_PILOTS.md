# Independent full-Fe56 weighted-method pilots

These validation adapters propagate the complete Fe56 PD+decay genealogy,
including both photon sources and the eight lepton count/energy rewards.
They implement importance sampling and root-PD splitting without modifying
the nuclear tables, physical hazards, product kernels or kinematics. This is
the pilot stage of [v2 preparation](ACCEPTANCE_V2.md), not final acceptance.

The source revision must descend from the published ICS/pair-sampling fix
`54b24e0fa7e6aa55c89804b51afd7611f632945d`. That fix remains in the build;
EM transport is disabled in these isolated nuclear pilots. The original v1
protocol, 1%/5% tolerances, fixed-I normalizations and 99% final confidence
requirement remain unchanged. No pilot output is pooled with final samples.

## Sampling contract

Let p be the production conditional channel law and s the topology score.
Over channels with p>0, the sampling law is
`q = (1-eta)*p + eta*s/sum(s)`; if all active scores vanish, q=p.
Every selected event contributes its exact p/q factor. Eta is strictly below
one, preserving every physical path. Weights are not self-normalized or
clipped. The same topology-only guide is used for both programs; it contains
no deterministic reference yields. Scores favor channels producing, or
leading to, parents of the selected beta sign.

- mcray uses its original PropagationEngine with the original full PD and
  decay rates, waiting times and process selection. Only conditional channel
  sampling is wrapped; products come from the production PD/decay kernels.
- CRPropa uses its pinned production product kernels and reproduces its z=0
  PD-then-decay module order, within-step repetition and next-step limits.
  PD hazards retain linear grid interpolation and conversion from Mpc^-1
  to m^-1. Branches retain nearest-node cumulative selection, including the
  residual probability of the last branch. The original independent decay
  exponential race sets the minimum waiting time; its label is resampled
  with p/q. Decay hazards are checked against the official meanFreePath API.
- At the first A56 PD event, four alternative channel/kernel continuations
  replace the parent, each carrying one fourth of its incoming weight.
  A strictly decreases in the frozen PD data, so this split occurs at most
  once along a lineage. All daughter nuclei continue through both processes.
  Previously emitted products retain their existing weights.
- Continuations use subsequent non-reused random draws in the process-local
  RNG stream. RNG states are not cloned. Descendants share upstream history:
  all rewards and variance estimates aggregate inside the original primary.

The sampler estimates additive particle counts and energies. It does not
establish probabilities of coincidences or arbitrary nonlinear observables.
The output's energy column already includes the particle weight; photon bins
use the physical energy obtained by dividing this column by the weight.

## Fixed independent pilot matrix

Both backgrounds (CMB and CMB+Inoue12), both programs and all three frozen
numerical settings are included. Each cell has these six configurations:

| Profile | Guide | Eta | Root-PD copies |
|---|---|---:|---:|
| ordinary | unused | 0 | 1 |
| split | unused | 0 | 4 |
| minus-moderate | beta-minus | 0.5 | 4 |
| minus-strong | beta-minus | 0.8 | 4 |
| plus-moderate | beta-plus | 0.5 | 4 |
| plus-strong | beta-plus | 0.8 | 4 |

The ordinary profiles use the original native scenario and retained official
CRPropa probe/library. Weighted CRPropa has a separate new build at the same
pinned upstream commit; the retained official control is never overwritten.

`pilot`: eight independent batches of 500 primaries per cell, 576 jobs and
288000 primary executions. `smoke`: one batch of two primaries, main numerical
setting only, 24 jobs and 48 primary executions. Smoke and pilot seed
namespaces are disjoint, as are program/background/profile/batch seeds.
Only numerical refinements of a given batch are paired. Both namespaces
are separate from v1, v2 ordinary preparation and development controls.

Each execution first runs 48 Fe56 zero-bias identity-control primary
executions: both programs/backgrounds/settings, comparing original outputs
with the new PD sampling path on identical seeds. It also runs 1024 Mn54
finite-path method-control primaries (two mixture strengths, both programs).
Those controls exercise a rare beta-minus channel against a competing
beta-plus channel, using the actual rates, production kinematics and the
unchanged 51-node beta spectrum. Their fixed six-empirical-SE sanity check
is an engineering test, not a calibrated final 99% confidence assertion.
The final statistical method-check family remains to be frozen after pilots.

Total reservation: **289072 primary executions for pilot**, or 1120 for smoke,
including controls. Root split alternatives are counted separately: at most
720000 additional alternatives in the pilot matrix; actual counts are in
each output and report. Prior v1 reservations remain explicit in the ledger.
These sample counts are pilot budgets, not permission for unlimited reruns.

## Server preparation in /home2/ok

Reuse the existing input bundle, original prepared pilot references and
Python environment. NumPy 2.2.6, SciPy 1.15.3 and GitPython 3.2.0 are supported;
the build uses the existing C/C++/Fortran, Boost/GSL, yaml-cpp/OpenSSL and
CRPropa FFTW dependencies. Ubuntu 22.04 and CMake 3.22 remain supported.
No new generated nuclear data or deterministic-cache transfer is required.

The tools must first be committed and pushed to public mcray. Replace
`PUBLISHED_COMMIT` by the exact new commit containing these adapters. The
launcher resolves and records that revision, builds and runs all CTest
checks, clones CRPropa at its pinned commit, builds a separate weighted probe,
verifies input fingerprints, reserves the campaign and archives results.

Using an existing public checkout for the launcher:

```bash
source /home2/ok/transport-acceptance-venv/bin/activate
git -C /home2/ok/transport-acceptance-tools fetch origin talys
git -C /home2/ok/transport-acceptance-tools checkout --detach PUBLISHED_COMMIT

bash /home2/ok/transport-acceptance-tools/tests/crpropa/run_external_rare_transport_pilots.sh \
  /home2/ok/rate-inputs \
  /home2/ok/acceptance-references/20261016T102021Z-03eb6c57/transport-photons \
  /home2/ok/mcray-test-data \
  PUBLISHED_COMMIT 16 pilot
```

An optional first launch can use `smoke` instead of `pilot`; it has a separate
ledger and seed namespace. WORKERS is the fifth argument, 16 above.
OMP_NUM_THREADS, OPENBLAS_NUM_THREADS and MKL_NUM_THREADS are set to one in
each worker process. The long pilots must run on external native Linux amd64.

A launcher archive contains `tests/server/run_tests.py` and the external
wrapper. If using it instead of a checkout, verify its supplied checksum,
extract into `/home2/ok/rare-transport-launcher`, and use that wrapper path
with the same arguments. The archive bootstraps the launcher; the server
still fetches the exact published public code revision.

## Outputs, recovery and interpretation

The usual timestamped `.tar.gz` and `.json` are written under
`/home2/ok/mcray-test-data`. Keep the persistent
`rare-transport-v2-pilot-ledger.json` there; do not delete it to start another
campaign. The archive includes jobs/seeds, original-input and generated-guide
hashes, source/Python/dependency identities, all probe/library hashes,
control outputs, raw genealogical products and completed-batch markers.

`rare-transport-pilots/pilot-report.json` records batch means/SEs, per-primary
variance, nonzero-history counts, maximum primary contributions, baryon and
charge normalization diagnostics, clone counts and variance/cost ratios.
Zero variance is unresolved and never produces a claimed efficiency gain.
Recovered outputs have no usable original timing and are excluded from time
efficiency claims. No pass/fail physical acceptance is emitted.

Resume only from the retained exact source checkout and original output,
repeating the recorded workload command with `--resume --output ORIGINAL_OUTPUT`
and the same ledger. Complete marked batches are hash-checked and skipped.
A launcher start marker is written before spawning every batch. An unmarked
batch can be recovered only with matching start provenance and validated
complete histories/footer. A started job with no output is also unknown
consumption, because the native driver can fail before opening its final
output. Neither this case nor a partial file is automatically rerun. Incomplete control
directories stop and require an audit of consumed histories. Re-running the
external wrapper creates another checkout/output and is blocked by the ledger.
Any resumed results must be archived separately; preserve the original archive.

Choose final proposal parameters and budgets from variance/cost/tail diagnostics
and method checks, not from closeness of a pilot mean to its target. The
subsequent final seeds, full method-check inventory and confidence allocation
must be frozen before new acceptance observations. If waiting-time rarity or
weight tails remain limiting, keep that finding explicit; these pilots do
not authorize changing the physical hazards or acceptance tolerances.
