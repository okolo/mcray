# Bounded rejection: reproducible external server test

Python 3.10 is supported on Ubuntu 22.04; upgrading Python is unnecessary.
The portable package contains a Git bundle, this adapter and pinned Python
requirements. NR is fetched from blackstonep/Numerical-Recipes at the original
pinned revision, and all four headers are SHA-256 checked. An optional
`--nr /path/to/nr` uses matching local headers instead.

From the published branch on Ubuntu 22.04 / Python 3.10:

```bash
sudo apt-get update
sudo apt-get install -y git cmake make gcc g++ gfortran libgsl-dev \
  libboost-filesystem-dev libboost-system-dev python3.10-venv ca-certificates

git clone --branch codex/bounded-rejection-sampling \
  https://github.com/okolo/mcray.git mcray-bounded-rejection
cd mcray-bounded-rejection
python3.10 -m venv .venv
.venv/bin/python -m pip install -r tests/sampling/secondary-requirements-py310.txt

# Freeze the revision of these checked-out tools for the server experiment.
revision=$(git rev-parse HEAD)
.venv/bin/python tests/sampling/run_secondary_server.py \
  --server --repo https://github.com/okolo/mcray.git --ref "$revision" \
  --baseline both --hours 24 --jobs 4 --raw-gib 20 \
  --data-dir /home3/ok/mcray-test-data
```

Use a fresh checkout directory, or update an existing clean checkout to the
published branch before freezing `revision`. No recursive submodule checkout
is needed for this test. NR is fetched and checked by the runner. Installing
system dependencies and the Python environment happens before the experiment's
24-hour clock. For an unattended run, execute the final command inside tmux
(or an equivalent persistent server session).

From the unpacked package:

```bash
PYTHON=python3.10 bash sampling/run_fixed_server.sh \
  --baseline both --data-dir /home3/ok/mcray-test-data \
  --hours 24 --jobs 4 --raw-gib 20
```

The default `--baseline both` samples three independent streams: frozen
rejection from `1cc28f3`, the unchanged numerical CDF implementation from
`fa562d2`, and the new bounded method. Conditional secondary laws are shared.
Single-control alternatives are `--baseline rejection`, `--baseline cdf`,
and `--baseline nr` (the original NR invariant sampler).

The 24-hour budget begins when the Python experiment starts, after the wrapper
creates its venv and installs requirements. It includes source/NR fetch,
build, unit tests, deterministic controls, reference integration, pilots,
sampling, analysis, hashing and compression. One hour is reserved for final
work. Build parallelism is controlled by jobs; Monte Carlo streams are run
sequentially. Run this full experiment on the external server.

The original eleven ICS/GammaPP cases, nonuniform energy bins, cumulative
tails, per-interaction multiplicity statistics, and seed namespace are retained.
Pilots time all three methods; N and every main seed/chunk are frozen before
main sampling. The raw-data plan accounts for all three streams. Insufficient
statistics, zero observed tail events, incomplete chunks and smoke runs cannot
establish equivalence. Counts are never used to stop or merge bins.

Precision targets remain **1% in core bins and 10% in tails**, with global
significance **1%**. Bounded-versus-rejection and bounded-versus-CDF each
receive alpha=0.5%. Each uses the original exact binomial/Fisher tests, Holm
adjustment and simultaneous exact confidence intervals. Bonferroni across
the two families controls the joint level even though they share the bounded
sample. A single-control run retains alpha=1%. The inference implementation
in secondary_statistics.py is unchanged.

The aggregate report is report.json/report.md. Detailed bins, confidence
intervals and tests are in report-old and report-cdf (JSON, TSV and Markdown).
Reanalysis uses the saved frozen input:

```bash
python3 sampling/analyze_bounded_secondary.py RESULT/analysis-input.json RESULT/reanalysis.json
```

An independent optional audit reconstructs all saved multiplicity bins directly
from binary events, checks types/kinematic support/energy and verifies hashes:

```bash
python3 sampling/audit_secondary_raw.py RESULT /tmp/raw-audit.json --seconds 600
```

Deterministic controls audit global/local envelopes throughout the support,
table-knot neighborhoods, and tail proposals corresponding to probabilities
down to 1e-12. CDF quantiles are checked against independent refined integration.
Timing controls include full per-event construction for all methods. The raw
bounded stream records proposal counts, maximum cell count/gap, and fallback
count; unexpected EM fallback or an envelope violation aborts with diagnostics.

The retained directory and tar.gz include exact source/history, resolved
commit, adapter, executable, build configuration, dependency versions, pinned
NR provenance, reference convergence, frozen plan/seeds, all raw events,
per-bin counts, timings, logs, reports and SHA-256 hashes. Private NR headers
remain in the user's result archive and are not distributed in Git/package.
Deadline and failure paths retain partial results and attempt compression.

Exit codes: 0 within target; 2 material deviation; 3 inconclusive; 124 deadline;
1 execution error. For an infrastructure check only:

```bash
PYTHON=python3.10 bash sampling/run_fixed_server.sh \
  --smoke --hours 0.25 --raw-gib 0.01 \
  --data-dir /home3/ok/mcray-test-data
```

Smoke exit code 3 is expected. A full 24-hour acceptance result remains a
separate server measurement, regardless of successful local checks.
