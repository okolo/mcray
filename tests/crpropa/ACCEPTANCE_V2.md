# Transport/photon v2 preparation tools

Implementation status, 2026-10-24: ordinary batch preparation and rare-channel
weight/kernel foundations are implemented. **There is no final v2 acceptance
runner or fully frozen rare-method protocol yet.** None of these tools starts
the planned 20-hour acceptance workload. The immutable
[v1 protocol](protocols/transport_photons_v1.json) and its
[inconclusive results](../../docs/validation/transport_photon_acceptance_results.md)
remain unchanged.

## Ordinary batch preparation

`prepare_transport_photon_v2.py` verifies the completed v1 manifest, sidecar,
all payload hashes, job identities and usable unrecovered timing records.
It uses execution times only; it does not select batches using measured
acceptance differences. Physics, numerical controls, reference values and
photon intervals are inherited from the hash-checked v1 contract.

The preparer writes a new `ordinary-preparation.json` and checksum, with:

- 24 ordinary and eight rare observables, partitioning all original outputs;
- 128 independent main batches and 64 at each refinement, per program/background;
- fresh seeds, independent across backgrounds/programs and paired across settings;
- fixed histories per batch, sized before any v2 acceptance observations;
- 50%/25%/25% planned CPU allocation across numerical settings, equal allocation
  across backgrounds, and equal history counts for mcray and CRPropa within a cell;
- a default 16-worker, 20-hour resource estimate with a 1.15 timing safety factor;
- explicit pending rare budget, proposal, method-check inventory and final freeze.

Equal program CPU allocations would waste many more histories on the much
faster table-based reference adapter. Paired program sample sizes avoid that
imbalance. Timing estimates include per-process overhead from the previous
server; they do not guarantee performance on another host. Twenty hours is
not a statistical stopping rule or permission to discard unfinished batches.

Run from the public repository root, using Python >=3.10 and the existing
NumPy/SciPy/GitPython environment:

```bash
python3 tests/crpropa/prepare_transport_photon_v2.py \
  --previous-results /home2/ok/acceptance-v1/20261020T153324Z-d619dbcb \
  --previous-summary /home2/ok/mcray-test-data/20261020T153324Z-d619dbcb.json \
  --output /home2/ok/v2-preparation --hours 20 --workers 16
```

Use the extracted original result directory containing `manifest.json`;
preserve its entire payload. The output directory must not exist. This
command is safe to run locally because it only verifies data and plans jobs.

The external wrapper builds and tests an exact published revision, runs
preparation only and archives its output using the standard server harness:

```bash
bash tests/crpropa/run_external_transport_photon_v2_preparation.sh \
  /home2/ok/acceptance-v1/20261020T153324Z-d619dbcb \
  /home2/ok/mcray-test-data/20261020T153324Z-d619dbcb.json \
  /home2/ok/mcray-test-data PUBLISHED_COMMIT 16
```

Replace `PUBLISHED_COMMIT` only after the tools have been committed and pushed.
This preparation wrapper is not the future acceptance launch command.

## Rare-channel foundations

`rare_channel_sampling.py` provides a full-support mixture proposal, exact
conditional-channel p/q factors, splitting weights and aggregation of clone
rewards inside one independent primary history. It includes a constant-hazard
analytic fixture with a competing process and a finite path boundary. This
fixture tests the estimator; it is not a nuclear transport backend and does
not replace production rates or kinematics by a deterministic reference.

Holding times must retain their original law for a channel-only p/q factor
to be valid. Changing a decay rate, conditioning on an event before the end
of the path, or skipping competing processes needs additional likelihood
factors. The helpers do not authorize those changes. No self-normalization,
weight clipping or independent-clone variance estimates are implemented.

`NuclearDecay::GetChannelSecondaries` now invokes the same product kernel as
ordinary `GetSecondaries`, from the original complete dataset. Ordinary
channel selection and RNG ordering are preserved. Its caller is responsible
for the original branch probability, transport scheduler and weights.
The native probe exposes this as a conditional-only diagnostic:

```bash
build/Release/tests/nuclear_decay_probe \
  --data /absolute/path/to/original/nuclei --mode channel-events \
  --a 11 --charge 4 --channel 10010 --gamma 1e9 \
  --primaries 1024 --seed 3004000001 --output /absolute/path/to/be11-events.tsv
```

The output records original channel probability and total lab rate. These
are conditional event executions, not unconditional Fe56 histories. Record
their consumption in the development/pilot ledger. This probe currently
has decay photons disabled; the shared native kernel retains its existing
photon support. The pinned CRPropa public
`NuclearDecay::performInteraction(candidate, channel)` already supplies its
conditional kernel through the existing reference probe's `--mode events`.
Do not regenerate or overwrite the retained v1 reference binary.

## Remaining freeze gates

The full Fe56 weighted adapters still need to connect conditional kernels to
production-compatible competing-process scheduling, biased branch selection,
upstream-state splitting, independent RNG continuations and genealogical
rewards at every numerical setting. Rare-parent production and finite-path
survival must remain part of the measured unconditional yield. Native PD
channel-kernel access and CRPropa stepping require explicit implementation
checks; a correct forced decay alone does not close those gates.

Independent method pilots must determine proposal/splitting parameters,
stable batch uncertainty, rare sample counts and the finite method-check
inventory. Freeze these and all source/binary/data hashes before final
samples. Ordinary rare estimates remain diagnostic; do not choose an
estimator after examining which one agrees.

The preparation records failure-probability allocations 0.008 for all 384
original primary comparisons and 0.002 for predeclared statistical method
checks. Integral/spectral tolerances remain 1%/5%, with the original fixed-I
normalizations and agreement formula. Zero/constant or unstable estimates
remain inconclusive. Student intervals need pilot applicability checks;
unbiased weights alone do not guarantee interval coverage.

## Confirmed short checks

On macOS arm64, 2026-10-24: all 145 CTest checks pass after the native kernel
refactor; seven v1 acceptance tests and six new planning/weight tests pass.
The new native test compares seeded ordinary and conditional products and
energies exactly, preserves rates and weights, and rejects invalid channels.
The weight tests include exhaustive two-transition expectation identities and
finite-distance competing-hazard fixtures. Sixteen real-data conditional
Be11 beta-proton events produce the expected products using the original
8.3e-6 branch, without changing input tables. This is development evidence,
not rare-yield acceptance. Linux verification and rare-method pilots remain
pending. No long transport calculation was run locally.
