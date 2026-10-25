# Full weighted nuclear transport: development preflight

Date: 2026-10-24. Implementation continues on talys above sampling-fix commit
`54b24e0fa7e6aa55c89804b51afd7611f632945d`. The ICS/pair correction is preserved;
the isolated PD/decay scenarios leave EM transport off. Exact development
source/binary fingerprints and raw results are retained outside Git. These
checks prepare independent external method pilots, not final acceptance.

## Implemented scope

- Native full Fe56 genealogy in the original PropagationEngine, with original
  PD full rates, exclusive-channel probabilities, decay rates, process timing,
  production conditional kernels and both photon sources.
- Pinned-CRPropa full genealogy with its PD-then-decay step convention,
  original product kernels, nearest-node branch selection and original
  interpolation/hazard laws. The official unweighted control is retained
  separately with its original binary/library fingerprints.
- Defensive topology-density proposals with exact p/q, first-A56-PD splitting,
  prefix weights for previously emitted rewards and primary-level variance
  accounting. No physical rates or tables are tilted, weights clipped,
  reference means inserted, or photon energy accounting changed.
- Fixed independent pilot planning, fresh seed namespaces, persistent history
  reservation, complete-output recovery guards and archived diagnostics.

The [public pilot contract](../../tests/crpropa/RARE_TRANSPORT_PILOTS.md)
contains the runnable settings and /home2/ok commands.

## Confirmed checks

On native macOS arm64:

- All 155 CTest entries pass, including the merged EM-spectrum tests, a new
  seeded PD conditional-kernel comparison and three weight-law controls.
- Eighteen Python checks pass: seven v1 acceptance, two v2 preparation,
  four weight fixtures and five pilot planning/aggregation/recovery checks.
- Forced use of the new zero-bias PD path reproduces original particle
  identities, weights, energies and order exactly on identical seeds for
  both programs/backgrounds and all three numerical settings: 12 paired
  controls, 48 primary executions. This exercises the new PD table/rate
  handling; it does not merely delegate to the original process method.
- The complete latest smoke workflow finishes all 24 jobs, 48 Fe56 primary
  executions and both sets of automatic preflight controls. Actual root
  split alternatives in that smoke: 51. No unmarked pilot outputs remain.
- The 1024 finite-path Mn54 controls exercise both programs at eta=0.5 and
  eta=0.8, including the rare beta-minus branch and its competing beta-plus
  branch. All eight lepton count/energy first moments pass the fixed
  engineering sanity check against the original-rate/51-node-sampler means.

| Program | Eta | Largest absolute difference / empirical SE over eight moments |
|---|---:|---:|
| mcray | 0.5 | 0.753 |
| mcray | 0.8 | 1.195 |
| CRPropa | 0.5 | 0.883 |
| CRPropa | 0.8 | 0.873 |

These SE comparisons are method-development diagnostics. They neither
replace the 1%/5% physical tolerances nor constitute the final joint 99%
statistical method-check family. The smoke has insufficient batches and rare
events to estimate production efficiency or establish physical acceptance.

The external full plan contains 576 jobs, 288000 Fe56 primary executions,
48 identity-control executions and 1024 finite-decay-control executions.
The long pilot matrix was not run on the workstation. Native Linux amd64
execution, weight-tail diagnostics, useful variance reduction at the original
Fe56 point and the final rare-method freeze remain pending.

Publication preflight on 2026-10-25 repeated all 155 native CTest entries
and the complete `tests/crpropa/test_*.py` discovery in the existing scientific
Python environment: all 32 Python tests passed. The system Python lacks SciPy;
use the documented validation environment for these checks.

## Development findings and retained evidence

An initial reference-adapter prototype omitted division by Mpc when loading
the PD rate file. Inspection of the root split counts exposed the unit error;
the production upstream loader was checked and its per-node conversion was
restored. Zero-bias identity controls were strengthened to force execution
through the new PD path. The invalid prototype outputs are retained and
excluded from validated evidence and all physical estimates.

Multiplying tiny physical branching fractions by a modest relative topology
guide did not resolve the Mn54 rare branch in a short development control.
The fixed pilot implementation instead mixes the physical law with a
separate topology proposal over positive physical support. Exact p/q keeps
the same target expectation. The source model and acceptance tolerances
were not changed. A provenance-field delimiter in the original native output
was also accommodated by the parser, without changing any output histories.

All development attempts, including rejected prototypes, remain outside Git;
their complete-output counts and uncertainty about any incomplete records
are retained separately from the server pilot reservation. Repeated fixed-seed
development regressions are not independent acceptance observations. Pilot
parameters must be frozen before final independent acceptance samples, after
the external variance/cost and batch-applicability analysis.
