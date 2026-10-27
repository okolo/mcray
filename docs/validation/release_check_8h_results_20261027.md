# Eight-hour release check: retained results, 2026-10-27

The user selected release after one resource-limited ordinary transport check,
with unresolved validation disclosed and no additional research iterations.
This report records that check without extending or rerunning its samples.

## Source and integrity

- Run: `20261026T210541Z-03385791` on native Linux amd64, 16 worker processes.
- Tested source: `18b28b51ed24945f3a9638fb55b804897159f1b6` on talys.
- Archive SHA-256:
  `6d892bd176de81d8e0c367050c76aa0bbf9ecc15235039e1b3f36e54c7e6a892`.
- All 4123 manifest-listed payload hashes and the archive checksum match.
- Build succeeded and all 155 CTest entries passed in 284.37 seconds.
- All ten current-production coefficient controls matched the frozen results.
- Original physical protocol, input manifests, retained official runtime and
  the effective eight-hour protocol remain identified by their recorded hashes.

The [check contract](../../tests/crpropa/RELEASE_CHECK_8H.md) fixes the scope:
initial Fe56, E/A=10^18.6 eV, 1 Mpc, CMB and CMB+Inoue12, both programs,
photodisintegration plus decays and both nuclear photon sources. This is
ordinary production sampling, without importance sampling or splitting.
GZK, PPP, EM transport, cosmological evolution and free-neutron decay are
outside this isolated comparison. Their software support and previous checks
must not be confused with this ensemble's scope.

## Deadline and incomplete work

The outer harness records `time_budget_exhausted`, exit code 124. Its logged
interval is 2026-10-26 21:05:41 UTC through 2026-10-27 05:05:40 UTC: the agreed
eight-hour calculation limit. Compression follows that limit.

Of 1024 fixed jobs, 1017 completed, supplying 7860369 complete primary
histories from a reservation of 7935872. All refinement groups and both
CRPropa main groups finished. The mcray main groups finished 125/128 batches
on CMB and 124/128 on CMB+Inoue12. The seven missing jobs reserve 75503
histories. Their files contain provenance headers but no completed particle
ensemble/footer. Actual consumed histories are unknown; those reservations
remain spent and must not be inferred to be zero.

The application was interrupted by the controller deadline; this is not
reported as a completed ensemble or an application-physics exception. The
outer deadline also stopped final report generation. Neither original
`acceptance-report.json` nor `release-check-report.json` exists in the archive.
This absence is a reporting limitation of the bounded workflow, not missing
raw results for the 1017 completed jobs.

An offline audit verified every completed job's start/output hashes and
recomputed all stored statistics from its raw TSV. The original parser checks
complete genealogy, unit weights, baryon number 56 and charge 26 per primary.
No complete-history normalization failure or cached-summary mismatch was
found. The audit does not recover, replay or count any incomplete history.

Reproduce after extracting the original archive into a retained directory:

```bash
python3 tests/crpropa/audit_release_check.py \
  --results /path/to/20261026T210541Z-03385791 \
  --summary /path/to/20261026T210541Z-03385791.json \
  --archive /path/to/20261026T210541Z-03385791.tar.gz \
  --output /path/to/new-offline-audit
```

The output is a separately identified `release-check-audit.json`, not a
replacement file inserted into the original archive. No transport is run.

## Ordinary comparisons

The original 32 outputs comprise 24 ordinary and eight rare quantities. With
three resolutions, two backgrounds and two programs there are 288 ordinary
and 96 rare rows. The audit retains the fixed independent normalizations,
1% integral / 5% spectral tolerances and nominal 0.008 confidence allocation
over the original 384 rows. It does not widen the tolerances or retune means.

Of 240 ordinary comparisons in fully completed groups, 108 meet the existing
agreement rule and 132 remain inconclusive; none meets its disagreement rule.
The remaining 48 ordinary comparisons belong to the two incomplete native
main groups and remain inconclusive regardless of their nominal intervals.
For reference only, 34 of those 48 nominal intervals would lie inside their
bands. Runtime censoring can depend on histories, so those descriptive
intervals are not promoted to frozen-ensemble acceptance.

Main-level photon-count estimates illustrate the achieved bulk precision:

| Background | Program | Complete primaries | Photons / primary | Own deterministic coefficient reference | Relative mean difference using fixed independent normalization |
|---|---|---:|---:|---:|---:|
| CMB | mcray | 1558625 | 0.41821221 | 0.41834338 | -0.0314% |
| CMB | CRPropa | 1596032 | 0.43302453 | 0.43275482 | +0.0645% |
| CMB+Inoue12 | mcray | 1180976 | 0.48895744 | 0.48957666 | -0.1265% |
| CMB+Inoue12 | CRPropa | 1219072 | 0.50384555 | 0.50370806 | +0.0281% |

These compare each Monte Carlo implementation with its own deterministic
coefficient representation, not mcray directly with CRPropa. The native main
numbers remain descriptive because their full fixed batch inventories did
not finish. Photon-energy intervals on the combined background remain
inconclusive even where central estimates are close to their references.

## Numerical resolution diagnostics

All 64 predeclared main/refinement batch pairs are available for both
backgrounds and programs. Pairing uses the same seeds; sample sizes differ
between settings according to the fixed resource plan. There is no empirical
mean correction or pooling refinements as independent replicas.

Main-to-finest photon differences, relative to the fixed independent
normalization, are:

| Program/background | Count difference +/- paired empirical SE | Energy difference +/- paired empirical SE |
|---|---:|---:|
| mcray / CMB | +0.0193% +/- 0.3740% | -0.0588% +/- 0.4380% |
| mcray / CMB+Inoue12 | -0.0995% +/- 0.3765% | +0.2156% +/- 0.4807% |
| CRPropa / CMB | -0.8350% +/- 0.5100% | -0.6891% +/- 0.6723% |
| CRPropa / CMB+Inoue12 | -0.0117% +/- 0.5221% | +0.6550% +/- 0.6587% |

The earlier 4000-history pilot's approximately 11% normalized CMB photon
shift is not reproduced at that size in this larger check. The intervals
above still do not establish a simultaneous 1% step-error bound. Record the
improved diagnostic evidence without declaring all numerical gates closed.

## Rare outputs and release interpretation

Main-level electron counts are only 1 and 4 for mcray, and 8 and 8 for
CRPropa on CMB and CMB+Inoue12 respectively. Main positron counts are
5210/6043 for mcray and 5802/6875 for CRPropa, but simultaneous 1% count/energy
acceptance remains unresolved. All 96 rare rows retain diagnostic status.

Three exploratory antineutrino-energy intervals nominally meet the old
disagreement formula: native CMB main, and both native CMB+Inoue12 refinements.
Each is based on a single observed antineutrino. Empirical Student intervals
in these extremely sparse samples lack established coverage. These signals
are disclosed, neither suppressed as successful acceptance nor promoted to
validated model rejection. No rare estimator or interval is selected after
looking for agreement.

Release proceeds under the user's resource-limit decision with these explicit
limits. The campaign is incomplete; it does not establish full simultaneous
99% physical acceptance, and it is not extended to obtain that claim.
See [release notes](nuclear_release_notes.md) for the supported functionality,
metadata migration, data requirements and retained physical approximations.
