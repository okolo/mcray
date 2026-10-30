# Historical charged-particle transport diagnostic

On 2026-10-31, the approved deterministic transport-first stage compared the
actual historical deflection functions against modern code. This is a
numerical diagnostic, not Monte Carlo or physical-equivalence acceptance.
The initial baseline below changed no production physics or integration
algorithm. The subsequently approved scale-rounding correction and its
unchanged-matrix results are recorded at the end of this document.

## Protocol and provenance

[The runnable protocol](../../tests/fork_compat/TRANSPORT.md) specifies the
756-case matrix, extraction, assumptions, outputs and reproduction command.
It covers electrons, positrons and free protons; free motion; axial/oblique
constant fields; mildly relativistic constant-field motion; and a nonuniform
complete grid with six boundary/inside/outside starts and four extension modes.
Accuracy 10/20/40 is crossed with 1/4/16 equal outer calls over 0.25 Mpc of
elapsed c*time. All combinations are retained; no ensemble is sampled.

Both paths use the explicitly agreed H0=70 km/s/Mpc and OmegaLambda=0.7,
constant energy and no interactions or observer. The audited CosmoTime class
and unit implementation are shared after source-identity checks; cosmology
itself and the propagation engine's outer-step scheduling are not under test.
The original double beta formula and particle mass/charge functions are
extracted with the historical transport, not replaced with modern arithmetic.

Production source, adapter, twelve historical inputs, generated data and
linked-library hashes are retained in the
[machine-readable evidence](fork_transport_comparison.json). The executed
implementation is `8de5cce`; all three executions recorded a clean working
tree. Generated complete raw trajectories, extracted source, fragment hashes,
manifests and compiler/execution logs remain in ignored `build/checks/transport-*`
directories. Reproduction accepts an explicit read-only fork path and a current
compatible modern build; it does not depend on an enclosing private workspace.

## Execution checks

| Environment | Completed cases | Maximum old/new position difference, Mpc | Maximum direction-vector difference |
| --- | ---: | ---: | ---: |
| macOS arm64, Clang 22.1.5, Release | 756 | 1.829804718473232e-4 | 1.108845366742202e-3 |
| Ubuntu 22.04 Linux amd64, GCC 11.4.0, Release | 756 | 1.829804718473953e-4 | 1.108845366742162e-3 |
| Same Linux toolchain with ASan+UBSan and leak detection | 756 | 1.829804718473953e-4 | 1.108845366742162e-3 |

Linux executions used amd64 emulation on the arm64 workstation, with the
frozen offline image recorded in the JSON evidence. They are not native Linux
performance measurements. Sanitizer flags from the configured instrumented
build applied to the driver, extracted reference and explicitly compiled
modern kernels too; runtime options enabled leak detection and fail-fast
address/undefined diagnostics, with no suppression. The runs completed without
sanitizer diagnostics. Compilation and execution were bounded below 30 minutes.

Six self-contained Python harness tests passed. They cover frozen case
coverage, deterministic inputs, reference-hash rejection, duplicate/nonfinite
report rejection, preservation of all refinement groups, and the unchanged
default compiler-helper path. The existing field/cone probe still compiles
through that default path; a field-only smoke query ran successfully, with a
2.07e-25 G component difference. No cone ensemble was rerun for this stage.
Before each transport run, independent analytic checks cover free motion
and quarter-turns for both charge signs.

## Free and constant-field controls

All 108 free/constant-field cases agree to rounding. On Linux, maximum
old/new position differences are 2.14e-19 Mpc (free), 1.49e-19 Mpc
(ultrarelativistic constant field), and 2.64e-18 Mpc (gamma=10 control).
Directions agree exactly in these controls. The mildly relativistic beta
values differ by at most 1.05e-17; each analytic position comparison uses its
implementation's own beta, separating that effect from integration error.

Both propagators have the same constant-field integration errors. Maxima
over constant-field cases and outer-call counts on Linux are:

| Accuracy | Position error versus analytic helix, Mpc | Direction-vector error versus analytic helix |
| ---: | ---: | ---: |
| 10 | 6.32754e-5 | 5.13515e-4 |
| 20 | 3.18485e-5 | 2.56741e-4 |
| 40 | 1.59770e-5 | 1.28366e-4 |

The approximately first-order error decrease is a property of the retained
normalized direction update; it is not evidence for exact helical propagation.
Maximum direction-norm errors on Linux are 1.09e-19 for both paths. Final
elapsed-time errors are zero at the recorded precision on Linux; macOS has
maximum 4.32e-14 Mpc from finite-precision absolute-time accumulation.

## Grid difference and convergence

The largest difference occurs in the zero-extension near-outside start,
positron/proton, accuracy 40, one outer call. The old global scale is
0.13333333333333336 Mpc; the modern scale is 0.13333333333333333 Mpc.
The old implementation computes extents from double internal-coordinate
bounds, while modern code retains wider source coordinates until narrowing
the final scale. This approximately one-bit rounding difference lies at an
integer step-count boundary: the fork takes 75 internal field steps and the
modern implementation 76 (76/77 field queries including the initial Rate call).

When the particle crosses from zero field into the grid, those partitions
sample the first nonzero field at different positions. The accumulated
position difference is about 0.183 kpc and the unit-direction difference
about 0.00111. A tiny scale difference therefore need not produce a tiny
trajectory difference at a discontinuous field boundary.

Maximum old/new differences across all 648 grid cases, grouped by numerical
controls, are:

| Accuracy | Outer calls | Position difference, Mpc | Direction-vector difference |
| ---: | ---: | ---: | ---: |
| 10 | 1 | 1.30e-17 | 1.08e-16 |
| 10 | 4 | 1.28e-17 | 1.08e-16 |
| 10 | 16 | 1.39e-17 | 1.13e-16 |
| 20 | 1 | 1.41e-17 | 1.15e-16 |
| 20 | 4 | 1.34e-17 | 1.09e-16 |
| 20 | 16 | 1.48e-17 | 1.17e-16 |
| 40 | 1 | 1.83e-4 | 1.11e-3 |
| 40 | 4 | 1.42e-17 | 1.13e-16 |
| 40 | 16 | 1.39e-17 | 1.12e-16 |

The nonmonotonic old/new difference is retained explicitly. Increasing
accuracy alone does not guarantee matching partitions; the outer-call
refinement is an independent numerical control. No cases or limits were
retuned after observing this difference. The fork's position-dependent
zero-scale shortcut was taken 225 times and is counted separately; this
bounded matrix did not demonstrate an additional endpoint error from it.

## Consequences and remaining scope

The rotation and direction/position update formulas agree. Remaining
deterministic differences are beta precision, scale evaluation/partition
rounding and the documented outside-zero-field shortcut. Modern true-minimum
scale selection also retains the previously documented correction for tied
axis scales; this unequal-spacing grid does not exercise that defect.

There is no blanket 100% trajectory-reproducibility claim. This baseline
motivated the explicitly approved production scale-arithmetic correction
recorded below, retaining the corrected true minimum and validation.
The results also motivate retaining independent outer-step refinement in later transport
or cascade comparisons. Source/observer distributions, interacting cascades,
statistical margins and whole-engine acceptance remain subsequent stages.

## Approved scale-rounding correction

After inspecting the baseline, the user approved matching the fork's scale
arithmetic. Production commit `5ee4527` computes each scale from the existing
double internal-coordinate box: convert both bounds, subtract, then divide
by node count. It retains a true minimum across all three axes, including
ties, and the finite-positive scale, regular-grid, coordinate and boundary
validation. It does not restore the historical zero-scale shortcut, alter
the integrator, or locate exact face crossings.

The new self-contained C++ regression
`GridFields.HistoricalScaleKeepsLegacyStepPartition` checks exact scale
identity and counts 75 field queries over a zero-field propagation interval.
The unchanged test was compiled against the previous implementation too:
it failed with 76 queries. It passes after the correction. The existing
true-minimum tie regression is retained unchanged.

All 756 cases were repeated with identical input hashes, reference hashes,
adapter hashes and numerical controls. The executions recorded clean source
revision `5ee4527`. [Post-fix evidence](fork_transport_scale_rounding.json)
preserves the original baseline hash, complete Linux convergence matrix,
raw-trajectory/executable hashes and the original worst case in each environment.
The generated post-fix files remain under ignored
`build/checks/scale-rounding-transport-*`; the baseline artifacts are retained.

| Environment | Completed cases | Maximum old/new position difference, Mpc | Maximum direction-vector difference |
| --- | ---: | ---: | ---: |
| macOS arm64 Release | 756 | 3.77490e-15 | 4.02456e-16 |
| Emulated Linux amd64 Release | 756 | 1.47720e-17 | 1.17769e-16 |
| Emulated Linux amd64 ASan+UBSan, leak detection | 756 | 1.47720e-17 | 1.17769e-16 |

All reported global scales now agree exactly with the fork in this matrix.
The original worst-case trajectory now takes 75 internal steps in both
implementations (76 field queries each including the initial Rate call).
The 0.183 kpc difference is removed; remaining endpoint differences are at
roundoff level for these controls. Sanitizer execution completed without
diagnostics. The input family, counts and settings were not tuned after
the baseline. No statistical or whole-cascade equivalence is inferred.

The complete rebuilt short CTest suites also passed after the production
change, including the new exact-partition regression and unchanged tie test:

| Environment | CTest passed | Elapsed |
| --- | ---: | ---: |
| macOS arm64 Release | 242/242 | 17.42 s |
| Emulated Linux amd64 Release | 242/242 | 123.83 s |
| Emulated Linux amd64 ASan+UBSan, leak detection | 242/242 | 173.14 s |

The six self-contained adapter tests passed again. The post-fix JSON evidence
also records hashes of these logs and the deliberately failing old-code
regression control. All build, temporary and diagnostic artifacts stayed
inside the project; external reference files were read-only.
