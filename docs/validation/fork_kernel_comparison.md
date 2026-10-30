# Direct historical-fork field and cone comparison

On 2026-10-30, the hash-pinned historical `CRbeam_v0` field and uniform-cone
functions were compiled and run against the modern implementations. The goal
is agreement of observables within Monte Carlo uncertainty, not equal seeded
trajectories. This record concerns these kernels, not complete cascade spectra.

The runnable [adapter and frozen design](../../tests/fork_compat/README.md)
explain extraction, reference hashes, input grids, seeds, sample counts,
observables and simultaneous diagnostic intervals. The reference functions
are extracted verbatim; they are not independently rewritten formulas.

## Results and implemented correction

Two synthetic grids, all four boundary modes and 16192 matched query points
exposed systematic field differences before correction. For fields of order
one Gauss, the largest component difference was 1.3333333333333333 Gauss;
strict interiors already agreed within 7.11e-15 Gauss.

The modern implementation now retains the fork's:

- Strict interior: a face invokes the outside handler.
- Nearest **node** on all axes when any coordinate is outside/on a face.
- Periodic/mirrored mapping of all axes and the `width/N` seam inset.
- Double-precision box/period arithmetic and subtraction-before-remainder
  order. Splitting the remainder can change which side of a seam is selected.
- Extent/node-count convention for the positive deflector scale.

Stable weighted-corner interpolation is retained: it is the same trilinear
polynomial, with better conditioning than expanded coefficients. The
[usage contract](../usage/grid_magnetic_fields.md) describes these rules and
their discontinuities. They replace the initially implemented closed-box,
clamp-and-interpolate behavior; existing runs using those extensions can change.

After correction, the largest component difference over all 16192 native
queries is 7.105427357601002e-15 Gauss, including boundary/outside controls.
On Linux amd64 the maximum is 1.709743457922741e-14 Gauss.
The two tested step scales also agree (1/3 and 1/2 Mpc). Bitwise equality is
neither required nor generally obtained.

For the uniform cone, the six ordinary cases (20, 90 and 180 degrees, each
on two axes) used 100000 primary draws **per implementation per case**.
All weights were one. None of their 156 tested differences excluded zero
under the joint 95% intervals allocated over all 260 comparisons. The
maximum standardized difference is 2.590468 on both platforms; the Bonferroni
diagnostic cutoff is 3.728913. For the normalized polar mean, the simultaneous
interval half-width is approximately 0.004814. This states the precision
reached; it does not exclude smaller systematic errors or constitute a
physical equivalence margin.

At 1e-6 radians the fork loses 13 of 100000 weights on each tested axis
(mean weight 0.99987), while the modern code preserves every weight. Their
joint intervals still contain zero, with maximum standardized difference
3.60634, but source inspection and the measured zero weights identify a
systematic rounding defect. A non-significant diagnostic must not conceal
that known algorithmic difference. At 1e-10 radians the loss becomes total.
These four narrow-cone cases are bug controls, not a claim of exact source
equivalence.

The maintained RNG and orthonormal source frame are unchanged. Their random
sequences and azimuth origins differ from the fork, but the ordinary cone
law is uniform in solid angle in both implementations.

## Deliberate exceptions

- At half-angle 1e-10 radians, the compiled fork gives every primary weight
  zero; the modern code preserves weight one. In the fork, double `cos(angle)`
  rounds to one and `coneZ_weight` rejects equality. This is an actual flux
  loss, not Monte Carlo noise. Stable `2*sin(angle/2)^2` sampling and removal
  of this spurious weight gate remain intentional numeric corrections.
- The fork's 1e-14-radian activation cutoff and failure to orient a zero-width
  source are not restored. Modern explicit pencil/axis and zero-opening
  contracts remain as approved.
- Actual grid extrema and strict input validation remain. The fork initializes
  extrema to zero and can read uninitialized nodes for positive-only or
  negative-only boxes. Undefined behavior is not a physical reference.
- Scale ties take the true minimum; the fork's strict comparisons can select
  a larger z scale. A separate compiled 2x2x2 constant-field control with
  extents (2,2,4) Mpc returns 2 Mpc in the fork and 1 Mpc in the modern code.
  No zero sentinel is returned outside zero-filled boxes:
  the modern deflector needs a positive subdivision scale to handle entry.
- Overflow guards and stable interpolation remain. Extremely large finite
  coordinates use a split remainder only if subtraction overflows.

The tests do not establish agreement of Gaussian jets, the fork's differently
weighted vMF branch, observers, full interaction cascades or final energy/time
spectra. No claim of whole-program “100% reproducibility” follows from these
kernel checks. In particular, discontinuous fields still need propagation
step-convergence checks for any intended physical calculation.

## Regression checks and provenance

| Final check | Result |
| --- | --- |
| Native macOS arm64 Release, complete CTest | 234/234 passed; 18.13 s. |
| Ubuntu 22.04 Linux amd64 Release, complete CTest | 234/234 passed; 113.84 s. |
| Linux amd64 ASan+UBSan, complete CTest | 234/234 passed; 164.03 s; no sanitizer/leak report. |

Native compiler: Homebrew Clang 22.1.5, C++17/libc++. Linux compiler:
GCC 11.4.0. The offline Ubuntu image is
`sha256:37682238eeaf9f01f7eddea17d4e55551d2e317307dfa5c8ea8dd0dc8084baf6`.
Linux runs were emulated on the macOS arm64 workstation; these are not
native Linux performance measurements. Each Docker check was bounded below
30 minutes, with build/temp/cache files inside the project.

The final sanitizer build instruments the modern C/C++/Fortran project with
`-O1 -g -fsanitize=address,undefined -fno-sanitize-recover=all
-fno-omit-frame-pointer -fno-pie` and links non-PIE executables. It uses
`ASAN_OPTIONS=detect_leaks=1:halt_on_error=1` and
`UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1`, without suppressions.
The historical reference probe was compared in Release; it is not an
ASan/UBSan audit of the entire fork or external libraries.

Compact measured results and source/input hashes are retained in
[fork_kernel_comparison.json](fork_kernel_comparison.json). Full query rows,
compiler commands, extracted sources and logs remain local ignored artifacts
under `build/fork-comparison`; they are not public build dependencies. The
before-correction source is the grid kernel at `c5f86a4`; measured working-tree
hashes identify the corrected source independently of its eventual commit ID.

The optional reference adapter is not required by ordinary CI. Thirteen
self-contained `unit.magnetic.*` cases cover the current behavior, including
all zero faces, nearest-node extension, seam insets, translated boxes, scale
ties, extreme-coordinate clamping, robust loading, concurrency and constant-field electron/Fe56 transport.
