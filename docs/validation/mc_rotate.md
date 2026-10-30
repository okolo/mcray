# Stable mc_rotate delay verification

On 2026-10-30, tools/postprocess/mc_rotate.py was updated after the unchanged
historical copy committed as 37370fb. The subsequent block optimization keeps
those numerical formulas and adds the rotate_block array interface. External reference scripts were not
modified. The change preserves the sphere-extrapolation model and output
column order; output numbers now retain scalar precision. Materially invalid
sphere positions are diagnosed instead of being projected silently.

The current 18 deterministic unittest cases in tests/test_mc_rotate.py passed on:

| Platform | Python | NumPy | Result |
| --- | --- | --- | --- |
| macOS arm64 | 3.12.14 | 2.5.3 | 18/18 |
| Linux amd64, emulated | 3.12.3 | 1.26.4 | 18/18 |
| Linux arm64 | 3.11.16 | 1.26.4 | 18/18 |

The initial scalar version passed ten checks; the current block version
passes 18. Additional coverage checks block sizes 1, 2, 7, 64 and 32768,
byte-identical output for axial/tiny-delay and fixed-seed random histories,
comments and errors at block boundaries, consistent column counts, empty
blocks/selections, immutable input arrays, mixed root branches in one block,
and the exact 180-degree full-sky boundary on extended-precision platforms.
The latter avoids rounding pi*180/180 to a value just below pi.

The tests include:

- Axial trajectories with the analytic delay/T=d, down to 1e-100, including
  a write/read round trip. Positive delays below machine epsilon survive.
- Outward, inward and tangent directions and small non-axial deficits,
  compared against an independent 100-digit Decimal quadratic solution.
- 200 fixed-seed (724106) geometries compared against the same reference.
- Finite unit endpoints/directions; compensated summation; 1e-20 radian
  angles and both pole rotations; azimuth and PSF/jet selection.
- Arithmetic-roundoff clamping, explicit lower-text-precision input tolerance,
  rejection of material sphere violations and malformed rows, unchanged tail
  values, and rejection of using the input path as the output path.

Numerical tolerances scale with the platform's longdouble epsilon (typically
128 eps for relative root comparisons, 256 eps for the fixed-seed cases).
These are arithmetic regression checks, not statistical physical acceptance.
NumPy's reported storage size does not identify significand precision:
Linux amd64 uses the x87 extended format stored in 16 bytes; macOS arm64
longdouble has the same precision as double. No conversion through Python
float or math.fsum is used in the delay/formatting path.

The Linux amd64 run used a read-only container root, no network, a writable
project bind mount, and the existing image
`sha256:8d38c4e97825896aa9f2ae1d4bc66f2c5e5131925bb39fd3c544136fb8ba2f74`.
NumPy 1.26.4 was extracted inside the project's diagnostic directory from
numpy-1.26.4-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl,
SHA-256 `675d61ffbfa78604709862923189bad94014bef562cc35cf61d3a07bba02a7ed`.
The Linux arm64 image was
`sha256:5c787caa4ac270f92b4cc33b90188d69d348cdaaaababc35cc2fd4e9f612b6dd`.
Local logs and the downloaded wheel are retained outside Git under
build/mc-rotate-verification. All test temporary files stay inside the project.

The routine still assumes the original fixed-time, origin-source, +z-pencil
input and straight speed-c continuation. It does not account for future
interactions/fields, a massive particle's speed correction, or a cosmological
arrival-time transformation. Accuracy is limited by the input values and
cancellation in the radial deficit when its terms nearly cancel. Explicitly
widening the input tolerance does not improve that accuracy. The implemented
contract and invocation are in [the utility manual](../../tools/postprocess/mc_rotate.md).

## Block performance check

The file reader defaults to 32768 rows per block; rotate_block performs the
geometry with component arrays and compensated sums across each row. It does
not allocate a matrix for each particle. Runtime workspace is O(B*C) for block
size B and column count C. Python/NumPy allocator retention also contributes
to the process high-water resident size.

Synthetic ten-column inputs contained four repeated geometries, including
a tiny axial delay and non-axial directions. All rows were accepted in the
full-sky measurements. Sequential runs used macOS arm64 / NumPy 2.5.3, one
OpenBLAS thread, and local text input/output. Times exclude input generation:

| Rows | Scalar seconds | Block seconds | Scalar peak MiB | Block peak MiB |
| --- | --- | --- | --- | --- |
| 20,000 | 1.266 | 0.146 | 28.3 | 43.5 |
| 100,000 | 6.150 | 0.704 | 28.8 | 84.2 |
| 1,000,000 | Not run | 6.789 | Not measured | 96.9 |

For 100,000 rows this is approximately 8.7x faster, with increased bounded
workspace. With a narrow angular/PSF cut, 20,000 rows took 0.090 seconds
versus 1.302 seconds for the scalar version on the same input. These timings
are workload-specific local measurements, not native Linux performance claims.
No C++ or propagation physics changed. Source snapshots, benchmark driver,
JSON results and retained input/output are under build/mc-rotate-benchmark.
Final platform logs are under build/mc-rotate-verification/*-block-tests.log.
