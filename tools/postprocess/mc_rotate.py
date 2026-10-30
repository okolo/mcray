"""Rotate legacy +z-pencil histories and extend them to the sphere R=cT.

The source is at the origin, all primaries have the same birth redshift > 0,
and all histories end at redshift zero. T is the actual common flight time.
The extension assumes speed c and no further interactions or deflections.
Delays are dimensionless Delta_t/T. File processing uses bounded-size blocks.
"""

import io
import operator
from pathlib import Path

import numpy as np

_REAL = np.longdouble
_PI = np.arccos(_REAL(-1))
_EPS = np.finfo(_REAL).eps
_DEFAULT_BLOCK_SIZE = 32768


class _RowError(ValueError):
    """A numerical validation error with a zero-based row within the block."""

    def __init__(self, message, row):
        super().__init__(message)
        self.row = row


def _require(valid, message):
    """Find the first failing row; the file reader supplies its real line number."""
    invalid = np.flatnonzero(~np.atleast_1d(valid))
    if invalid.size:
        raise _RowError(message, int(invalid[0]))


def _sum(values):
    """Neumaier summation for scalars or arrays, independently for each row.

    Only the small number of terms is traversed in Python. NumPy processes
    all rows together, preserving longdouble throughout the compensation.
    """
    total = correction = _REAL(0)
    for value in values:
        value = np.asarray(value, dtype=_REAL)
        updated = total + value
        correction = correction + np.where(
            np.abs(total) >= np.abs(value),
            (total - updated) + value,
            (value - updated) + total)
        total = updated
    return total + correction


def _norm(vector):
    """Overflow-safe norm along the last axis (one vector or a block)."""
    vector = np.asarray(vector, dtype=_REAL)
    return np.hypot(np.hypot(vector[..., 0], vector[..., 1]), vector[..., 2])


def _unit(vector):
    vector = np.asarray(vector, dtype=_REAL)
    if vector.ndim not in (1, 2) or vector.shape[-1] != 3:
        raise ValueError("Expected three-component vectors")
    _require(np.isfinite(vector).all(axis=-1), "Expected finite vectors")
    scale = np.max(np.abs(vector), axis=-1, keepdims=True)
    _require(scale[..., 0] > 0, "Cannot normalize a zero vector")
    scaled = vector / scale
    return scaled / _norm(scaled)[..., np.newaxis]


def _angle(first, second):
    """Stable small angles without the overhead of general-purpose np.cross."""
    first, second = np.asarray(first, dtype=_REAL), np.asarray(second, dtype=_REAL)
    x = first[..., 1] * second[..., 2] - first[..., 2] * second[..., 1]
    y = first[..., 2] * second[..., 0] - first[..., 0] * second[..., 2]
    z = first[..., 0] * second[..., 1] - first[..., 1] * second[..., 0]
    cross_norm = np.hypot(np.hypot(x, y), z)
    dot = _sum(first[..., i] * second[..., i] for i in range(3))
    return np.arctan2(cross_norm, dot)


def _sphere_extension(x, y, d, direction, input_rtol=0):
    """Return endpoints, unit momenta, delay/T and corrected 1-z/T.

    Arguments can describe a single particle or an entire block. Positive
    radial deficits are retained even when smaller than machine epsilon.
    """
    x, y, d = (np.asarray(value, dtype=_REAL) for value in (x, y, d))
    input_rtol = _REAL(input_rtol)
    if not np.isfinite(input_rtol) or input_rtol < 0:
        raise ValueError("Nonnegative input_rtol must be finite")
    _require(np.isfinite(x) & np.isfinite(y) & np.isfinite(d),
             "Coordinates must be finite")
    tolerance = input_rtol + 8 * _EPS
    _require((d >= -tolerance) & (d <= 2 + tolerance),
             "1-z/T is outside the physical range [0, 2]")
    d = np.clip(d, _REAL(0), _REAL(2))
    _require((np.abs(x) <= 1 + tolerance) & (np.abs(y) <= 1 + tolerance),
             "Recorded position is outside the flight sphere")
    n = _unit(direction)

    # Compute 1-|r|^2 from the saved complement d, rather than subtracting
    # a reconstructed radius from one. This retains very small delays.
    a = d * (2 - d)
    delta = _sum((a, -x * x, -y * y))
    bound = tolerance * (2 + np.abs(a) + x * x + y * y)
    _require(delta >= -bound, "Recorded position is outside the flight sphere")
    delta = np.maximum(delta, _REAL(0))
    q = _sum((x * n[..., 0], y * n[..., 1], n[..., 2], -d * n[..., 2]))
    root = np.hypot(q, np.sqrt(delta))

    # Rationalize the outward root. A masked divide avoids evaluating 0/0
    # in tangent/inward rows, where the original root is already stable.
    delay = np.asarray(root - q).copy()
    np.divide(delta, root + q, out=delay, where=q > 0)
    endpoint = _unit(np.stack([
        _sum((x, delay * n[..., 0])),
        _sum((y, delay * n[..., 1])),
        _sum((_REAL(1), -d, delay * n[..., 2]))], axis=-1))
    return endpoint, n, delay, d


def _rotated_transverse(r, n):
    """Return x/y of momentum after rotating the endpoint r to +z.

    Rodrigues' rotation is evaluated by components. No (rows, 3, 3) matrix
    array is allocated. Only accepted particles reach this calculation.
    """
    r, n = np.asarray(r, dtype=_REAL), np.asarray(n, dtype=_REAL)
    transverse = np.hypot(r[..., 0], r[..., 1])
    off_axis = transverse > 0
    ux, uy = np.zeros_like(transverse), np.zeros_like(transverse)
    np.divide(r[..., 1], transverse, out=ux, where=off_axis)
    np.divide(-r[..., 0], transverse, out=uy, where=off_axis)
    cosine = r[..., 2]
    difference = np.asarray(1 - cosine).copy()
    # Retain 1-cos(theta) near the north pole; avoid division at the south pole.
    np.divide(transverse * transverse, 1 + cosine,
              out=difference, where=cosine >= 0)
    axis_dot = _sum((ux * n[..., 0], uy * n[..., 1]))
    final_x = _sum((cosine * n[..., 0], ux * axis_dot * difference,
                    uy * transverse * n[..., 2]))
    final_y = _sum((cosine * n[..., 1], uy * axis_dot * difference,
                    -ux * transverse * n[..., 2]))
    # At a pole, choose the original identity / pi-about-x convention.
    return np.where(off_axis, final_x, n[..., 0]), final_y


def _parameters(jet_half_size_deg, jet_direction_deg, psf_deg, input_rtol):
    half, direction, psf, tolerance = map(
        _REAL, (jet_half_size_deg, jet_direction_deg, psf_deg, input_rtol))
    if not np.isfinite([half, direction, psf, tolerance]).all() or tolerance < 0:
        raise ValueError("Angles and nonnegative input_rtol must be finite")
    if not 0 <= direction <= 180 or not 0 <= psf <= 180:
        raise ValueError("Jet direction and PSF must be in [0, 180] degrees")
    # Keep full-sky and antipodal directions exact: pi*180/180 can differ
    # from pi by one ulp on extended-precision platforms.
    half = half * _PI / 180 if 0 < half < 180 else _PI
    direction = _PI if direction == 180 else direction * _PI / 180
    jet_axis = np.array([0, np.sin(direction), np.cos(direction)], dtype=_REAL)
    return half, jet_axis, psf, tolerance


def rotate_block(particles, jet_half_size_deg, jet_direction_deg, psf_deg=180,
                 *, input_rtol=0):
    """Transform a numeric (rows, columns) block, retaining accepted row order.

    Input columns follow the legacy normalized layout; at least eight are
    required. The returned longdouble array has two fewer columns: energy,
    weight, Theta, Phi, raw delay/T, calculated delay/T, then the input tail.
    An empty selection returns a correctly shaped empty array. The input is
    never modified. Numerical errors carry a zero-based block row index.
    """
    half, jet_axis, psf, tolerance = _parameters(
        jet_half_size_deg, jet_direction_deg, psf_deg, input_rtol)
    particles = np.asarray(particles, dtype=_REAL)
    if particles.ndim != 2 or particles.shape[1] < 8:
        raise ValueError("Expected a two-dimensional block with at least eight columns")
    _require(np.isfinite(particles).all(axis=1), "Expected finite numeric columns")
    if not len(particles):
        return np.empty((0, particles.shape[1] - 2), dtype=_REAL)
    momentum = particles[:, 5:8].copy()
    momentum[:, 2] = 1 - momentum[:, 2]
    r, n, delay, d = _sphere_extension(
        particles[:, 2], particles[:, 3], particles[:, 4], momentum, tolerance)

    keep = np.ones(len(particles), dtype=bool)
    if half < _PI:
        # For the original rotation, R @ (0,0,1) = (-r_x,-r_y,r_z).
        # The endpoint is already unit-normalized; no rotation matrix is needed.
        initial = r * np.array([-1, -1, 1], dtype=_REAL)
        keep = _angle(initial, jet_axis) < half
    indices = np.flatnonzero(keep)
    beta = _angle(r[indices], n[indices]) * 180 / _PI
    psf_keep = beta <= psf
    indices, beta = indices[psf_keep], beta[psf_keep]
    final_x, final_y = _rotated_transverse(r[indices], n[indices])
    phi = np.where(beta > 0, np.arctan2(final_y, final_x) % (2 * _PI), _REAL(0))
    return np.column_stack((particles[indices, :2], beta, phi * 180 / _PI,
                            d[indices], delay[indices], particles[indices, 8:]))


def _parse_block(rows, line_numbers, filename, column_count):
    """Parse in NumPy; inspect individual lines only to diagnose a bad block."""
    try:
        data = np.loadtxt(io.StringIO("".join(rows)), dtype=_REAL,
                          ndmin=2, comments=None)
        if data.shape[1] == column_count and column_count >= 8:
            return data
    except ValueError:
        pass
    # This slower path runs only on invalid input, preserving precise diagnostics.
    for line, number in zip(rows, line_numbers):
        try:
            values = np.array(line.split(), dtype=_REAL)
            if values.size != column_count or column_count < 8:
                raise ValueError(f"Expected {column_count} columns (at least eight required)")
        except ValueError as error:
            raise ValueError(f"{filename}:{number}: {error}") from error
    raise ValueError(f"{filename}:{line_numbers[0]}: Invalid numeric block")


def _write_block(rows, line_numbers, stream, filename, column_count, parameters):
    data = _parse_block(rows, line_numbers, filename, column_count)
    half_deg, direction_deg, psf_deg, tolerance = parameters
    try:
        result = rotate_block(data, half_deg, direction_deg, psf_deg, input_rtol=tolerance)
    except _RowError as error:
        raise ValueError(f"{filename}:{line_numbers[error.row]}: {error}") from error
    # Format row by row to avoid building a second block-sized list of strings.
    # NumPy formats the actual scalar type without a float64 conversion.
    stream.writelines(" ".join(np.format_float_scientific(value, unique=True, trim="-")
                               for value in row) + "\n" for row in result)


def mc_rotate(mc_file, jet_half_size_deg, jet_direction_deg, psf_deg=180,
              rotated_file=None, *, input_rtol=0, block_size=_DEFAULT_BLOCK_SIZE):
    """Stream a legacy file through rotate_block, keeping only one block in RAM.

    block_size bounds the number of buffered data rows (default 32768).
    Numerical validation and output precision are independent of block size.
    input_rtol widens roundoff bounds for older lower-precision input; it
    cannot restore information lost when that input was written.
    """
    if rotated_file is None:
        rotated_file = f"{mc_file}_{jet_half_size_deg}_{jet_direction_deg}_{psf_deg}"
    if Path(mc_file).resolve() == Path(rotated_file).resolve():
        raise ValueError("Input and output paths must differ")
    if isinstance(block_size, bool):
        raise ValueError("block_size must be a positive integer")
    try:
        block_size = operator.index(block_size)
    except TypeError as error:
        raise ValueError("block_size must be a positive integer") from error
    if block_size <= 0:
        raise ValueError("block_size must be a positive integer")
    parameters = (jet_half_size_deg, jet_direction_deg, psf_deg, input_rtol)
    _parameters(*parameters)  # Validate before opening the output file.
    rows, line_numbers = [], []
    column_count = None

    with open(mc_file, "r") as lines, open(rotated_file, "w") as rotated:
        def flush():
            if rows:
                _write_block(rows, line_numbers, rotated, mc_file, column_count, parameters)
                rows.clear()
                line_numbers.clear()

        for number, line in enumerate(lines, 1):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                # Flush first so comments embedded in a file keep their position.
                flush()
                columns = stripped.split()
                if "weight" in columns:
                    names = columns[1:3] + ["Theta", "Phi", "dT_raw/T", "dT_calc/T"] + columns[9:]
                    rotated.write("# " + "\t".join(names) + "\n")
                else:
                    rotated.write(stripped + "\n")
                continue
            if column_count is None:
                column_count = len(stripped.split())
            rows.append(line)
            line_numbers.append(number)
            if len(rows) == block_size:
                flush()
        flush()
    return rotated_file
