"""Numerical and file-contract checks for legacy sphere-extrapolated delays."""
from decimal import Decimal, localcontext
import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("mc_rotate", ROOT / "tools/postprocess/mc_rotate.py")
rotate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rotate)
REAL = np.longdouble
EPS = np.finfo(REAL).eps


def decimal_value(value):
    return Decimal(np.format_float_scientific(REAL(value), unique=True, trim="-"))


def reference_delay(x, y, d, direction):
    """Independent 100-digit quadratic solution using normalized input direction."""
    with localcontext() as context:
        context.prec = 100
        x, y, d = map(decimal_value, (x, y, d))
        n = list(map(decimal_value, direction))
        length = sum(value * value for value in n).sqrt()
        n = [value / length for value in n]
        position = [x, y, 1 - d]
        # Deliberately form the original equation at high precision.
        b = sum(value * component for value, component in zip(position, n))
        c = sum(value * value for value in position) - 1
        return (b * b - c).sqrt() - b


class McRotateNumericsTest(unittest.TestCase):
    def assertRelative(self, actual, expected, factor=128):
        self.assertLessEqual(abs(decimal_value(actual) - expected),
                             abs(expected) * decimal_value(factor * EPS))

    def test_axial_delays_remain_positive_far_below_epsilon(self):
        for text in ["0.2", "1e-12", "1e-18", "1e-30", "1e-100"]:
            with self.subTest(delay=text):
                d = REAL(text)
                r, n, delay, raw = rotate._sphere_extension(0, 0, d, [0, 0, 1])
                self.assertGreater(delay, 0)
                self.assertRelative(delay, decimal_value(d))
                self.assertEqual(raw, d)
                np.testing.assert_allclose(r, [0, 0, 1], rtol=0, atol=8 * EPS)

    def test_both_root_branches_and_tangent_match_decimal_reference(self):
        cases = [
            (".1", ".2", ".1", [".3", "-.1", "1"]),
            (".1", "0", ".2", [".2", "0", "-1"]),
            ("0", "0", ".2", ["1", "0", "0"]),
            ("1e-10", "0", "1e-20", [".2", "-.3", "1"]),
            ("1e-20", "2e-20", "1e-30", ["0", "1", "1e-10"]),
        ]
        for x, y, d, n in cases:
            with self.subTest(x=x, d=d, n=n):
                r, unit, delay, _ = rotate._sphere_extension(x, y, d, n)
                self.assertRelative(delay, reference_delay(x, y, d, n))
                self.assertLessEqual(abs(rotate._norm(r) - 1), 8 * EPS)
                self.assertLessEqual(abs(rotate._norm(unit) - 1), 8 * EPS)

    def test_fixed_seed_cases_match_independent_high_precision_roots(self):
        rng = np.random.default_rng(724106)
        for _ in range(200):
            x, y = rng.uniform(-.2, .2, 2).astype(REAL)
            d = REAL(rng.uniform(.1, .8))
            direction = rng.normal(size=3).astype(REAL)
            _, _, delay, _ = rotate._sphere_extension(x, y, d, direction)
            self.assertRelative(delay, reference_delay(x, y, d, direction), factor=256)

    def test_roundoff_bounds_reject_real_violations_and_preserve_positive_deficit(self):
        for x, y, d in [(0, 0, "-.001"), (0, 0, "2.001"), (".1", 0, 0)]:
            with self.subTest(x=x, d=d), self.assertRaises(ValueError):
                rotate._sphere_extension(x, y, d, [0, 0, 1])
        _, _, delay, d = rotate._sphere_extension(0, 0, -4 * EPS, [0, 0, 1])
        self.assertEqual(delay, 0)
        self.assertEqual(d, 0)
        # Rounded transverse coordinates can put a true boundary point just outside.
        rotate._sphere_extension(np.sqrt(REAL(3)) / 2 * (1 + 2 * EPS), 0, ".5", [0, 0, 1])
        with self.assertRaises(ValueError):
            rotate._sphere_extension(".1", 0, ".004999", [0, 0, 1])
        rotate._sphere_extension(".1", 0, ".004999", [0, 0, 1], input_rtol="1e-4")
        for direction in [[0, 0, 0], [0, np.nan, 1]]:
            with self.assertRaises(ValueError):
                rotate._sphere_extension(0, 0, ".2", direction)

    def test_compensated_sum_keeps_scalar_precision(self):
        tiny = REAL("1e-30")
        self.assertEqual(rotate._sum([REAL(1), tiny, REAL(-1)]), tiny)

    def test_tiny_angle_and_pole_rotations(self):
        r, n, _, _ = rotate._sphere_extension("1e-20", 0, "1e-30", [0, 0, 1])
        angle = rotate._angle(r, n)
        self.assertRelative(angle, Decimal("1e-20"))
        final_x, final_y = rotate._rotated_transverse(r, n)
        self.assertLess(final_x, 0)  # The tiny rotation survives cos(theta)==1.
        self.assertLessEqual(abs(final_x / REAL("-1e-20") - 1), 128 * EPS)
        self.assertEqual(final_y, 0)
        for point in [[0, 0, 1], [0, 0, -1]]:
            final_x, final_y = rotate._rotated_transverse(point, [0, 0, 1])
            self.assertEqual(final_x, 0)
            self.assertEqual(final_y, 0)


def reference_rotation(r):
    """Explicit scalar matrix, independent of the block component formulas."""
    transverse = np.hypot(r[0], r[1])
    if transverse == 0:
        return np.diag(np.array([1, 1 if r[2] >= 0 else -1,
                                 1 if r[2] >= 0 else -1], dtype=REAL))
    ux, uy = r[1] / transverse, -r[0] / transverse
    cosine, sine = r[2], transverse
    difference = sine * sine / (1 + cosine) if cosine >= 0 else 1 - cosine
    return np.array([
        [cosine + ux * ux * difference, ux * uy * difference, uy * sine],
        [ux * uy * difference, cosine + uy * uy * difference, -ux * sine],
        [-uy * sine, ux * sine, cosine]], dtype=REAL)


class McRotateBlockTest(unittest.TestCase):
    def particles(self):
        rng = np.random.default_rng(724106)
        rows = np.zeros((200, 10), dtype=REAL)
        rows[:, 0], rows[:, 1] = 1e12, 1
        rows[:, 2:4] = rng.uniform(-.2, .2, (200, 2))
        rows[:, 4] = rng.uniform(.1, .8, 200)
        rows[:, 5:8] = rng.normal(size=(200, 3))
        rows[:, 7] = 1 - rows[:, 7]
        rows[:, 9] = np.arange(200)
        return rows

    def test_vector_roots_match_decimal_and_input_is_unchanged(self):
        rows = self.particles()
        original = rows.copy()
        result = rotate.rotate_block(rows, 180, 0)
        self.assertEqual(result.shape, (200, 8))
        np.testing.assert_array_equal(rows, original)
        np.testing.assert_array_equal(result[:, -1], rows[:, -1])
        for row, output in zip(rows, result):
            direction = row[5:8].copy()
            direction[2] = 1 - direction[2]
            expected = reference_delay(*row[2:5], direction)
            self.assertLessEqual(abs(decimal_value(output[5]) - expected),
                                 abs(expected) * decimal_value(256 * EPS))

    def test_component_rotation_and_selection_match_scalar_matrices(self):
        rows = self.particles()
        for half, tilt, psf in [(180, 0, 180), (50, 20, 60), (1, 0, .01)]:
            expected = []
            axis = np.array([0, np.sin(REAL(tilt) * rotate._PI / 180),
                             np.cos(REAL(tilt) * rotate._PI / 180)])
            for row in rows:
                direction = row[5:8].copy()
                direction[2] = 1 - direction[2]
                r, n, delay, d = rotate._sphere_extension(*row[2:5], direction)
                matrix = reference_rotation(r)
                initial = matrix @ np.array([0, 0, 1], dtype=REAL)
                initial /= rotate._norm(initial)
                start = np.arctan2(rotate._norm(np.cross(initial, axis)), np.dot(initial, axis))
                if half < 180 and start >= REAL(half) * rotate._PI / 180:
                    continue
                beta = np.arctan2(rotate._norm(np.cross(r, n)), np.dot(r, n))
                if beta * 180 / rotate._PI > psf:
                    continue
                final = matrix @ n
                phi = np.arctan2(final[1], final[0]) % (2 * rotate._PI) if beta > 0 else REAL(0)
                expected.append([*row[:2], beta * 180 / rotate._PI,
                                 phi * 180 / rotate._PI, d, delay, *row[8:]])
            result = rotate.rotate_block(rows, half, tilt, psf)
            expected = np.asarray(expected, dtype=REAL).reshape((-1, 8))
            np.testing.assert_allclose(result, expected, rtol=1024 * EPS, atol=1024 * EPS)

    def test_mixed_small_delays_and_root_branches_in_one_block(self):
        rows = np.array([
            [1, 1, 0, 0, REAL("1e-100"), 0, 0, 0],
            [1, 1, 0, 0, REAL(".2"), 1, 0, 1],
            [1, 1, 0, 0, REAL(".2"), 0, 0, 2],
            [1, 1, 0, 0, 0, 0, 0, 0],
        ], dtype=REAL)
        with np.errstate(divide="raise", invalid="raise"):
            result = rotate.rotate_block(rows, 180, 0)
        np.testing.assert_allclose(result[:, 5],
                                   np.array(["1e-100", ".6", "1.8", "0"], dtype=REAL),
                                   rtol=128 * EPS, atol=0)

    def test_empty_blocks_and_rejected_blocks_keep_column_shape(self):
        self.assertEqual(rotate.rotate_block(np.empty((0, 10)), 180, 0).shape, (0, 8))
        rows = self.particles()
        self.assertEqual(rotate.rotate_block(rows, 1, 0, .01).shape, (0, 8))
        with self.assertRaises(ValueError):
            rotate.rotate_block(np.empty((2, 7)), 180, 0)


class McRotateFileTest(unittest.TestCase):
    def setUp(self):
        # Keep every test artifact inside the development project.
        directory = ROOT / "build/mc-rotate-tests"
        directory.mkdir(parents=True, exist_ok=True)
        self.temporary = tempfile.TemporaryDirectory(dir=directory)
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.source = self.root / "input.tsv"
        self.target = self.root / "rotated.tsv"

    def write(self, rows):
        self.source.write_text("# photon\tT/Mpc = 1\n#\n"
                               "# E/eV\tweight\tx/T\ty/T\t1-(z/T)\tPx/P\tPy/P\t1-(Pz/P)\tz_production\tId\n" + rows)

    def run_rotation(self, **kwargs):
        rotate.mc_rotate(self.source, 180, 0, rotated_file=self.target, **kwargs)
        lines = self.target.read_text().splitlines()
        rows = [line.split() for line in lines if line and not line.startswith("#")]
        return lines, np.array(rows, dtype=REAL)

    def test_output_round_trips_small_delays_and_tail_columns(self):
        delays = [".2", "1e-12", "1e-18", "1e-30", "1e-100"]
        self.write("".join(f"1.234567890123456789e12 2 0 0 {d} 0 0 0 .001 123\n" for d in delays))
        lines, data = self.run_rotation()
        self.assertIn("Theta\tPhi\tdT_raw/T\tdT_calc/T\tz_production\tId", lines[2])
        self.assertEqual(data.shape, (len(delays), 8))
        for row, d in zip(data, delays):
            self.assertGreater(row[5], 0)
            self.assertLessEqual(abs(row[5] / REAL(d) - 1), 128 * EPS)
            self.assertEqual(row[0], REAL("1.234567890123456789e12"))
            self.assertEqual(row[6], REAL(".001"))
            self.assertEqual(row[7], 123)

    def test_tiny_angle_phi_and_psf_cut(self):
        self.write("1 1 1e-20 0 1e-30 0 0 0 0 1\n")
        _, data = self.run_rotation()
        self.assertGreater(data[0, 2], 0)
        self.assertLessEqual(abs(data[0, 2] / (REAL("1e-20") * 180 / rotate._PI) - 1), 128 * EPS)
        self.assertLessEqual(abs(data[0, 3] - 180), 128 * EPS * 180)
        _, cut = self.run_rotation(psf_deg=REAL("1e-20"))
        self.assertEqual(cut.size, 0)

    def test_non_axial_reference_and_jet_selection(self):
        self.write("1 1 .1 0 .2 0 0 0 0 1\n")
        _, data = self.run_rotation()
        with localcontext() as context:
            context.prec = 100
            expected = Decimal(".99").sqrt() - Decimal(".8")
        self.assertLessEqual(abs(decimal_value(data[0, 5]) - expected), 128 * decimal_value(EPS))
        rotate.mc_rotate(self.source, 1, 0, rotated_file=self.target)
        self.assertTrue(all(line.startswith("#") for line in self.target.read_text().splitlines()))

    def test_malformed_input_and_same_path_are_diagnosed(self):
        for row in ["1 1 .1 0 0 0 0 0 0 1\n", "1 1 nan 0 .2 0 0 0\n", "1 2\n"]:
            self.write(row)
            with self.assertRaisesRegex(ValueError, "input.tsv:4"):
                self.run_rotation()
        self.write("1 1 0 0 0 0 0 0 0 1\n")
        original = self.source.read_bytes()
        with self.assertRaises(ValueError):
            rotate.mc_rotate(self.source, 180, 0, rotated_file=self.source)
        self.assertEqual(self.source.read_bytes(), original)

    def test_block_size_does_not_change_output_or_row_order(self):
        self.write("".join(f"1 1 0 0 {d} 0 0 0 0 {i}\n"
                           for i, d in enumerate([".2", "1e-12", "1e-18", "1e-30", "1e-100"] * 5)))
        expected = None
        for size in [1, 2, 7, 32768]:
            self.run_rotation(block_size=size)
            text = self.target.read_bytes()
            if expected is None:
                expected = text
            self.assertEqual(text, expected)

    def test_random_rows_have_identical_output_for_different_blocks(self):
        rows = McRotateBlockTest().particles()
        self.write("".join(" ".join(np.format_float_scientific(value, unique=True, trim="-")
                                   for value in row) + "\n" for row in rows))
        expected = None
        for size in [1, 7, 64, 32768]:
            self.run_rotation(block_size=size)
            text = self.target.read_bytes()
            if expected is None:
                expected = text
            self.assertEqual(text, expected)

    def test_comments_empty_files_and_errors_across_blocks(self):
        self.write("1 1 0 0 .2 0 0 0 0 1\n\n# middle\n1 1 0 0 1e-18 0 0 0 0 2\n")
        outputs = []
        for size in [1, 4]:
            self.run_rotation(block_size=size)
            outputs.append(self.target.read_bytes())
        self.assertEqual(outputs[0], outputs[1])
        self.assertIn(b"# middle\n", outputs[0])
        for bad in ["1 1 .1 0 0 0 0 0 0 2\n", "1 1 0 0 nan 0 0 0 0 2\n",
                    "1 2\n", "invalid numbers\n"]:
            for size in [1, 2, 32768]:
                self.write("1 1 0 0 .2 0 0 0 0 1\n\n" + bad)
                with self.assertRaisesRegex(ValueError, "input.tsv:6"):
                    self.run_rotation(block_size=size)
        self.write("")
        self.run_rotation(block_size=2)
        self.assertTrue(all(line.startswith("#") for line in self.target.read_text().splitlines()))

    def test_invalid_block_size_is_rejected_before_output_creation(self):
        self.write("1 1 0 0 .2 0 0 0 0 1\n")
        for size in [0, -1, True, 1.5]:
            with self.assertRaises(ValueError):
                self.run_rotation(block_size=size)
            self.assertFalse(self.target.exists())


if __name__ == "__main__":
    unittest.main()
