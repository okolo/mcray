"""Loader and geometric delay regression checks; requires NumPy and PyYAML."""
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("file_helpers", ROOT / "tools/postprocess/file_helpers.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
CRbeamFile = module.CRbeamFile


class SpatialLoaderTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        (self.root / "z0").mkdir()
        self.file = self.root / "z0/photon"

    def write(self, scale="1e-05", rows=None, flight=True, metadata=""):
        columns = "E/eV\tweight\tx/T\ty/T\t1-(z/T)\tPx/P\tPy/P\t1-(Pz/P)\tz_production\tin_target\tsource_x\tsource_y\tsource_z\tPathMpc"
        if flight:
            columns += "\tFlightMpc"
        self.file.write_text("# photon\tT/Mpc = " + scale + metadata + "\n#\n# " + columns + "\n" + (rows or ""))

    def test_exponents_single_row_and_initial_direction(self):
        self.write(rows="1e12 2 3 4 -4 0 0 0 0 1 0 0 1 7 8\n")
        data = CRbeamFile(self.file)
        self.assertEqual(len(data), 1)
        self.assertAlmostEqual(data.Norm, .01)
        np.testing.assert_allclose([data.x[0], data.y[0], data.z[0]], [.03, .04, .05])
        self.assertEqual(data["source_z"][0], 1)
        self.assertEqual(data.in_target_num, 1)
        self.write(scale="1.234e+06")
        self.assertEqual(CRbeamFile(self.file).scale_Mpc, 1234000)

    def test_yaml_and_legacy_cli_units(self):
        self.write(rows="1 1 0 0 1 0 0 0 0 0 0 0 1 0 0\n")
        (self.root / "params").write_text("--source-x\t1200\n--source-y\t-3\n--observer-radius\t5\n")
        legacy = CRbeamFile(self.file)
        np.testing.assert_allclose(legacy.source_position_Mpc, [1.2, -.003, 0])
        (self.root / "config.yaml").write_text("schema_version: 2\nsource:\n  redshift: 0.01\n  position_Mpc: [2, -1, 3]\n  angular_distribution: {half_angle_deg: 30}\nobserver:\n  center_Mpc: [0, 0, 4]\n  radius_Mpc: 0.5\n")
        data = CRbeamFile(self.file)
        np.testing.assert_allclose(data.source_position_Mpc, [2, -1, 3])
        self.assertEqual(data["observer-z"], 4000)
        self.assertEqual(data["observer-radius"], 500)
        self.assertAlmostEqual(data["jet-angle"], np.pi / 6)

    def test_delay_uses_each_flight_and_displacement_from_source(self):
        self.write(scale="10", metadata="\tnormalization = fallback\tspatial_format = 1\tsource_x_Mpc = 2\tsource_y_Mpc = 0\tsource_z_Mpc = 0",
                   rows="1 1 .5 0 .6 0 0 0 0 1 0 0 1 1 6\n1 1 .5 0 .6 0 0 0 0 1 0 0 1 1 7\n")
        data = CRbeamFile(self.file)
        np.testing.assert_allclose(data.delay_Mpc(), [1, 2])
        np.testing.assert_allclose(data.delay_Mpc("axis", [0, 0, 5]), [2, 3])
        np.testing.assert_allclose(data["source_z"], [1, 1])
        with self.assertRaises(ValueError):
            data.delay_Mpc("axis", [0, 0, 0])
        result = subprocess.run([sys.executable, str(ROOT / "tools/postprocess/flight_delay.py"), str(self.file)],
                                check=True, capture_output=True, text=True)
        self.assertIn("delay_Mpc", result.stdout)

    def test_legacy_scale_is_never_assumed_to_be_particle_flight(self):
        self.write(flight=False, rows="1 1 0 0 1 0 0 0 0 0 0 0 1 0\n")
        with self.assertRaisesRegex(ValueError, "FlightMpc"):
            CRbeamFile(self.file).delay_Mpc()

    @unittest.skipUnless(os.environ.get("MCRAY_CRBEAM"), "Set MCRAY_CRBEAM for the executable round trip")
    def test_production_birth_and_exit_files_round_trip(self):
        binary = os.environ["MCRAY_CRBEAM"]
        for stop in ["inside", "outside"]:
            destination = self.root / stop
            command = [binary, "--particle", "proton", "--z", ".001", "--emax", "1e12",
                       "--nparticles", "2", "--seed", "42", "--batch", "2", "--threads", "1",
                       "--background", "0", "--EGMF", "0", "--noEMcascade",
                       "--source-position-Mpc", "3,0,-2", "--observer", "sphere",
                       "--observer-center-Mpc", "3,0,-2", "--observer-radius-Mpc", ".001",
                       "--observer-stop-when", stop, "--output-primary-direction",
                       "--output", str(destination)]
            subprocess.run(command, cwd=ROOT / "bin", check=True, capture_output=True, text=True)
            data = CRbeamFile(destination / "z0/proton")
            self.assertEqual(len(data), 2)
            self.assertEqual(data.in_target_num, 2)
            if stop == "inside":
                self.assertEqual(data.scale_Mpc, 3)
                np.testing.assert_allclose(data.delay_Mpc(), 0, atol=1e-15)
                np.testing.assert_allclose(data["FlightMpc"], 0)
            else:
                self.assertTrue(np.all(data["FlightMpc"] > 0))
                self.assertTrue(np.all(data.delay_Mpc() > 0))
            np.testing.assert_allclose(data["source_z"], 1)

    def test_empty_and_malformed_tables(self):
        self.write()
        self.assertEqual(len(CRbeamFile(self.file)), 0)
        for scale, row in [("0", ""), ("nan", ""), ("1", "1 2\n"),
                           ("1", "1 1 nan 0 1 0 0 0 0 0 0 0 1 0 0\n")]:
            with self.subTest(scale=scale, row=row):
                self.write(scale, row)
                with self.assertRaises(ValueError):
                    CRbeamFile(self.file)


if __name__ == "__main__":
    unittest.main()
