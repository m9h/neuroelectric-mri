"""Reproducible tank-release contract, distinct from numerical solver tests."""

import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

AVAILABLE = all(importlib.util.find_spec(n) is not None for n in ("numpy", "scipy"))
if AVAILABLE:
    import numpy as np
    from examples.tank_benchmark import write_benchmark


@unittest.skipUnless(AVAILABLE, "Install requirements-tank.txt")
class TankBenchmarkTests(unittest.TestCase):
    def test_release_has_si_fields_conservation_and_provenance(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder)
            write_benchmark(output, convergence_sizes=(4, 8))
            manifest = json.loads((output / "manifest.json").read_text())
            self.assertEqual(manifest["schema"], "neuroelectric-mri.synthetic-tank/0.1")
            self.assertEqual(manifest["data_kind"], "synthetic_numerical")
            self.assertFalse(manifest["scanner_data_available"])
            self.assertFalse(manifest["complete_circuit_field_available"])
            for filename, digest in manifest["files_sha256"].items():
                self.assertEqual(hashlib.sha256((output / filename).read_bytes()).hexdigest(), digest)
            self.assertIn("solver", manifest["source_sha256"])
            parameters = json.loads((output / "parameters.json").read_text())
            self.assertEqual(parameters["coordinate_frame"], "right-handed Cartesian, axes x/y/z")
            with np.load(output / "fields.npz", allow_pickle=False) as fields:
                sigma = fields["inclusion_conductivity_s_m"]
                j = fields["inclusion_current_density_a_m2"]
                e = fields["inclusion_electric_field_v_m"]
                self.assertEqual(j.shape, (*sigma.shape, 3))
                np.testing.assert_allclose(j, sigma[..., None] * e)
                balance = sum(np.diff(fields[f"inclusion_face_current_{axis}_a"], axis=i)
                              for i, axis in enumerate("xyz"))
                np.testing.assert_allclose(balance, 0, atol=1e-14)
                np.testing.assert_allclose(fields["scaled_inclusion_current_density_a_m2"], j,
                                            atol=1e-11)
                np.testing.assert_array_equal(fields["current_off_current_density_a_m2"], np.zeros_like(j))
            with (output / "convergence.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual([int(r["cells_per_axis"]) for r in rows], [4, 8])
            self.assertLess(float(rows[1]["potential_relative_l2_error"]),
                            float(rows[0]["potential_relative_l2_error"]))

    def test_generated_arrays_and_metrics_repeat(self):
        with tempfile.TemporaryDirectory() as folder:
            a, b = Path(folder) / "a", Path(folder) / "b"
            for output in (a, b):
                write_benchmark(output, convergence_sizes=(4, 8))
            for name in ("metrics.csv", "convergence.csv", "parameters.json", "manifest.json"):
                self.assertEqual((a / name).read_bytes(), (b / name).read_bytes())
            with np.load(a / "fields.npz") as left, np.load(b / "fields.npz") as right:
                self.assertEqual(set(left.files), set(right.files))
                for key in left.files:
                    np.testing.assert_array_equal(left[key], right[key])

    def test_invalid_convergence_grids_are_rejected(self):
        for sizes in ((8,), (8, 4), (4, 4), (0, 8), (4, 8.5)):
            with tempfile.TemporaryDirectory() as folder, self.subTest(sizes=sizes):
                with self.assertRaises(ValueError):
                    write_benchmark(Path(folder), convergence_sizes=sizes)


if __name__ == "__main__":
    unittest.main()
