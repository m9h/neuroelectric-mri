import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

AVAILABLE = all(importlib.util.find_spec(p) is not None for p in ("dolfinx", "gmsh", "scipy"))
if AVAILABLE:
    import numpy as np
    from examples.fem_benchmark import write_benchmark


@unittest.skipUnless(AVAILABLE, "Optional DOLFINx/Gmsh environment required")
class FEMReleaseTests(unittest.TestCase):
    def test_release_oracles_refinement_fields_and_hashes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_benchmark(root)
            manifest = json.loads((root / "manifest.json").read_text())
            self.assertTrue(manifest["synthetic"])
            self.assertEqual(manifest["solver"], "DOLFINx CEM / SciPy serial")
            for path, digest in manifest["files_sha256"].items():
                self.assertEqual(hashlib.sha256((root / path).read_bytes()).hexdigest(), digest)
            rows = list(csv.DictReader((root / "metrics.csv").open()))
            self.assertEqual(len(rows), 10)
            patches = [r for r in rows if r["case"] == "finite_patch"]
            for r in rows:
                self.assertLess(float(r["relative_residual"]), 1e-10)
                self.assertLess(float(r["relative_power_error"]), 1e-9)
                self.assertLess(float(r["max_electrode_current_error_a"]), 1e-11)
                if r["case"] != "finite_patch":
                    self.assertLess(float(r["voltage_relative_error"]), 1e-9)
            # Nested P1/P2 spaces on each mesh give a variational resistance bound.
            for p1, p2 in zip(patches[::2], patches[1::2]):
                self.assertGreater(float(p2["voltage_v"]), float(p1["voltage_v"]))
            p2 = [float(r["voltage_v"]) for r in patches if r["degree"] == "2"]
            self.assertLess(abs(p2[-1] - p2[-2]) / p2[-1], .01)
            self.assertLess(abs(p2[-1] - p2[-2]), abs(p2[1] - p2[0]))
            with np.load(root / "fields.npz", allow_pickle=False) as data:
                self.assertEqual(data["potential_cell_dofs"].shape[1], 10)
                self.assertEqual(data["vector_cell_dofs"].shape[1], 4)
                self.assertEqual(data["cell_tags"].size, data["potential_cell_dofs"].shape[0])
                np.testing.assert_allclose(data["potential_reversed_v"], -data["potential_v"], atol=1e-12)
                np.testing.assert_array_equal(data["potential_off_v"], 0.)
                self.assertTrue(np.all(np.isfinite(data["slice_j_a_m2"])))
                self.assertEqual(data["slice_j_a_m2"].shape[-1], 3)

    def test_repeated_release_is_numerically_reproducible(self):
        with tempfile.TemporaryDirectory() as directory:
            first, second = Path(directory) / "a", Path(directory) / "b"
            for output in (first, second):
                write_benchmark(output, mesh_sizes_m=(.03, .02))
            self.assertEqual((first / "metrics.csv").read_bytes(), (second / "metrics.csv").read_bytes())
            with np.load(first / "fields.npz") as a, np.load(second / "fields.npz") as b:
                self.assertEqual(set(a.files), set(b.files))
                for key in a.files:
                    np.testing.assert_array_equal(a[key], b[key])

    def test_invalid_refinement_sizes(self):
        for sizes in ((), (.02,), (.02, .03), (.03, 0), (.03, np.nan)):
            with tempfile.TemporaryDirectory() as directory:
                with self.subTest(sizes=sizes), self.assertRaises(ValueError):
                    write_benchmark(directory, mesh_sizes_m=sizes)
